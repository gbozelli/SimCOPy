"""Montagem do sistema IMDD a partir da configuracao (docs/sistema_imdd.md, secao 4).

A montagem devolve uma lista de (id, bloco, entradas), na ordem de execucao. `entradas` sao os ids
dos estagios cujas saidas o bloco recebe: nenhuma na fonte de bits, uma na maioria, uma por canal no
mux. Com um canal: transmissor -> enlace -> receptor. Com N canais: N transmissores -> mux ->
enlace (campo agregado) -> N receptores, cada um comecando pelo demux do seu canal.

As condicionais da cascata (modulacao direta ou externa, potencia de lancamento, amplificadores,
WDM) sao escolhas da montagem, nunca `if` dentro de um bloco (arquitetura, secao 6.2). Os ids sao a
chave das sementes dos blocos novos e o rotulo de cada estagio no resultado.

Ainda nao entram, e a montagem recusa com ConfigError: PAM4, compensacao de linha e o canal
coerente.
"""
from __future__ import annotations
from dataclasses import replace
from ..components.optical.amplifier import Amplifier
from ..components.optical.fiber import Fiber
from ..components.optical.power import PowerControl
from ..components.optical.wdm import Demultiplexer, Multiplexer
from ..components.rx.photodiode import Photodiode
from ..components.tx.bits import BitSource
from ..components.tx.driver import NRZDriver
from ..components.tx.laser import DirectModulatedLaser
from ..components.tx.modulator import ExternalModulator
from ..config.resolve import Resolved, ResolvedChannel
from ..config.simulation import SimulationConfig
from ..config.system import ConfigError, DirectModulation, IMDDChannel, SystemConfig

Chain = list[tuple[str, object]]                      # estagios em serie: (id, bloco)
Pipeline = list[tuple[str, object, tuple[str, ...]]]  # (id, bloco, ids das entradas)

def check_supported(system: SystemConfig) -> None:
    """Recusa, com a causa, o que a montagem ainda nao cobre."""
    for c in system.channels:
        if c.kind != "imdd":
            raise ConfigError(f"Canal {c.index}: kind = {c.kind!r}. A montagem cobre so IMDD por "
                              f"enquanto; o receptor coerente ainda nao existe.")
        if c.modulation.name != "OOK":
            raise ConfigError(f"Canal {c.index}: modulation = {c.modulation.name}. So OOK por "
                              f"enquanto (duvida 51).")
        if system.n_channels > 1 and c.rx.optical_filter is None:
            raise ConfigError(f"Canal {c.index}: com {system.n_channels} canais o receptor precisa "
                              f"do filtro do demux (rx.optical_filter: bandwidth, filter_type).")
    for k, span in enumerate(system.link.spans):
        if span.line_compensation:
            raise ConfigError(f"link.spans[{k}].line_compensation = true: compensacao de linha "
                              f"ainda nao implementada (duvida 52). A fibra compensadora entra "
                              f"como mais um elemento de link.spans.")

def channel_frequency(system: SystemConfig, rc: ResolvedChannel) -> float:
    """Frequencia absoluta do canal: centro da banda + desvio do canal."""
    return system.center_frequency + rc.f0

def build_imdd_tx(channel: IMDDChannel, rc: ResolvedChannel, system: SystemConfig) -> Chain:
    """Bits -> sinal eletrico -> campo optico, em banda base do canal, com o desvio do canal em
    Waveform.f0 (quem desloca o canal para a grade e o mux)."""
    p = f"ch{channel.index}.tx"
    freq = channel_frequency(system, rc)          # h*nu do laser
    chain: Chain = [
        (f"{p}.bits", BitSource(rc.geometry, channel.seed, id=f"{p}.bits")),
        (f"{p}.driver", NRZDriver(channel.symbol_rate, channel.tx.driver, channel.index,
                                  id=f"{p}.driver")),
    ]
    fe = channel.tx.frontend
    if isinstance(fe, DirectModulation):
        chain.append((f"{p}.laser", DirectModulatedLaser(fe.laser, freq, fe.bias_ratio,
                                                         fe.modulation_ratio,
                                                         channel_offset=rc.f0)))
    else:
        chain.append((f"{p}.modulator", ExternalModulator(fe.laser, fe.mzm, freq,
                                                          fe.drive_amplitude,
                                                          id=f"{p}.modulator",
                                                          channel_offset=rc.f0)))
    if channel.launch_power is not None:
        chain.append((f"{p}.power", PowerControl(channel.launch_power, id=f"{p}.power")))
    return chain

def build_link(system: SystemConfig, simulation: SimulationConfig) -> Chain:
    """Booster (opcional) e, para cada elemento de link.spans, a fibra e o amplificador (opcional).
    Fibra padrao + compensadora + amplificador sao dois elementos: [SMF] e [DCF + amplificador]."""
    fc, n = system.center_frequency, system.n_channels
    chain: Chain = []
    if system.link.booster is not None:
        chain.append(("link.booster", Amplifier(system.link.booster, fc, n, id="link.booster")))
    for k, span in enumerate(system.link.spans):
        s = f"link.span{k}"
        chain.append((f"{s}.fiber", Fiber(span.fiber, fc, id=f"{s}.fiber",
                                          max_step=simulation.ssfm_max_step,
                                          max_phase_deg=simulation.ssfm_max_phase_change_deg)))
        if span.amplifier is not None:
            chain.append((f"{s}.amplifier", Amplifier(span.amplifier, fc, n,
                                                      id=f"{s}.amplifier")))
    return chain

def build_imdd_rx(channel: IMDDChannel, rc: ResolvedChannel,
                  simulation: SimulationConfig) -> Chain:
    """Demux (se configurado) e fotodetector. As chaves de ruido tem dois donos (duvida 6): o ruido
    entra se o dispositivo e a simulacao o pedem."""
    p = f"ch{channel.index}.rx"
    chain: Chain = []
    if channel.rx.optical_filter is not None:
        chain.append((f"{p}.demux", Demultiplexer(channel.rx.optical_filter, rc.f0,
                                                  channel.index, rc.geometry, id=f"{p}.demux")))
    pd = channel.rx.photodiode
    pd = replace(pd, shot_noise=pd.shot_noise and simulation.shot_noise,
                 thermal_noise_density=pd.thermal_noise_density if simulation.thermal_noise
                 else 0.0)
    chain.append((f"{p}.photodiode", Photodiode(pd, channel.symbol_rate,
                                                id=f"{p}.photodiode")))
    return chain

def _in_series(chain: Chain, first_input: str | None) -> Pipeline:
    out, prev = [], first_input
    for sid, block in chain:
        out.append((sid, block, (prev,) if prev is not None else ()))
        prev = sid
    return out

def build(system: SystemConfig, simulation: SimulationConfig, resolved: Resolved) -> Pipeline:
    """Transmissores -> (mux) -> enlace -> receptores, na ordem de execucao."""
    check_supported(system)
    pipeline: Pipeline = []
    tx_ends = []
    for c, rc in zip(system.channels, resolved.channels):
        tx = _in_series(build_imdd_tx(c, rc, system), None)
        pipeline += tx
        tx_ends.append(tx[-1][0])
    if system.n_channels > 1:
        offsets = tuple(rc.f0 for rc in resolved.channels)
        pipeline.append(("mux", Multiplexer(offsets, system.mux, id="mux"), tuple(tx_ends)))
        link_input = "mux"
    else:
        link_input = tx_ends[0]
    link = _in_series(build_link(system, simulation), link_input)
    pipeline += link
    link_end = link[-1][0] if link else link_input
    for c, rc in zip(system.channels, resolved.channels):
        pipeline += _in_series(build_imdd_rx(c, rc, simulation), link_end)
    return pipeline
