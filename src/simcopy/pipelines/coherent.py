"""Montagem do sistema coerente digital a partir da configuracao (docs/sistema_coerente.md, secoes 5
e 11).

Mesmo formato da montagem do IMDD (pipelines/imdd.py): lista de (id, bloco, entradas), na ordem de
execucao. Por canal:

    bits -> acionamento (Nyquist) -> laser CW -> modulador IQ -> (potencia)
    N canais -> mux -> enlace (booster, fibras, amplificadores)
    (filtro do canal) -> oscilador local -> hibrida de 90 graus + deteccao balanceada
    -> com rx.dsp: compensacao de dispersao -> (desvio de frequencia) -> (filtro casado)
    -> com rx.dsp.mode = "blind": ADC -> receptor cego

O laser CW e um estagio proprio: o modulador recebe a portadora e o sinal eletrico, e o oscilador
local travado (rx.lo_source = "tx_laser") recebe a mesma saida do laser. O laser so usa a entrada
para saber a grade e o canal. O oscilador local livre (rx.lo_source = "free") e outro laser CW,
com o papel de semente "lo_laser_phase".

Com mais de um canal, o mux leva cada canal para a sua posicao e cada receptor comeca pelo filtro
optico do seu canal. O oscilador local fica na frequencia do canal mais rx.lo_frequency_offset.

A parte do DSP que usa os simbolos transmitidos (instante, polarizacao, fase, equalizador) fica na
medida (kernels/demodulation.py). No modo cego a medida so casa a saida com os dados. Ainda nao
entram, e a montagem recusa com ConfigError: 32QAM, dupla polarizacao sem rx.dsp e compensacao de
linha. O ADC so existe no modo cego: no modo com simbolos conhecidos a reamostragem para 2
amostras por simbolo acontece na medida, sem quantizacao.
"""
from __future__ import annotations
from dataclasses import replace
from ..components.optical.power import PowerControl
from ..components.optical.wdm import Demultiplexer, Multiplexer
from ..components.rx.adc import ADC, BlindReceiver
from ..components.rx.coherent import CoherentFrontEnd, LockedLocalOscillator
from ..components.rx.dsp import DispersionCompensator, FrequencyOffsetCompensator, MatchedFilter
from ..components.tx.bits import BitSource
from ..components.tx.iq_driver import NyquistDriver
from ..components.tx.iq_modulator import IQModulator
from ..components.tx.laser import CWLaser
from ..config.resolve import Resolved, ResolvedChannel
from ..config.simulation import SimulationConfig
from ..config.system import (CoherentChannel, ConfigError, OpticalFilterConfig, Pol,
                             SystemConfig)
from ..kernels.mapping import BITS_PER_SYMBOL, SQUARE_QAM
from .imdd import Pipeline, _in_series, build_link, channel_frequency

LO_RNG_ROLE = "lo_laser_phase"

def measured_stage(channel: CoherentChannel) -> str:
    """Sufixo do id do estagio cuja saida a medida usa: o ultimo do receptor do canal."""
    dsp = channel.rx.dsp
    if dsp is None:
        return "rx.frontend"
    if dsp.mode == "blind":
        return "rx.dsp"
    if dsp.matched_filter and channel.tx.pulse == "root_raised_cosine":
        return "rx.filter"
    if dsp.frequency_offset:
        return "rx.frequency"
    return "rx.cd" if dsp.cd_compensation else "rx.frontend"

def check_supported(system: SystemConfig) -> None:
    """Recusa, com a causa, o que a montagem ainda nao cobre."""
    for c in system.channels:
        if c.kind != "dcs":
            raise ConfigError(f"Canal {c.index}: kind = {c.kind!r} numa montagem coerente.")
        if c.modulation.name not in BITS_PER_SYMBOL:
            raise ConfigError(f"Canal {c.index}: modulation = {c.modulation.name}. So "
                              f"{sorted(BITS_PER_SYMBOL)} por enquanto.")
        if c.rx.dsp is None:
            if c.pol is not Pol.SP:
                raise ConfigError(f"Canal {c.index}: pol = {c.pol.name} pede rx.dsp (quem separa as "
                                  f"polarizacoes e o DSP). Use `dsp: {{}}` para os padroes.")
            if c.modulation.name in SQUARE_QAM and c.tx.pulse == "root_raised_cosine":
                raise ConfigError(f"Canal {c.index}: tx.pulse = 'root_raised_cosine' pede o filtro "
                                  f"casado de rx.dsp.")
    for k, span in enumerate(system.link.spans):
        if span.line_compensation:
            raise ConfigError(f"link.spans[{k}].line_compensation = true: compensacao de linha "
                              f"ainda nao implementada (duvida 52).")

def build_coherent_tx(channel: CoherentChannel, rc: ResolvedChannel,
                      system: SystemConfig) -> Pipeline:
    """Bits -> I(t) + i*Q(t) -> campo optico, em banda base do canal (Waveform.f0 = desvio)."""
    p = f"ch{channel.index}.tx"
    tx = channel.tx
    pipeline: Pipeline = [
        (f"{p}.bits", BitSource(rc.geometry, channel.seed, id=f"{p}.bits"), ()),
        (f"{p}.driver", NyquistDriver(channel.symbol_rate, channel.modulation.name, tx.rolloff,
                                      channel.index, id=f"{p}.driver", pulse=tx.pulse),
         (f"{p}.bits",)),
        (f"{p}.laser", CWLaser(tx.laser, channel_frequency(system, rc), channel_offset=rc.f0),
         (f"{p}.driver",)),
        (f"{p}.modulator", IQModulator(tx.mzm, tx.iq, tx.drive_depth, tx.predistortion,
                                       id=f"{p}.modulator"), (f"{p}.laser", f"{p}.driver")),
    ]
    if channel.launch_power is not None:
        pipeline.append((f"{p}.power", PowerControl(channel.launch_power, id=f"{p}.power"),
                         (f"{p}.modulator",)))
    return pipeline

def build_coherent_rx(channel: CoherentChannel, rc: ResolvedChannel, system: SystemConfig,
                      simulation: SimulationConfig, resolved: Resolved, signal: str) -> Pipeline:
    """Filtro do canal (com mais de um canal), oscilador local, receptor e os blocos do DSP.
    `signal` e o id do estagio que entrega o sinal optico. As chaves de ruido tem dois donos, como
    no IMDD: o dispositivo e a simulacao."""
    p = f"ch{channel.index}.rx"
    rx = channel.rx
    if rx.dsp is not None and rx.dsp.mode == "blind":
        # recusado na montagem, antes de rodar o enlace: o equalizador cego trabalha a T/2
        if rx.adc.samples_per_symbol != 2:
            raise ConfigError(f"Canal {channel.index}: rx.dsp.mode = 'blind' pede "
                              f"rx.adc.samples_per_symbol = 2 (chegou {rx.adc.samples_per_symbol}).")
        if rc.samples_per_symbol % 2:
            raise ConfigError(f"Canal {channel.index}: {rc.samples_per_symbol} amostras por simbolo "
                              f"na grade nao dao as 2 do ADC; use um numero par.")
    pipeline: Pipeline = []
    if system.n_channels > 1:
        optical = OpticalFilterConfig(bandwidth=rx.optical_bandwidth, filter_type="supergaussian",
                                      filter_order=rx.optical_filter_order)
        pipeline.append((f"{p}.demux", Demultiplexer(optical, rc.f0, channel.index, rc.geometry,
                                                     id=f"{p}.demux"), (signal,)))
        signal = f"{p}.demux"
    # depois do mux o campo esta referenciado ao centro da banda: o oscilador local vai para a
    # frequencia do canal. Com um canal so, rc.f0 = 0.
    lo_frequency = rc.f0 + rx.lo_frequency_offset
    if rx.lo_source == "tx_laser":
        lo = (f"{p}.lo", LockedLocalOscillator(rx.lo.power, lo_frequency, id=f"{p}.lo"),
              (f"ch{channel.index}.tx.laser",))
    else:
        lo = (f"{p}.lo", CWLaser(rx.lo, channel_frequency(system, rc),
                                 frequency_offset=lo_frequency, rng_role=LO_RNG_ROLE,
                                 channel_offset=rc.f0), (signal,))
    pd = rx.photodiode
    pd = replace(pd, shot_noise=pd.shot_noise and simulation.shot_noise,
                 thermal_noise_density=pd.thermal_noise_density if simulation.thermal_noise
                 else 0.0)
    pipeline += [lo, (f"{p}.frontend", CoherentFrontEnd(pd, channel.symbol_rate, id=f"{p}.frontend",
                                                        n_pol=channel.pol.n_pol),
                      (signal, f"{p}.lo"))]
    dsp, last = rx.dsp, f"{p}.frontend"
    if dsp is None:
        return pipeline
    if dsp.cd_compensation:
        beta2_length = sum(f.beta2 * s.fiber.length
                           for f, s in zip(resolved.fibers, system.link.spans))
        beta3_length = sum(f.beta3 * s.fiber.length
                           for f, s in zip(resolved.fibers, system.link.spans))
        pipeline.append((f"{p}.cd", DispersionCompensator(beta2_length, lo_frequency,
                                                          id=f"{p}.cd",
                                                          beta3_length=beta3_length), (last,)))
        last = f"{p}.cd"
    if dsp.frequency_offset:
        pipeline.append((f"{p}.frequency", FrequencyOffsetCompensator(channel.symbol_rate,
                                                                      id=f"{p}.frequency"),
                         (last,)))
        last = f"{p}.frequency"
    if dsp.matched_filter and channel.tx.pulse == "root_raised_cosine":
        pipeline.append((f"{p}.filter", MatchedFilter(channel.symbol_rate, channel.tx.rolloff,
                                                      id=f"{p}.filter"), (last,)))
        last = f"{p}.filter"
    if dsp.mode == "blind":
        # DSP cego: ADC (2 amostras por simbolo, quantizado) e o receptor cego, que nao usa os
        # dados. A medida so casa a saida com os dados (kernels/demodulation.py, measure_blind).
        pipeline.append((f"{p}.adc", ADC(rx.adc, channel.symbol_rate, id=f"{p}.adc"), (last,)))
        pipeline.append((f"{p}.dsp", BlindReceiver(
            channel.modulation.name, taps=dsp.equalizer_taps, step_size=dsp.step_size,
            cma_iterations=dsp.cma_iterations, rde_iterations=dsp.rde_iterations,
            dd_iterations=dsp.dd_iterations, test_angles=dsp.test_angles,
            phase_window=dsp.phase_window, stokes_demux=dsp.stokes_demux, id=f"{p}.dsp"),
            (f"{p}.adc",)))
    return pipeline

def build(system: SystemConfig, simulation: SimulationConfig, resolved: Resolved) -> Pipeline:
    """Transmissores -> (mux) -> enlace -> receptores, na ordem de execucao."""
    check_supported(system)
    pipeline: Pipeline = []
    tx_ends = []
    for c, rc in zip(system.channels, resolved.channels):
        tx = build_coherent_tx(c, rc, system)
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
        pipeline += build_coherent_rx(c, rc, system, simulation, resolved, link_end)
    return pipeline
