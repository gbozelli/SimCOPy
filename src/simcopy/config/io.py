"""Carregamento de configuracao em YAML."""
from __future__ import annotations
from typing import Any
import yaml
from ..core.units import C_LIGHT
from .system import (SystemConfig, LinkConfig, SpanConfig, FiberConfig,
                     BirefringenceConfig, AmplifierConfig, CoherentChannel,
                     IMDDChannel, CoherentTransmitter, CoherentReceiver,
                     IMDDTransmitter, IMDDReceiver, DirectModulation,
                     ExternalModulation, LaserConfig, MZMConfig, IQImpairments,
                     PhotodiodeConfig, ADCConfig, Pol, DCS_FORMATS, IMDD_FORMATS,
                     ConfigError)
from .simulation import SimulationConfig

_INT_KEYS = {"index", "seed", "s21_order", "filter_order", "samples_per_symbol",
             "resolution_bits", "n_symbols", "root_seed"}

# variáveis com _ n são chamadas externamente, 
# mas preciso chamá-las para testes
def _num(x: Any, key=None) -> float:
    """
    PyYAML (YAML 1.1) nao reconhece 59.84e9 como float, so 59.84e+9.
    Então precisamos tratá-lo aqui
    """
    if isinstance(x, dict):
        return {k: _num(v, k) for k, v in x.items()}
    if isinstance(x, list):
        return [_num(v, key) for v in x]
    if isinstance(x, str):
        try:
            x = float(x)
        except ValueError:
            return x
    if key in _INT_KEYS and isinstance(x, float):
        return int(x)
    return x

# o q são essas funções? Elas são chamadas externamente?
# vc só colocou pq sim?
def _laser(d): return LaserConfig(**(d or {}))
def _mzm(d):   return MZMConfig(**(d or {}))
def _pd(d):    return PhotodiodeConfig(**(d or {}))
def _adc(d):   return ADCConfig(**(d or {}))

# já podemos extrair channel daqui?
# n seria melhor fazer isso em outro lugar? tipo em system.py?

def _coherent(channel: dict) -> CoherentChannel: 
    # precisamos nomear essas variáveis de forma mais clara, esse era um problema do código antigo
    """
    Carrega as configurações d eum canal corente digital
    Sistemas coerentes já carregam pré-definicções por si só,
    então podemos carregar parâmetros de tx e rx.

    
    """
    tx, rx = channel.get("tx") or {}, channel.get("rx") or {}
    return CoherentChannel(
        index=channel["index"], symbol_rate=channel["symbol_rate"], seed=channel["seed"],
        # Obrigatoriamnete a polarização precisa ser dupla
        modulation=DCS_FORMATS[channel["modulation"]], pol=Pol[channel.get("pol", "DP")],

        tx=CoherentTransmitter(laser=_laser(tx.get("laser")), mzm=_mzm(tx.get("mzm")),
                               iq=IQImpairments(**(tx.get("iq") or {})),
                               # porque estamos definindo rolloff aqui? Ela não é derivada de YAML?
                               rolloff=tx.get("rolloff", 0.18)),

        rx=CoherentReceiver(lo=_laser(rx.get("lo")),
                            # o mesmo aqui para frequência e banda óptica
                            lo_frequency_offset=rx.get("lo_frequency_offset", 500e6),
                            photodiode=_pd(rx.get("photodiode")), adc=_adc(rx.get("adc")),
                            optical_bandwidth=rx.get("optical_bandwidth", 70e9)))

def _imdd(channel: dict) -> IMDDChannel:
    """
    Carrega as configurações de um canal IMDD.
    Para IMDD, pré-definições também podem ser carregadas,
    como por exemplo, a modulação.
    """
    tx, rx = channel.get("tx") or {}, channel.get("rx") or {}
    # fe significa o que exatamente?
    fe = tx.get("frontend") or {"scheme": "direct"} # o q sisgnifica schema mesmo?
    if fe.get("scheme") == "external":
        # já estamos mexendo com frontend óptico? daria pra simular isso então
        frontend = ExternalModulation(laser=_laser(fe.get("laser")), mzm=_mzm(fe.get("mzm")))
    else:
        frontend = DirectModulation(laser=_laser(fe.get("laser")),
                                    chirp_alpha=fe.get("chirp_alpha", 0.0))
    if "pol" in channel and channel["pol"] != "SP":
        raise ConfigError(f"Canal {c['index']}: 'pol' nao se aplica.\n"
                          f"  Causa: channel.kind = IMDD, que e sempre SP ([D-15]).")
    return IMDDChannel(index=c["index"], symbol_rate=c["symbol_rate"], seed=c["seed"],
                       modulation=IMDD_FORMATS[c["modulation"]],
                       tx=IMDDTransmitter(frontend=frontend, rolloff=tx.get("rolloff", 0.2)),
                       rx=IMDDReceiver(photodiode=_pd(rx.get("photodiode")),
                                       adc=_adc(rx.get("adc"))))

def _center_frequency(s: dict) -> float:
    """
    Aceita center_frequency OU center_wavelength; os dois so se forem coerentes.
    Ambos precisam respeitar a relação v = lambda * f
    """
    f, lam = s.get("center_frequency"), s.get("center_wavelength")
    if f is None and lam is None:
        return 193.7e12
    if f is not None and lam is not None:
        if abs(C_LIGHT / f - lam) / lam > 1e-6:
            raise ConfigError(
                f"center_frequency = {f/1e12:.4f} THz e center_wavelength = "
                f"{lam*1e9:.3f} nm nao sao coerentes: c/f = {C_LIGHT/f*1e9:.3f} nm.\n"
                f"  Informe apenas um dos dois; o outro e derivado.")
        return f
    return f if f is not None else C_LIGHT / lam

def load(path: str) -> tuple[SystemConfig, SimulationConfig]:
    with open(path) as fh:
        raw = _num(yaml.safe_load(fh))
    s, sim = raw["system"], raw["simulation"]

    chans = []
    for c in s["channels"]:
        if c["kind"] == "dcs":
            chans.append(_coherent(c))
        elif c["kind"] == "imdd":
            chans.append(_imdd(c))
        else:
            raise ConfigError(f"kind desconhecido: {c['kind']!r} (use 'imdd' ou 'dcs')")

    lk = s["link"]
    spans = []
    for sp in lk["spans"]:
        f = dict(sp["fiber"])
        b = f.pop("birefringence", {}) or {}
        bir = BirefringenceConfig(
            pmd_coefficient=b.get("pmd_coefficient_ps_sqrt_km", 0.1) * 1e-12 / 31.62,
            correlation_length=b.get("correlation_length", 50.0))
        spans.append(SpanConfig(fiber=FiberConfig(birefringence=bir, **f),
                                line_compensation=sp.get("line_compensation", False)))
    b = lk.get("booster")
    booster = None
    if b:
        booster = AmplifierConfig(
            mode=b.get("mode", "target_output_power"),
            output_power=10 ** ((b["output_power_dbm"] - 30) / 10)
                         if "output_power_dbm" in b else None,
            noise_figure_db=b.get("noise_figure_db", 5.0))

    system = SystemConfig(channels=tuple(chans),
                          link=LinkConfig(spans=tuple(spans), booster=booster),
                          center_frequency=_center_frequency(s),
                          channel_spacing=s.get("channel_spacing", 75e9))
    simulation = SimulationConfig.from_symbols(
        symbol_rate=chans[0].symbol_rate,
        samples_per_symbol=sim["samples_per_symbol"], n_symbols=sim["n_symbols"],
        root_seed=sim.get("root_seed", 0),
        force_power_of_two=sim.get("force_power_of_two", True),
        ssfm_max_step=sim.get("ssfm_max_step", 5e3),
        ssfm_max_phase_change_deg=sim.get("ssfm_max_phase_change_deg", 0.5))
    return system, simulation
