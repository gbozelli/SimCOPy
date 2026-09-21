"""Carregamento de configuracao em YAML.

Estrutura: transmissor e receptor ficam DENTRO de cada canal, abaixo da decisao
de `kind` (IMDD ou DCS). Para nao repetir blocos iguais em N canais, use ancoras
YAML (&nome / *nome) e chaves de merge (<<: *nome) com sobrescrita local.
"""
from __future__ import annotations
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

def _num(x, key=None):
    """PyYAML (YAML 1.1) nao reconhece 59.84e9 como float, so 59.84e+9."""
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

def _laser(d): return LaserConfig(**(d or {}))
def _mzm(d):   return MZMConfig(**(d or {}))
def _pd(d):    return PhotodiodeConfig(**(d or {}))
def _adc(d):   return ADCConfig(**(d or {}))

def _coherent(c):
    tx, rx = c.get("tx") or {}, c.get("rx") or {}
    return CoherentChannel(
        index=c["index"], symbol_rate=c["symbol_rate"], seed=c["seed"],
        modulation=DCS_FORMATS[c["modulation"]], pol=Pol[c.get("pol", "DP")],
        tx=CoherentTransmitter(laser=_laser(tx.get("laser")), mzm=_mzm(tx.get("mzm")),
                               iq=IQImpairments(**(tx.get("iq") or {})),
                               rolloff=tx.get("rolloff", 0.18)),
        rx=CoherentReceiver(lo=_laser(rx.get("lo")),
                            lo_frequency_offset=rx.get("lo_frequency_offset", 500e6),
                            photodiode=_pd(rx.get("photodiode")), adc=_adc(rx.get("adc")),
                            optical_bandwidth=rx.get("optical_bandwidth", 70e9)))

def _imdd(c):
    tx, rx = c.get("tx") or {}, c.get("rx") or {}
    fe = tx.get("frontend") or {"scheme": "direct"}
    if fe.get("scheme") == "external":
        frontend = ExternalModulation(laser=_laser(fe.get("laser")), mzm=_mzm(fe.get("mzm")))
    else:
        frontend = DirectModulation(laser=_laser(fe.get("laser")),
                                    chirp_alpha=fe.get("chirp_alpha", 0.0))
    if "pol" in c and c["pol"] != "SP":
        raise ConfigError(f"Canal {c['index']}: 'pol' nao se aplica.\n"
                          f"  Causa: channel.kind = IMDD, que e sempre SP ([D-15]).")
    return IMDDChannel(index=c["index"], symbol_rate=c["symbol_rate"], seed=c["seed"],
                       modulation=IMDD_FORMATS[c["modulation"]],
                       tx=IMDDTransmitter(frontend=frontend, rolloff=tx.get("rolloff", 0.2)),
                       rx=IMDDReceiver(photodiode=_pd(rx.get("photodiode")),
                                       adc=_adc(rx.get("adc"))))

def _center_frequency(s):
    """Aceita center_frequency OU center_wavelength; os dois so se forem coerentes."""
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

def load(path):
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
