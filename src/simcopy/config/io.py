"""Carregamento de configuracao em YAML ou JSON (o JSON e escolhido pela extensao .json)."""
from __future__ import annotations
import json
from typing import Any
import yaml
from ..core.units import C_LIGHT, dbm_to_w
from .system import (SystemConfig, LinkConfig, SpanConfig, FiberConfig,
                     BirefringenceConfig, AmplifierConfig, CoherentChannel,
                     IMDDChannel, CoherentTransmitter, CoherentReceiver,
                     IMDDTransmitter, IMDDReceiver, DirectModulation,
                     ExternalModulation, CWLaserConfig, MZMConfig, IQImpairments,
                     PhotodiodeConfig, ADCConfig, RateEquationLaserConfig,
                     ElectricalFilterConfig, OpticalFilterConfig, MuxConfig, Pol,
                     DCS_FORMATS, IMDD_FORMATS,
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
def _laser(d): return CWLaserConfig(**(d or {}))
def _mzm(d):   return MZMConfig(**(d or {}))
def _pd(d):    return PhotodiodeConfig(**(d or {}))
def _adc(d):   return ADCConfig(**(d or {}))
def _driver(d): return ElectricalFilterConfig(**d) if d else None

def _optical_filter(d, where: str):
    if not d:
        return None
    _only(d, {"bandwidth", "filter_type", "filter_order"}, where)
    return OpticalFilterConfig(**d)

def _only(d: dict, allowed: set, where: str) -> dict:
    """Recusa chave desconhecida em vez de descarta-la em silencio (issue #9)."""
    extra = sorted(set(d) - allowed)
    if extra:
        raise ConfigError(f"{where}: chave(s) desconhecida(s) {extra}.\n"
                          f"  Chaves aceitas: {sorted(allowed)}.")
    return d

def _required(d: dict, key: str, where: str):
    if key not in d:
        raise ConfigError(f"{where}: falta '{key}' (campo obrigatorio, sem valor padrao).")
    return d[key]

def _power(d: dict, key: str):
    """Potencia em W (`key`) ou em dBm (`key_dbm`); nunca as duas."""
    if key in d and f"{key}_dbm" in d:
        raise ConfigError(f"Informe '{key}' (W) ou '{key}_dbm', nao os dois.")
    if f"{key}_dbm" in d:
        return dbm_to_w(d[f"{key}_dbm"])
    return d.get(key)

def _rate_laser(d: dict, where: str) -> RateEquationLaserConfig:
    try:
        return RateEquationLaserConfig(**d)
    except TypeError as e:          # campo faltando ou desconhecido: nomeia o lugar
        raise ConfigError(f"{where}: {e}") from None

def _amp(d: dict, where: str) -> AmplifierConfig:
    _only(d, {"mode", "gain_db", "output_power", "output_power_dbm", "noise_figure_db",
              "saturation_power", "saturation_power_dbm"}, where)
    kw = {k: d[k] for k in ("mode", "gain_db", "noise_figure_db") if k in d}
    for key in ("output_power", "saturation_power"):
        if _power(d, key) is not None:
            kw[key] = _power(d, key)
    return AmplifierConfig(**kw)

def _mux(d) -> MuxConfig:
    if not d:
        return MuxConfig()
    _only(d, {"insertion_loss_db", "optical_filter"}, "system.mux")
    kw = {"insertion_loss_db": d["insertion_loss_db"]} if "insertion_loss_db" in d else {}
    return MuxConfig(optical_filter=_optical_filter(d.get("optical_filter"),
                                                    "system.mux.optical_filter"), **kw)

def _frontend(fe: dict, where: str):
    """Variante do transmissor IMDD. Sem `scheme`, vale o padrao do IMDDTransmitter (externa)."""
    scheme = fe.get("scheme", "external")
    if scheme == "external":
        _only(fe, {"scheme", "laser", "mzm", "drive_amplitude"}, where)
        kw = {"drive_amplitude": fe["drive_amplitude"]} if "drive_amplitude" in fe else {}
        return ExternalModulation(laser=_laser(fe.get("laser")), mzm=_mzm(fe.get("mzm")), **kw)
    if scheme == "direct":
        _only(fe, {"scheme", "laser", "bias_ratio", "modulation_ratio"}, where)
        return DirectModulation(
            laser=_rate_laser(_required(fe, "laser", where), f"{where}.laser"),
            bias_ratio=_required(fe, "bias_ratio", where),
            modulation_ratio=_required(fe, "modulation_ratio", where))
    raise ConfigError(f"{where}: scheme = {scheme!r}; use 'direct' ou 'external'.")

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
        # Obrigatoriamnete a polarização precisa ser dupla (NA VERDADE NÃO, ME CONFUNDI AQUI)
        modulation=DCS_FORMATS[channel["modulation"]], pol=Pol[channel.get("pol", "DP")],
        launch_power=_power(channel, "launch_power"),

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
    where = f"canal {channel.get('index')}: tx"
    _only(tx, {"frontend", "driver", "rolloff"}, where)
    # fe significa o que exatamente?
    # (frontend: a parte eletro-optica do transmissor, modulacao direta ou externa;
    #  scheme e o nome da variante, 'direct' ou 'external')
    fe = tx.get("frontend") or {} # o q sisgnifica schema mesmo?
    # já estamos mexendo com frontend óptico? daria pra simular isso então
    frontend = _frontend(fe, f"{where}.frontend")
    if "pol" in channel and channel["pol"] != "SP":
        raise ConfigError(f"Canal {channel['index']}: 'pol' nao se aplica.\n"
                          f"  Causa: channel.kind = IMDD, que e sempre SP ([D-15]).")
    kw = {"rolloff": tx["rolloff"]} if "rolloff" in tx else {}
    _only(rx, {"photodiode", "adc", "optical_filter"}, f"canal {channel.get('index')}: rx")
    return IMDDChannel(index=channel["index"], symbol_rate=channel["symbol_rate"],
                       seed=channel["seed"],
                       modulation=IMDD_FORMATS[channel["modulation"]],
                       launch_power=_power(channel, "launch_power"),
                       tx=IMDDTransmitter(frontend=frontend, driver=_driver(tx.get("driver")),
                                          **kw),
                       rx=IMDDReceiver(photodiode=_pd(rx.get("photodiode")),
                                       adc=_adc(rx.get("adc")),
                                       optical_filter=_optical_filter(
                                           rx.get("optical_filter"),
                                           f"canal {channel.get('index')}: rx.optical_filter")))

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

def _read(path) -> dict:
    with open(path, encoding="utf-8") as fh:
        if str(path).lower().endswith(".json"):
            return json.load(fh)
        return yaml.safe_load(fh)

def load(path: str) -> tuple[SystemConfig, SimulationConfig]:
    raw = _num(_read(path))
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
    for k, sp in enumerate(lk["spans"]):
        where = f"link.spans[{k}]"
        _only(sp, {"fiber", "amplifier", "line_compensation"}, where)
        f = dict(sp["fiber"])
        b = f.pop("birefringence", {}) or {}
        _only(b, {"pmd_coefficient_ps_sqrt_km", "correlation_length", "enabled"},
              f"{where}.fiber.birefringence")
        bir = BirefringenceConfig(
            pmd_coefficient=b.get("pmd_coefficient_ps_sqrt_km", 0.1) * 1e-12 / 31.62,
            correlation_length=b.get("correlation_length", 50.0),
            **({"enabled": b["enabled"]} if "enabled" in b else {}))
        amp = _amp(sp["amplifier"], f"{where}.amplifier") if sp.get("amplifier") else None
        spans.append(SpanConfig(fiber=FiberConfig(birefringence=bir, **f), amplifier=amp,
                                line_compensation=sp.get("line_compensation", False)))
    b = lk.get("booster")
    booster = _amp(b, "link.booster") if b else None

    system = SystemConfig(channels=tuple(chans),
                          link=LinkConfig(spans=tuple(spans), booster=booster),
                          center_frequency=_center_frequency(s),
                          channel_spacing=s.get("channel_spacing", 75e9),
                          mux=_mux(s.get("mux")))
    simulation = SimulationConfig.from_symbols(
        symbol_rate=chans[0].symbol_rate,
        samples_per_symbol=sim["samples_per_symbol"], n_symbols=sim["n_symbols"],
        root_seed=sim.get("root_seed", 0),
        force_power_of_two=sim.get("force_power_of_two", True),
        ssfm_max_step=sim.get("ssfm_max_step", 5e3),
        ssfm_max_phase_change_deg=sim.get("ssfm_max_phase_change_deg", 0.5),
        **{k: sim[k] for k in ("shot_noise", "thermal_noise", "broadening_margin") if k in sim})
    return system, simulation
