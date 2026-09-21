"""Invariantes checados na construcao, antes de alocar qualquer array (PR-3).
Todos os validadores rodam e devolvem TODAS as violacoes, nao so a primeira."""
from __future__ import annotations
import numpy as np
from ..config.system import SystemConfig
from ..config.simulation import SimulationConfig
from ..config.resolve import resolve
from .base import Violation, Severity, ValidationFailed

def _tx_rolloff(c) -> float:
    return getattr(c.tx, "rolloff", 0.2)

def nyquist_satisfied(system, simulation, res):
    v = []
    edge = max(abs(f0) + c.symbol_rate * (1 + _tx_rolloff(c)) / 2
               for c, f0 in zip(system.channels, system.channel_offsets()))
    need = edge * simulation.broadening_margin
    half = simulation.sampling_rate / 2
    if need > half:
        v.append(Violation("RequiresNyquistSatisfied", Severity.ERROR,
            f"A grade WDM ocupa +-{edge/1e9:.1f} GHz e, com margem de "
            f"alargamento {simulation.broadening_margin:g}x, precisa de "
            f"+-{need/1e9:.1f} GHz, mas a banda de simulacao e +-{half/1e9:.1f} GHz.",
            f"fs >= {2*need/1e9:.1f} GHz, ou reduzir o espacamento de canais."))
    return v

def integer_samples_per_symbol(system, simulation, res):
    v = []
    for c in system.channels:
        r = simulation.sampling_rate / c.symbol_rate
        if abs(r - round(r)) > 1e-9:
            v.append(Violation("IntegerSamplesPerSymbol", Severity.ERROR,
                f"Canal {c.index}: fs/Rs = {r:.6f} nao e inteiro.",
                f"[D-4] exige fs multiplo inteiro de todos os baud rates. "
                f"fs = {c.symbol_rate*round(r)/1e9:.2f} GHz resolveria."))
    return v

def channels_do_not_overlap(system, simulation, res):
    v = []
    offs = sorted(zip(system.channel_offsets(), system.channels), key=lambda x: x[0])
    for (f1, c1), (f2, c2) in zip(offs, offs[1:]):
        need = (c1.symbol_rate*(1+_tx_rolloff(c1)) + c2.symbol_rate*(1+_tx_rolloff(c2)))/2
        if f2 - f1 < need:
            v.append(Violation("ChannelsDoNotOverlap", Severity.WARNING,
                f"Canais {c1.index} e {c2.index}: espacamento {(f2-f1)/1e9:.1f} GHz "
                f"menor que a banda ocupada somada ({need/1e9:.1f} GHz).",
                "Aumentar channel_spacing ou reduzir o roll-off."))
    return v

def fft_size_is_efficient(system, simulation, res):
    n = res.n_samples
    if n & (n - 1):
        return [Violation("FFTSizeIsEfficient", Severity.WARNING,
            f"N = {n} nao e potencia de 2.",
            "GreatestPrimeFactorLimit = 2 no VPI; ativar force_power_of_two.")]
    return []

def ssfm_step_resolves_nonlinearity(system, simulation, res):
    v = []
    p = max((c.launch_power or 0.0) for c in system.channels) or None
    if p is None:
        return v
    for i, (span, rf) in enumerate(zip(system.link.spans, res.fibers)):
        l_nl = 1.0 / (rf.gamma * p)
        if simulation.ssfm_max_step > 0.1 * l_nl:
            v.append(Violation("SSFMStepResolvesNonlinearity", Severity.WARNING,
                f"Span {i}: passo maximo {simulation.ssfm_max_step/1e3:.1f} km vs "
                f"comprimento nao linear {l_nl/1e3:.1f} km.",
                f"ssfm_max_step <= {0.1*l_nl/1e3:.2f} km."))
    return v

def window_fits_symbols(system, simulation, res):
    v = []
    for c, rc in zip(system.channels, res.channels):
        if rc.n_symbols < 1024:
            v.append(Violation("WindowFitsSymbols", Severity.WARNING,
                f"Canal {c.index}: apenas {rc.n_symbols} simbolos na janela.",
                "Para medir BER ~1e-3 com 10% de erro sao precisos ~3300 simbolos."))
    return v

def dsp_rate_divides_simulation_rate(system, simulation, res):
    """A saida do ADC (taxa do DSP) tem que ser obtida da grade por dizimacao
    inteira. Com SpS da grade = 16 e do DSP = 2, o fator e 8. Se o SpS do DSP nao
    divide o da grade, o ADC precisaria de reamostragem fracionaria, que [D-4]
    deixou fora da v1."""
    v = []
    for c, rc in zip(system.channels, res.channels):
        dsp = c.rx.adc.samples_per_symbol
        sim = rc.samples_per_symbol
        if sim % dsp:
            v.append(Violation("DSPRateDividesSimulationRate", Severity.ERROR,
                f"Canal {c.index}: o DSP opera a {dsp} amostras/simbolo e a grade a "
                f"{sim}; {sim}/{dsp} nao e inteiro.",
                f"Use um SpS do DSP que divida {sim} (por exemplo "
                f"{', '.join(str(d) for d in range(1, sim + 1) if sim % d == 0 and d <= 8)})."))
    return v

ALL = (nyquist_satisfied, integer_samples_per_symbol, channels_do_not_overlap,
       fft_size_is_efficient, ssfm_step_resolves_nonlinearity, window_fits_symbols,
       dsp_rate_divides_simulation_rate)

def validate(system: SystemConfig, simulation: SimulationConfig, raise_on_error=False):
    res = resolve(system, simulation)
    out = []
    for rule in ALL:
        out.extend(rule(system, simulation, res))
    if raise_on_error and any(x.severity is Severity.ERROR for x in out):
        raise ValidationFailed(out)
    return out
