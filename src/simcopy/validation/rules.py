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
        if rf.gamma == 0.0:              # sem Kerr: comprimento nao linear infinito
            continue
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

MIN_MEASURED_SYMBOLS = 100

def edge_guard_fits_window(system, simulation, res):
    """A medida de um canal IMDD descarta edge_guard_symbols simbolos em cada borda da janela,
    contaminados pela emenda da janela periodica (config/resolve.py; docs/sistema_imdd.md,
    secao 3). Sem simbolos sobrando nao ha medida (erro). Com menos de MIN_MEASURED_SYMBOLS,
    cerca de 50 por nivel, o desvio de cada nivel sai com erro relativo acima de 1/sqrt(2*50) =
    10 % (aviso)."""
    v = []
    for c, rc in zip(system.channels, res.channels):
        if c.kind != "imdd":
            continue
        g, n = rc.edge_guard_symbols, rc.n_symbols
        measured = n - 2 * g
        if measured < 2:
            v.append(Violation("EdgeGuardFitsWindow", Severity.ERROR,
                f"Canal {c.index}: a guarda de {g} simbolos em cada borda (dispersao acumulada "
                f"do enlace e memoria do transmissor) nao deixa simbolos para medir numa janela "
                f"de {n}.",
                f"n_symbols >= {2 * g + MIN_MEASURED_SYMBOLS}."))
        elif measured < MIN_MEASURED_SYMBOLS:
            v.append(Violation("EdgeGuardFitsWindow", Severity.WARNING,
                f"Canal {c.index}: so {measured} dos {n} simbolos entram na medida (guarda de {g} "
                f"em cada borda); o fator Q sai com erro estatistico acima de ~10 %.",
                f"n_symbols >= {2 * g + MIN_MEASURED_SYMBOLS}."))
    return v

def channel_seeds_are_distinct(system, simulation, res):
    """Canais com a mesma semente transmitem os mesmos bits: era o defeito do rascunho, em que os
    4 canais WDM carregavam a mesma sequencia (correcao C-1)."""
    seen = {}
    v = []
    for c in system.channels:
        if c.seed in seen:
            v.append(Violation("ChannelSeedsAreDistinct", Severity.WARNING,
                f"Canais {seen[c.seed]} e {c.index} tem a mesma semente ({c.seed}): transmitem "
                f"os mesmos bits.",
                "Uma semente por canal."))
        seen.setdefault(c.seed, c.index)
    return v

def coherent_needs_carrier_recovery(system, simulation, res):
    """O receptor coerente ainda nao tem recuperacao de portadora (docs/sistema_coerente.md,
    secao 6): a medida so tira uma rotacao constante da constelacao. Com oscilador local livre
    (ruido de fase proprio) ou com desvio de frequencia, a constelacao gira e a BER medida nao
    representa o enlace. Aviso, nao erro: o sinal continua valido para ver o batimento."""
    v = []
    for c in system.channels:
        if c.kind != "dcs" or c.rx.dsp is not None:      # com rx.dsp a portadora e recuperada
            continue
        causes = []
        if c.rx.lo_source == "free":
            causes.append("oscilador local livre (rx.lo_source = 'free')")
        if c.rx.lo_frequency_offset != 0.0:
            causes.append(f"desvio de frequencia de {c.rx.lo_frequency_offset/1e6:g} MHz")
        if causes:
            v.append(Violation("CoherentNeedsCarrierRecovery", Severity.WARNING,
                f"Canal {c.index}: {' e '.join(causes)}, e o receptor ainda nao recupera a "
                f"portadora. A constelacao gira e a BER medida nao representa o enlace.",
                "Para o homodino ideal: rx.lo_source = 'tx_laser' e rx.lo_frequency_offset = 0."))
    return v

def coherent_mzm_at_null(system, simulation, res):
    """No coerente cada MZM do modulador em quadratura trabalha no nulo, V_b = V_pi (Agrawal,
    secao 10.1.2; o professor usa Vpi = Vbias no `DPQAM_transmitter`). `vbias` continua sendo um
    parametro com esse significado, e esta regra exige vbias = vpi (duvida 10, opcao a, resposta
    de 2026-11-16). O padrao de MZMConfig, vbias = vpi/2, e a quadratura do IMDD."""
    v = []
    for c in system.channels:
        if c.kind != "dcs":
            continue
        m = c.tx.mzm
        if abs(m.vbias - m.vpi) > 1e-9 * abs(m.vpi):
            v.append(Violation("CoherentMZMAtNull", Severity.ERROR,
                f"Canal {c.index}: tx.mzm.vbias = {m.vbias:g} V e tx.mzm.vpi = {m.vpi:g} V. No "
                f"coerente cada MZM trabalha no nulo, em que o campo e bipolar.",
                f"tx.mzm.vbias = {m.vpi:g} (igual a vpi)."))
    return v

SBS_THRESHOLD = 5e-3          # W
SBS_LONG_FIBER = 10e3         # m

def imdd_power_below_brillouin(system, simulation, res):
    """Nota da aula "IMDD < 5 mW" (duvida 42, resposta de 2026-11-16: e potencia, aviso na
    montagem). Agrawal, secoes 2.6.1 e 6.3: o limiar do espalhamento Brillouin estimulado e da
    ordem de 5 mW por canal em fibras longas (> 10 km), e sobe para perto de 50 mW em 1 km. O
    simulador nao tem Brillouin: acima do limiar a potencia que chega ao fim da fibra seria menor
    do que a simulada. So confere as potencias que a configuracao impoe (launch_power do canal e
    output_power dos amplificadores, por canal)."""
    if not any(c.kind == "imdd" for c in system.channels) or not system.link.spans:
        return []
    v = []
    def check(power, where, k):
        fiber = system.link.spans[k].fiber
        if power is not None and power > SBS_THRESHOLD and fiber.length > SBS_LONG_FIBER:
            v.append(Violation("IMDDPowerBelowBrillouin", Severity.WARNING,
                f"{where} = {1e3 * power:.3g} mW por canal entra em {fiber.length / 1e3:g} km de "
                f"fibra, acima do limiar de Brillouin (~5 mW em fibras de mais de 10 km). O "
                f"simulador nao tem esse efeito.",
                "Agrawal, secoes 2.6.1 e 6.3. Abaixo de 5 mW por canal o aviso some."))
    booster = system.link.booster
    if booster is not None and booster.mode == "target_output_power":
        check(booster.output_power, "link.booster.output_power", 0)
    elif booster is None:
        for c in system.channels:
            if c.kind == "imdd":
                check(c.launch_power, f"channels[{c.index}].launch_power", 0)
    for k, span in enumerate(system.link.spans[:-1]):
        a = span.amplifier
        if a is not None and a.mode == "target_output_power":
            check(a.output_power, f"link.spans[{k}].amplifier.output_power", k + 1)
    return v

ALL = (nyquist_satisfied, integer_samples_per_symbol, channels_do_not_overlap,
       fft_size_is_efficient, ssfm_step_resolves_nonlinearity, window_fits_symbols,
       dsp_rate_divides_simulation_rate, edge_guard_fits_window, channel_seeds_are_distinct,
       coherent_needs_carrier_recovery, coherent_mzm_at_null, imdd_power_below_brillouin)

def validate(system: SystemConfig, simulation: SimulationConfig, raise_on_error=False):
    res = resolve(system, simulation)
    out = []
    for rule in ALL:
        out.extend(rule(system, simulation, res))
    if raise_on_error and any(x.severity is Severity.ERROR for x in out):
        raise ValidationFailed(out)
    return out
