"""Medida do receptor coerente: sincronismo, EVM, BER contada, fator Q, MI e GMI
(docs/sistema_coerente.md, secoes 4, 11 e 12).

Sincronismo (decisao de 2026-11-14): o simulador conhece os simbolos transmitidos, entao o instante
dos simbolos sai da correlacao circular entre o sinal recebido e o trem de simbolos transmitido,
procurada na janela inteira. O codigo do professor (draft.py, `QAM_receiver`) correlaciona com uma
sequencia de sincronismo BPSK no inicio do quadro; essa alternativa esta descrita no guia.

Escala e fase: um ganho complexo g, de minimos quadrados entre os simbolos recebidos r e os
transmitidos a, g = sum(r*conj(a))/sum(|a|^2), tira a escala da corrente e uma rotacao constante
da constelacao. Nao corrige fase que varia no tempo: sem recuperacao de portadora, um oscilador
local com desvio de frequencia ou ruido de fase proprio aparece como EVM e erros.

EVM (rms, relativa a potencia media da constelacao): sqrt(sum|r/g - a|^2 / sum|a|^2). Com ruido
aditivo, SNR = 1/EVM^2.

Fator Q (Agrawal, Eq. 10.3.2 e 10.3.7): no BPSK e no QPSK cada bit e uma decisao binaria numa
quadratura, com niveis I1 = -I0. Juntando as decisoes de todos os bits,
    Q = (I1 - I0)/(s1 + s0),   BER = (1/2)*erfc(Q/sqrt(2)),
a mesma conta de kernels/decision.py. O limiar e zero (niveis simetricos).

QAM quadrado (16QAM, 64QAM): os bits nao sao decisoes binarias, entao o fator Q acima nao se
aplica. `ber_q` passa a ser a BER que a SNR medida daria com ruido gaussiano e codigo de Gray,
    BER = (4/k)*(1 - 1/sqrt(M))*(1/2)*erfc{sqrt[3*SNR/(2*(M - 1))]},   k = log2(M),
que para o 16QAM e (3/8)*erfc(sqrt(SNR/10)), a mesma do ajuste do pacote
(reference/vpi/ajuste_abc.py). E a aproximacao usual de vizinho mais proximo, nao esta no Agrawal.
`q` e o fator Q equivalente, o que daria essa BER pela Eq. 10.3.2.

MI e GMI (kernels/information.py): em bits por simbolo e por polarizacao, com a variancia do ruido
medida nos proprios simbolos.

Tres caminhos:
- sem rx.dsp (`measure_coherent`): uma polarizacao, so o instante dos simbolos e um ganho complexo;
- rx.dsp com os simbolos conhecidos (`measure_data_aided`): a cadeia de kernels/dsp.py;
- rx.dsp cego (`measure_blind`): os simbolos ja saem do receptor cego. A medida so os casa com os
  dados (polarizacao, atraso e giro), como o `SynchDSP` do VPI, e conta os saltos de fase.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import erfc, erfcinv
from .decision import ber_from_q, level_statistics, q_factor
from .dsp import data_aided_receiver
from .information import information_rates
from .mapping import (BINARY_PER_QUADRATURE, BITS_PER_SYMBOL, ROTATIONAL_SYMMETRY, constellation,
                      decision_components, demap, map_bits, nearest_points)

@dataclass(frozen=True)
class CoherentMeasurement:
    q: float                      # fator Q das decisoes binarias (Eq. 10.3.7), ou o equivalente
    ber_q: float                  # (1/2)*erfc(Q/sqrt(2)), Eq. 10.3.2
    ber_counted: float            # erros / bits
    errors: int
    n_bits: int                   # bits usados na medida (sem as guardas)
    evm: float                    # EVM rms (1 = 100 %)
    snr_db: float                 # -20*log10(EVM)
    gain: float                   # |g|: amplitude recebida de um nivel unitario [unidade do sinal]
    phase: float                  # arg(g) [rad]: rotacao constante da constelacao
    delay: int                    # amostras do centro do simbolo transmitido ao instante de decisao
    first_symbol: int             # primeiro simbolo medido
    last_symbol: int              # ultimo simbolo medido + 1
    ser: float = 0.0              # erros de simbolo / simbolos
    symbol_errors: int = 0
    n_symbols: int = 0            # simbolos usados na medida, somando as polarizacoes
    snr_db_pol: tuple = ()        # por polarizacao (X, Y)
    ber_pol: tuple = ()
    ser_pol: tuple = ()
    mi: float = 0.0               # informacao mutua [bits por simbolo, por polarizacao]
    gmi: float = 0.0              # informacao mutua generalizada [bits por simbolo, por polarizacao]
    mi_pol: tuple = ()
    gmi_pol: tuple = ()
    phase_slips: int = 0          # saltos de fase do receptor cego, contados contra os dados

def ber_square_qam(snr, modulation: str):
    """BER de um QAM quadrado com codigo de Gray e ruido gaussiano, para a SNR linear `snr`
    (Es/N0): aproximacao de vizinho mais proximo. No 16QAM, (3/8)*erfc(sqrt(snr/10))."""
    k = BITS_PER_SYMBOL[modulation]
    m = 2 ** k
    return (4.0 / k) * (1.0 - 1.0 / np.sqrt(m)) * 0.5 * erfc(
        np.sqrt(3.0 * np.asarray(snr, dtype=float) / (2.0 * (m - 1))))

def ber_16qam(snr):
    """BER do 16QAM com codigo de Gray e ruido gaussiano, para a SNR linear `snr` (Es/N0)."""
    return ber_square_qam(snr, "16QAM")

def symbol_metrics(received, bits, modulation: str, *, gain: float, phase: float, delay: int,
                   first_symbol: int, last_symbol: int, phase_slips: int = 0) -> CoherentMeasurement:
    """Metricas de simbolos ja normalizados (n_pol, n) contra os bits transmitidos (n_pol,
    n*bits_por_simbolo): EVM, erros de simbolo e de bit, MI e GMI, no total e por polarizacao."""
    y = np.atleast_2d(np.asarray(received, dtype=complex))
    b = np.atleast_2d(np.asarray(bits)).astype(np.uint8)
    a = np.array([map_bits(row, modulation) for row in b])
    noise = np.sum(np.abs(y - a) ** 2, axis=1)
    power = np.sum(np.abs(a) ** 2, axis=1)
    bit_errors = np.array([np.count_nonzero(demap(y[p], modulation) != b[p]) for p in range(y.shape[0])])
    sym_errors = np.array([np.count_nonzero(nearest_points(y[p], modulation) != a[p])
                           for p in range(y.shape[0])])
    evm = float(np.sqrt(noise.sum() / power.sum()))
    if modulation in BINARY_PER_QUADRATURE:
        parts = [decision_components(y[p], b[p], modulation) for p in range(y.shape[0])]
        values = np.concatenate([v for v, _ in parts])
        high = np.concatenate([h for _, h in parts])
        q = float("nan")
        if high.any() and not high.all():
            q = float(q_factor(*level_statistics(values, high)))
        ber_q = float(ber_from_q(q))
    else:
        ber_q = float(ber_square_qam(1.0 / evm ** 2, modulation)) if evm > 0 else 0.0
        q = float(np.sqrt(2.0) * erfcinv(2.0 * ber_q)) if ber_q > 0 else float("inf")
    with np.errstate(divide="ignore"):
        snr_pol = tuple(float(v) for v in -10.0 * np.log10(noise / power))
    points = constellation(modulation)
    rates = [information_rates(y[p], a[p], points) for p in range(y.shape[0])]
    return CoherentMeasurement(
        q=q, ber_q=ber_q, ber_counted=float(bit_errors.sum() / b.size), errors=int(bit_errors.sum()),
        n_bits=int(b.size), evm=evm,
        snr_db=float(-20.0 * np.log10(evm)) if evm > 0 else float("inf"),
        gain=float(gain), phase=float(phase), delay=int(delay), first_symbol=int(first_symbol),
        last_symbol=int(last_symbol), ser=float(sym_errors.sum() / a.size),
        symbol_errors=int(sym_errors.sum()), n_symbols=int(a.size), snr_db_pol=snr_pol,
        ber_pol=tuple(float(e / b.shape[1]) for e in bit_errors),
        ser_pol=tuple(float(e / a.shape[1]) for e in sym_errors),
        mi=float(np.mean([r[0] for r in rates])), gmi=float(np.mean([r[1] for r in rates])),
        mi_pol=tuple(float(r[0]) for r in rates), gmi_pol=tuple(float(r[1]) for r in rates),
        phase_slips=int(phase_slips))

def measure_data_aided(signal, bits, modulation: str, samples_per_symbol: int, guard: int, *,
                       taps: int, phase_window: int, block_length: int):
    """Medida com o DSP de simbolos conhecidos (kernels/dsp.py). `signal` (n_pol, N) ja passou
    pela compensacao de dispersao e, se houver, pelo filtro casado. `bits` (n_pol,
    n_sym*bits_por_simbolo). Devolve (CoherentMeasurement, simbolos equalizados)."""
    b = np.atleast_2d(np.asarray(bits)).astype(np.uint8)
    k_bits = BITS_PER_SYMBOL[modulation]
    a = np.array([map_bits(row, modulation) for row in b])
    y, delay = data_aided_receiver(signal, a, samples_per_symbol, guard, taps=taps,
                                   phase_window=phase_window, block_length=block_length)
    first, last = guard, a.shape[1] - guard
    m = symbol_metrics(y, b[:, first * k_bits:last * k_bits], modulation, gain=1.0, phase=0.0,
                       delay=delay, first_symbol=first, last_symbol=last)
    return m, y

# ------------------------------------------------------------------- receptor cego -------------
SYNC_BLOCK = 128          # simbolos por bloco da correlacao que casa a saida cega com os dados
SYNC_MAX_LAG = 16         # atraso maximo procurado [simbolos]: o equalizador desloca poucos simbolos

def align_to_transmitted(symbols, sent, modulation: str):
    """Casa as saidas do receptor cego (n_pol, n) com os simbolos transmitidos (n_pol, n).

    Para cada par (saida p, polarizacao transmitida q) e cada atraso de -SYNC_MAX_LAG a
    +SYNC_MAX_LAG simbolos, soma |correlacao|^2 em blocos de SYNC_BLOCK simbolos. A soma em blocos
    nao depende do giro de cada bloco, entao funciona mesmo com saltos de fase. Fica o par e o
    atraso de maior soma. Em cada bloco, o giro (multiplo de 2*pi/S) que alinha a saida aos dados:
    o giro mais frequente e tirado do bloco inteiro, e cada mudanca de giro entre blocos vizinhos
    conta como um salto de fase.
    Devolve (simbolos alinhados (n_pol, n), saltos de fase, atrasos por polarizacao).
    """
    y = np.atleast_2d(np.asarray(symbols, dtype=complex))
    a = np.atleast_2d(np.asarray(sent, dtype=complex))
    n_pol, n = a.shape
    block = min(SYNC_BLOCK, n)
    usable = (n // block) * block
    sector = 2.0 * np.pi / ROTATIONAL_SYMMETRY[modulation]
    lags = np.arange(-SYNC_MAX_LAG, SYNC_MAX_LAG + 1)
    score = np.zeros((n_pol, n_pol, lags.size))
    for p in range(n_pol):
        for j, lag in enumerate(lags):
            shifted = np.roll(y[p], -lag)[:usable]
            for q in range(n_pol):
                c = np.sum((shifted * np.conj(a[q, :usable])).reshape(-1, block), axis=1)
                score[p, q, j] = np.sum(np.abs(c) ** 2)
    best = score.max(axis=2)
    if n_pol == 2 and best[0, 1] + best[1, 0] > best[0, 0] + best[1, 1]:
        order = (1, 0)                    # a saida 1 e a polarizacao X transmitida
    else:
        order = tuple(range(n_pol))
    out = np.empty_like(a)
    slips, chosen = 0, []
    for q in range(n_pol):
        p = order[q]
        lag = int(lags[np.argmax(score[p, q])])
        row = np.roll(y[p], -lag)
        c = np.sum((row[:usable] * np.conj(a[q, :usable])).reshape(-1, block), axis=1)
        turns = np.round(np.angle(c) / sector).astype(int) % ROTATIONAL_SYMMETRY[modulation]
        slips += int(np.count_nonzero(np.diff(turns)))
        common = np.bincount(turns).argmax()
        out[q] = row * np.exp(-1j * common * sector)
        chosen.append(lag)
    return out, slips, tuple(chosen)

def measure_blind(symbols, bits, modulation: str, guard: int):
    """Medida do receptor cego. `symbols` (n_pol, n) sao a saida do receptor cego
    (kernels/dsp.py, `blind_receiver`), ja na escala da constelacao. Nao ha ganho nem fase
    ajustados com os dados: so a polarizacao, o atraso e um giro constante por polarizacao. Um
    salto de fase no meio do bloco aparece como erros e entra em `phase_slips`.
    Devolve (CoherentMeasurement, simbolos alinhados sem as guardas)."""
    b = np.atleast_2d(np.asarray(bits)).astype(np.uint8)
    k_bits = BITS_PER_SYMBOL[modulation]
    a = np.array([map_bits(row, modulation) for row in b])
    y, slips, lags = align_to_transmitted(symbols, a, modulation)
    first, last = guard, a.shape[1] - guard
    if last - first < 2:
        raise ValueError(f"guarda de {guard} simbolos por borda nao deixa simbolos para medir "
                         f"({a.shape[1]} simbolos na janela)")
    y = y[:, first:last]
    m = symbol_metrics(y, b[:, first * k_bits:last * k_bits], modulation, gain=1.0, phase=0.0,
                       delay=lags[0], first_symbol=first, last_symbol=last, phase_slips=slips)
    return m, y

# ------------------------------------------------------------------- sem DSP -------------------
def symbol_timing(signal, symbols, samples_per_symbol: int) -> int:
    """Atraso d [amostras] que maximiza |sum_k signal[k*sps + d]*conj(symbols[k])|, procurado em
    toda a janela (correlacao circular pela FFT). Sai entre -N/2 e N/2."""
    x = np.asarray(signal, dtype=complex)
    a = np.asarray(symbols, dtype=complex)
    sps = int(samples_per_symbol)
    if a.size * sps > x.size:
        raise ValueError(f"{a.size} simbolos x {sps} amostras nao cabem em {x.size} amostras")
    impulses = np.zeros(x.size, dtype=complex)
    impulses[:a.size * sps:sps] = a
    corr = np.fft.ifft(np.fft.fft(x) * np.conj(np.fft.fft(impulses)))
    d = int(np.argmax(np.abs(corr)))
    return d - x.size if d > x.size // 2 else d

def sample_symbols(signal, samples_per_symbol: int, delay: int, first_symbol: int,
                   last_symbol: int) -> np.ndarray:
    """Uma amostra por simbolo, no instante de decisao, dos simbolos first_symbol a last_symbol - 1."""
    x = np.asarray(signal)
    k = np.arange(first_symbol, last_symbol)
    return x[(k * int(samples_per_symbol) + delay) % x.size]

def measure_coherent(signal, bits, modulation: str, samples_per_symbol: int,
                     guard: int) -> CoherentMeasurement:
    """Mede um sinal I + i*Q (N,) de uma polarizacao contra os bits transmitidos, sem DSP: so o
    instante dos simbolos e um ganho complexo.

    Os `guard` simbolos de cada borda da janela ficam fora da medida, como no IMDD
    (config/resolve.py, `edge_guard_symbols`).
    """
    x = np.asarray(signal, dtype=complex)
    b = np.asarray(bits).astype(np.uint8)
    a = map_bits(b, modulation)
    k_bits = b.size // a.size
    first, last = guard, a.size - guard
    if last - first < 2:
        raise ValueError(f"guarda de {guard} simbolos por borda nao deixa simbolos para medir "
                         f"({a.size} simbolos na janela)")
    delay = symbol_timing(x, a, samples_per_symbol)
    r = sample_symbols(x, samples_per_symbol, delay, first, last)
    a_m = a[first:last]
    b_m = b[first * k_bits:last * k_bits]
    g = np.sum(r * np.conj(a_m)) / np.sum(np.abs(a_m) ** 2)
    if g == 0:
        raise ValueError("sinal recebido sem correlacao com os simbolos transmitidos (ganho nulo)")
    return symbol_metrics(r / g, b_m, modulation, gain=abs(g), phase=np.angle(g), delay=delay,
                          first_symbol=first, last_symbol=last)
