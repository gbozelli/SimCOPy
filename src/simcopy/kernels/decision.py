"""Decisao e fator Q de um sinal OOK (docs/sistema_imdd.md, secao 3).

Agrawal, Sistemas de Comunicacao por Fibra Optica, Eq. 4.6.6-4.6.11: com ruido gaussiano, niveis
I1, I0 e desvios s1, s0 nos bits 1 e 0,
    I_D = (s0*I1 + s1*I0)/(s0 + s1),   Q = (I1 - I0)/(s1 + s0),   BER = (1/2)*erfc(Q/sqrt(2)).
E a mesma conta do notebook 05 (secao 4), que a confere contra a BER contada.

Por que Q e nao contagem (anotacao do notebook 07): a BER de interesse em IMDD (1e-9 a 1e-12)
exigiria 1e10 bits ou mais; Q se estima com algumas centenas. Com interferencia entre simbolos
(dispersao, chirp), as distribuicoes deixam de ser gaussianas e Q mede a abertura do olho: a BER
pela Eq. 4.6.10 vira estimativa.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import erfc

@dataclass(frozen=True)
class QMeasurement:
    q: float
    ber_q: float                  # (1/2)*erfc(Q/sqrt(2)), Eq. 4.6.10
    ber_counted: float            # erros / bits, com o limiar otimo
    errors: int
    n_bits: int                   # bits usados na medida (sem as guardas)
    threshold: float              # I_D, Eq. 4.6.6
    level_one: float              # I1
    level_zero: float             # I0
    sigma_one: float              # s1
    sigma_zero: float             # s0
    delay: int                    # amostras do inicio do simbolo transmitido ao instante de decisao
    first_symbol: int             # primeiro simbolo medido
    last_symbol: int              # ultimo simbolo medido + 1

def ber_from_q(q):
    """BER = (1/2)*erfc(Q/sqrt(2)) (Agrawal, Eq. 4.6.10)."""
    return 0.5 * erfc(np.asarray(q, dtype=float) / np.sqrt(2.0))

def optimal_threshold(level_one, level_zero, sigma_one, sigma_zero):
    """Limiar que minimiza a BER com ruido gaussiano (Agrawal, Eq. 4.6.6)."""
    total = sigma_one + sigma_zero
    if total == 0.0:
        return 0.5 * (level_one + level_zero)
    return (sigma_zero * level_one + sigma_one * level_zero) / total

def level_statistics(samples, bits):
    """(I1, I0, s1, s0): media e desvio das amostras dos bits 1 e dos bits 0."""
    s = np.asarray(samples, dtype=float)
    b = np.asarray(bits).astype(bool)
    return s[b].mean(), s[~b].mean(), s[b].std(), s[~b].std()

def q_factor(level_one, level_zero, sigma_one, sigma_zero):
    """Q = (I1 - I0)/(s1 + s0) (Agrawal, Eq. 4.6.11); infinito sem ruido nem interferencia."""
    total = sigma_one + sigma_zero
    return np.inf if total == 0.0 else (level_one - level_zero) / total

def measure_ook(signal, bits, samples_per_symbol: int, guard: int,
                max_delay: int | None = None) -> QMeasurement:
    """Fator Q de um sinal OOK amostrado uma vez por bit, no instante que maximiza Q.

    O instante e procurado entre 0 e `max_delay` amostras depois do inicio de cada simbolo
    (padrao: 2 simbolos), com os bits transmitidos conhecidos, como no notebook 05: cobre o atraso
    de grupo dos filtros e a resposta do laser. Os `guard` simbolos de cada borda da janela ficam
    fora da medida (config/resolve.py, `edge_guard_symbols`).
    """
    x = np.asarray(signal, dtype=float)
    b = np.asarray(bits).astype(bool)
    sps = int(samples_per_symbol)
    if max_delay is None:
        max_delay = 2 * sps
    first, last = guard, b.size - guard
    if last - first < 2:
        raise ValueError(f"guarda de {guard} simbolos por borda nao deixa bits para medir "
                         f"({b.size} simbolos na janela)")
    k = np.arange(first, last)
    bk = b[k]
    if bk.all() or not bk.any():
        raise ValueError("os bits medidos tem um nivel so; Q nao e definido")
    best = None
    for d in range(max_delay):
        s = x[(k * sps + d) % x.size]
        stats = level_statistics(s, bk)
        q = q_factor(*stats)
        if best is None or q > best[0]:
            best = (q, d, s, stats)
    q, d, s, (i1, i0, s1, s0) = best
    threshold = optimal_threshold(i1, i0, s1, s0)
    errors = int(np.count_nonzero((s > threshold) != bk))
    return QMeasurement(q=float(q), ber_q=float(ber_from_q(q)), ber_counted=errors / bk.size,
                        errors=errors, n_bits=int(bk.size), threshold=float(threshold),
                        level_one=float(i1), level_zero=float(i0), sigma_one=float(s1),
                        sigma_zero=float(s0), delay=int(d), first_symbol=int(first),
                        last_symbol=int(last))

# ------------------------------------------------------------------- PAM de varios niveis ------
@dataclass(frozen=True)
class PAMMeasurement:
    """Medida de um sinal PAM de M niveis (M = 4 no PAM4). `q` e o do pior olho, para servir aos
    mesmos resumos e graficos do OOK; `eye_q` tem os M - 1."""
    q: float                      # menor fator Q entre os olhos
    ber_q: float                  # BER pelos fatores Q de todos os olhos (ver `measure_pam`)
    ber_counted: float            # erros / bits, com os limiares otimos
    errors: int
    n_bits: int
    ser: float                    # erros de simbolo / simbolos
    symbol_errors: int
    levels: tuple                 # media de cada nivel, do mais baixo ao mais alto
    sigmas: tuple                 # desvio de cada nivel
    thresholds: tuple             # M - 1 limiares
    eye_q: tuple                  # M - 1 fatores Q, um por olho
    delay: int                    # amostras do inicio do simbolo transmitido ao instante de decisao
    first_symbol: int
    last_symbol: int

def _pam_statistics(samples, index, n_levels: int):
    """(niveis, desvios, limiares, Q de cada olho) das amostras `samples` com nivel transmitido
    `index`. Cada olho, entre dois niveis vizinhos, e tratado como o olho do OOK: limiar otimo da
    Eq. 4.6.6 e fator Q da Eq. 4.6.11 do Agrawal."""
    levels = np.array([samples[index == i].mean() for i in range(n_levels)])
    sigmas = np.array([samples[index == i].std() for i in range(n_levels)])
    thresholds = np.array([optimal_threshold(levels[i + 1], levels[i], sigmas[i + 1], sigmas[i])
                           for i in range(n_levels - 1)])
    eye_q = np.array([q_factor(levels[i + 1], levels[i], sigmas[i + 1], sigmas[i])
                      for i in range(n_levels - 1)])
    return levels, sigmas, thresholds, eye_q

def measure_pam(signal, bits, bits_per_symbol: int, samples_per_symbol: int,
                guard: int) -> PAMMeasurement:
    """Fatores Q, limiares e erros de um sinal PAM de 2^bits_per_symbol niveis (PAM4: 2 bits).

    Os bits de cada simbolo escolhem o nivel pelo codigo de Gray de kernels/mapping.py (PAM4:
    00, 01, 11, 10, do mais baixo ao mais alto), como no acionamento (kernels/electrical.py).

    Instante de decisao: o atraso de maior correlacao entre o sinal e a sequencia de niveis,
    procurado na janela inteira, e depois o instante, a ate um simbolo dele, que maximiza o Q do
    pior olho. A procura na janela inteira cobre o atraso de um canal WDM fora do centro da banda.

    BER pelos fatores Q, com ruido gaussiano em cada nivel e codigo de Gray (um erro para o nivel
    vizinho troca um bit): cada olho i contribui com (1/2)*erfc(Q_i/sqrt(2)) para os erros dos
    seus dois niveis, entao
        SER = (1/M)*sum_i erfc(Q_i/sqrt(2)),   BER = SER/log2(M).
    Com M = 2 volta a Eq. 4.6.10 do Agrawal. A generalizacao para M niveis nao esta no Agrawal.
    """
    from .mapping import gray_to_index, index_to_gray
    x = np.asarray(signal, dtype=float)
    k = int(bits_per_symbol)
    b = np.asarray(bits).astype(np.uint8).reshape(-1, k)
    n_levels = 2 ** k
    index = gray_to_index(b)
    sps = int(samples_per_symbol)
    first, last = guard, index.size - guard
    if last - first < 4 * n_levels:
        raise ValueError(f"guarda de {guard} simbolos por borda nao deixa simbolos para medir "
                         f"({index.size} simbolos na janela)")
    sym = np.arange(first, last)
    idx = index[sym]
    if np.unique(idx).size < n_levels:
        raise ValueError("os simbolos medidos nao tem todos os niveis; Q nao e definido")
    # atraso grosseiro: correlacao circular do sinal com o trem de niveis (sem a media)
    impulses = np.zeros(x.size)
    impulses[sym * sps] = idx - idx.mean()
    corr = np.fft.ifft(np.fft.fft(x - x.mean()) * np.conj(np.fft.fft(impulses))).real
    coarse = int(np.argmax(corr))
    best = None
    for d in range(coarse - sps, coarse + sps + 1):
        s = x[(sym * sps + d) % x.size]
        stats = _pam_statistics(s, idx, n_levels)
        worst = float(np.min(stats[3]))
        if best is None or worst > best[0]:
            best = (worst, d, s, stats)
    worst, d, s, (levels, sigmas, thresholds, eye_q) = best
    decided = np.searchsorted(np.sort(thresholds), s)
    sym_errors = int(np.count_nonzero(decided != idx))
    bit_errors = int(np.count_nonzero(index_to_gray(decided, k) != b[sym]))
    ser_q = float(np.sum(2.0 * ber_from_q(eye_q)) / n_levels)
    d = d - x.size if d > x.size // 2 else d
    return PAMMeasurement(
        q=worst, ber_q=ser_q / k, ber_counted=bit_errors / (idx.size * k), errors=bit_errors,
        n_bits=int(idx.size * k), ser=sym_errors / idx.size, symbol_errors=sym_errors,
        levels=tuple(float(v) for v in levels), sigmas=tuple(float(v) for v in sigmas),
        thresholds=tuple(float(v) for v in thresholds), eye_q=tuple(float(v) for v in eye_q),
        delay=int(d), first_symbol=int(first), last_symbol=int(last))
