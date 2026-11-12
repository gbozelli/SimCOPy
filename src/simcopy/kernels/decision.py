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
