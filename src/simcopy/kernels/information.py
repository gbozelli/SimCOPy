"""Informacao mutua (MI) e informacao mutua generalizada (GMI) a partir dos simbolos recebidos
(docs/sistema_coerente.md, secao 12).

Nao estao no Agrawal. Referencia usual: A. Alvarado, E. Agrell, D. Lavery, R. Maher e P. Bayvel,
"Replacing the soft-decision FEC limit paradigm in the design of optical communication systems",
J. Lightw. Technol. 33, 4338 (2015).

Com simbolos equiprovaveis de uma constelacao de M = 2^k pontos, e um decodificador que supoe
ruido gaussiano circular de variancia s2 (a medida nos proprios simbolos), para N simbolos
recebidos y_n com x_n transmitido:

    MI  = k - (1/N) * sum_n log2 { sum_j exp[-(|y_n - c_j|^2 - |y_n - x_n|^2)/s2] }
    GMI = k - (1/N) * sum_n sum_i log2 { 1 + exp[-(1 - 2*b_ni)*L_ni] }
    L_ni = ln sum_{j: bit i de c_j = 0} exp(-|y_n - c_j|^2/s2)
         - ln sum_{j: bit i de c_j = 1} exp(-|y_n - c_j|^2/s2)

em bits por simbolo. A MI e a taxa alcancavel com decodificacao por simbolo, a GMI com
decodificacao bit a bit (codigo binario com entrelacamento), que e o caso dos codigos usados em
sistemas opticos. GMI <= MI <= k. Quando cada bit e uma decisao independente (BPSK, QPSK), as duas
coincidem. A GMI depende do mapeamento de bits: aqui e o de Gray de kernels/mapping.py.

O ruido real nao e gaussiano (imagem de quadratura, vazamento de portadora, nao linearidade): os
valores sao os que um decodificador gaussiano alcancaria com esse sinal, nao a capacidade do canal.
"""
from __future__ import annotations
import numpy as np
from scipy.special import logsumexp

def information_rates(received, sent, points, noise_variance: float | None = None,
                      chunk: int = 16384) -> tuple[float, float]:
    """(MI, GMI) em bits por simbolo, de simbolos recebidos `received` (n,) e transmitidos `sent`
    (n,), para a constelacao `points` na ordem do valor inteiro dos bits (kernels/mapping.py,
    `constellation`). `noise_variance`: variancia do ruido complexo; None = medida, mean|y - x|^2.
    """
    y = np.asarray(received, dtype=complex)
    x = np.asarray(sent, dtype=complex)
    c = np.asarray(points, dtype=complex)
    m = c.size
    k = int(round(np.log2(m)))
    if 2 ** k != m:
        raise ValueError(f"constelacao de {m} pontos: precisa ser potencia de 2")
    s2 = float(np.mean(np.abs(y - x) ** 2)) if noise_variance is None else float(noise_variance)
    if s2 <= 0.0:
        return float(k), float(k)
    labels = ((np.arange(m)[:, None] >> np.arange(k - 1, -1, -1)) & 1).astype(bool)     # (M, k)
    index = np.argmin(np.abs(x[:, None] - c[None, :]), axis=1) if x.size <= chunk else None
    mi_loss, gmi_loss = 0.0, 0.0
    for start in range(0, y.size, chunk):
        yb, xb = y[start:start + chunk], x[start:start + chunk]
        metric = -np.abs(yb[:, None] - c[None, :]) ** 2 / s2                              # (n, M)
        own = -np.abs(yb - xb) ** 2 / s2
        mi_loss += float(np.sum(logsumexp(metric, axis=1) - own))
        sent_index = index[start:start + chunk] if index is not None else \
            np.argmin(np.abs(xb[:, None] - c[None, :]), axis=1)
        bits = labels[sent_index]                                                          # (n, k)
        for i in range(k):
            zero = logsumexp(metric[:, ~labels[:, i]], axis=1)
            one = logsumexp(metric[:, labels[:, i]], axis=1)
            llr = np.where(bits[:, i], one - zero, zero - one)       # a favor do bit transmitido
            gmi_loss += float(np.sum(np.logaddexp(0.0, -llr)))
    scale = 1.0 / (y.size * np.log(2.0))
    return k - mi_loss * scale, k - gmi_loss * scale
