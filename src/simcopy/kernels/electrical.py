"""Sinal eletrico do transmissor IMDD (docs/sistema_imdd.md).

Aula (docs/04-class/context.md): "criar sequencia pseudoaleatoria, converter para sinal eletrico
(nao se usa Nyquist, e sim um filtro de Bessel)". O filtro e o de kernels/filters.py.
"""
from __future__ import annotations
import numpy as np

def pam_drive_levels(bits, bits_per_symbol: int) -> np.ndarray:
    """Niveis eletricos em [0, 1] de um sinal PAM de 2^bits_per_symbol niveis igualmente
    espacados: 0 e 1 no OOK; 0, 1/3, 2/3 e 1 no PAM4. Os bits de cada simbolo escolhem o nivel
    pelo codigo de Gray refletido de kernels/mapping.py (PAM4: 00, 01, 11, 10, do mais baixo ao
    mais alto), para que um erro para o nivel vizinho troque um bit so."""
    from .mapping import gray_to_index
    k = int(bits_per_symbol)
    b = np.asarray(bits).astype(int)
    if k == 1:
        return b.astype(float)
    if b.size % k:
        raise ValueError(f"{b.size} bits nao formam simbolos de {k} bits")
    return gray_to_index(b.reshape(-1, k)) / (2.0 ** k - 1.0)

def nrz(levels, samples_per_symbol: int, n_samples: int | None = None) -> np.ndarray:
    """Forma de onda NRZ: cada nivel repetido `samples_per_symbol` vezes (notebook 07).

    Com `n_samples` maior que len(levels)*samples_per_symbol (janela arredondada para potencia de
    2 no resolve), o ultimo nivel e mantido ate o fim; essas amostras ficam fora da medida.
    """
    x = np.repeat(np.asarray(levels, dtype=float), samples_per_symbol)
    if n_samples is None or n_samples == x.size:
        return x
    if n_samples < x.size:
        raise ValueError(f"{x.size} amostras de simbolos nao cabem em {n_samples}")
    return np.concatenate([x, np.full(n_samples - x.size, x[-1])])
