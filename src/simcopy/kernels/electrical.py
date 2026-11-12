"""Sinal eletrico do transmissor IMDD (docs/sistema_imdd.md).

Aula (docs/04-class/context.md): "criar sequencia pseudoaleatoria, converter para sinal eletrico
(nao se usa Nyquist, e sim um filtro de Bessel)". O filtro e o de kernels/filters.py.
"""
from __future__ import annotations
import numpy as np

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
