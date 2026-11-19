"""Formatacao de pulso de Nyquist do transmissor coerente (docs/sistema_coerente.md, secao 2).

Codigo do professor (draft.py, `DAC_Nyquist`): os simbolos viram um trem de impulsos, um a cada
`samples_per_symbol` amostras, e o espectro passa por uma janela plana ate (1 - rolloff)*Rs/2 que
cai em cosseno ate (1 + rolloff)*Rs/2. E o espectro de cosseno levantado:
    H(f) = 1                                             |f| <= (1 - r)*Rs/2
    H(f) = (1/2)*{1 + cos[pi*(|f| - (1 - r)*Rs/2)/(r*Rs)]}   ate (1 + r)*Rs/2
    H(f) = 0                                             acima
Como a soma das replicas de H espacadas de Rs e constante (criterio de Nyquist), o sinal vale
exatamente o simbolo no centro de cada simbolo: nao ha interferencia entre simbolos na saida do
transmissor. Entre os simbolos o sinal passa dos niveis da constelacao.

O filtro e aplicado na FFT, que trata a janela como periodica, como no professor.
"""
from __future__ import annotations
import numpy as np

def raised_cosine_response(freq, symbol_rate: float, rolloff: float) -> np.ndarray:
    """H(f) de cosseno levantado no eixo da FFT `freq` [Hz]; H(0) = 1, 0 <= rolloff <= 1."""
    if not 0.0 <= rolloff <= 1.0:
        raise ValueError(f"rolloff = {rolloff} fora de [0, 1]")
    f = np.abs(np.asarray(freq, dtype=float))
    f1, f2 = 0.5 * symbol_rate * (1.0 - rolloff), 0.5 * symbol_rate * (1.0 + rolloff)
    h = np.zeros_like(f)
    h[f <= f1] = 1.0
    edge = (f > f1) & (f < f2)
    if rolloff > 0.0:
        h[edge] = 0.5 * (1.0 + np.cos(np.pi * (f[edge] - f1) / (rolloff * symbol_rate)))
    return h

def root_raised_cosine_response(freq, symbol_rate: float, rolloff: float) -> np.ndarray:
    """Raiz do cosseno levantado: o filtro do transmissor do VPI (`NyquistResponse =
    squareRootRaisedCosine`). O criterio de Nyquist so fecha depois do filtro casado do receptor,
    que tem a mesma resposta: o produto dos dois e o cosseno levantado."""
    return np.sqrt(raised_cosine_response(freq, symbol_rate, rolloff))

def shape_pulses(symbols, samples_per_symbol: int, response) -> np.ndarray:
    """Simbolos complexos -> sinal I(t) + i*Q(t) com `samples_per_symbol` amostras por simbolo.

    `response` e H(f) no eixo da FFT do sinal de saida (len(symbols)*samples_per_symbol pontos).
    O trem de impulsos leva o fator `samples_per_symbol`, como no professor, para que o sinal
    valha o simbolo no centro do simbolo (amostra k*samples_per_symbol).
    """
    a = np.asarray(symbols, dtype=complex)
    sps = int(samples_per_symbol)
    response = np.asarray(response)
    if response.size != a.size * sps:
        raise ValueError(f"H(f) tem {response.size} pontos; {a.size} simbolos x {sps} amostras "
                         f"dao {a.size * sps}")
    impulses = np.zeros(a.size * sps, dtype=complex)
    impulses[::sps] = sps * a
    return np.fft.ifft(np.fft.fft(impulses) * response)
