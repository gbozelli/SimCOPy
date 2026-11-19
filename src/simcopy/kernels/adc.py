"""Conversor analogico-digital: amostragem e quantizacao (docs/sistema_coerente.md, secao 12).

VPI (documento de parametros, secao 7, `Rx_DualPol`): `SamplingRate = 2*SymbolRate`,
`TransferFunctionADC = Bipolar`, `Resolution = 12 bit`, `Range = 1.0`. O documento nao diz em que
escala a corrente chega ao conversor, entao o fundo de escala aqui e dado em multiplos do valor
rms do sinal (controle de ganho ideal): duvida 65.

Quantizador uniforme bipolar de b bits, com fundo de escala +-F: 2^b niveis de largura
D = 2F/2^b, no centro de cada intervalo, e ceifamento em +-F. Para um sinal de valor rms s que
nao ceifa, o ruido de quantizacao tem variancia D^2/12 e
    SQNR = 6,02*b + 4,77 - 20*log10(F/s)  [dB].
"""
from __future__ import annotations
import numpy as np

def quantize(x, bits: int, full_scale: float) -> np.ndarray:
    """Quantiza um sinal real em 2^bits niveis uniformes entre -full_scale e +full_scale."""
    if bits < 1:
        raise ValueError(f"bits = {bits}: o quantizador pede pelo menos 1 bit")
    levels = 2 ** int(bits)
    step = 2.0 * full_scale / levels
    index = np.clip(np.floor(np.asarray(x, dtype=float) / step), -levels // 2, levels // 2 - 1)
    return (index + 0.5) * step

def sample_and_quantize(signal, decimation: int, bits: int | None, full_scale_rms: float,
                        offset: int = 0) -> np.ndarray:
    """Amostra um sinal complexo (n_pol, N) a cada `decimation` amostras, a partir de `offset`, e
    quantiza a parte real e a imaginaria de cada polarizacao em separado, com fundo de escala
    `full_scale_rms` vezes o valor rms da componente. `bits` None ou 0: sem quantizacao."""
    x = np.atleast_2d(np.asarray(signal, dtype=complex))[:, offset::int(decimation)]
    if not bits:
        return x.copy()
    out = np.empty_like(x)
    for p, row in enumerate(x):
        parts = []
        for part in (row.real, row.imag):
            rms = np.sqrt(np.mean(part ** 2))
            parts.append(quantize(part, bits, full_scale_rms * rms) if rms > 0 else part)
        out[p] = parts[0] + 1j * parts[1]
    return out
