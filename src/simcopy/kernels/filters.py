"""Filtros eletricos passa-baixas, aplicados na FFT (docs/fotodetectores.md, docs/sistema_imdd.md).

Estavam em kernels/photodetector.py; vieram para ca para servir tambem ao acionamento do
transmissor e ao S21 do modulador com a mesma definicao de banda: -3 dB na banda pedida, resposta
analogica (Bessel normalizado em modulo). Aplicar na FFT trata a janela como periodica, como a
fibra, e nao tem o transiente de partida do lfilter.
"""
from __future__ import annotations
import numpy as np
from scipy import signal

def lowpass_response(freq, bandwidth, kind, order):
    """Resposta em frequencia H(f) de um passa-baixas analogico com -3 dB em `bandwidth` [Hz].

    kind = "bessel" (fase quase linear) ou "butterworth"; `freq` no eixo da FFT [Hz].
    """
    wn = 2 * np.pi * bandwidth
    if kind == "bessel":
        b, a = signal.bessel(order, wn, analog=True, norm="mag")
    elif kind == "butterworth":
        b, a = signal.butter(order, wn, analog=True)
    else:
        raise ValueError(f"filtro desconhecido: {kind!r}")
    return signal.freqs(b, a, 2 * np.pi * np.asarray(freq))[1]

def noise_bandwidth(response, freq):
    """Banda equivalente de ruido delta_f = integral de |H|^2 nas frequencias positivas (Eq. 4.4.4),
    normalizada por |H(0)|^2, a partir da resposta amostrada na grade da FFT."""
    df = abs(freq[1] - freq[0])
    h2 = np.abs(response) ** 2
    return 0.5 * np.sum(h2) * df / h2[0]

def apply_response(x, response):
    """Filtra um sinal real pela resposta H(f) dada no eixo da FFT."""
    return np.fft.ifft(np.fft.fft(x) * response).real
