"""Multiplexador e demultiplexador WDM (docs/mux_demux.md).

Codigo do professor (draft.py): o mux e a soma dos campos dos canais, cada um ja deslocado em
frequencia (`Freq_offset` do laser), e o demux e um filtro optico passa-banda (`optical_filter`,
mascara retangular na FFT). Agrawal, secao 6.3.1: com o filtro ajustado no canal m, chega ao
detector P = P_m + soma_n T_mn*P_n (interferencia linear fora da banda, Eq. 6.3.1).

Frequencias no eixo da FFT do NumPy (o sentido fisico depende da duvida M1 de
docs/moduladores.md). Os filtros sao aplicados no dominio da frequencia, como no professor.
"""
from __future__ import annotations
import numpy as np

def frequency_shift(field, shift, time):
    """Desloca o campo (n_pol, N) de `shift` [Hz] no eixo da FFT: campo*exp(+i*2*pi*shift*t)."""
    return np.asarray(field) * np.exp(2j * np.pi * shift * np.asarray(time))

def multiplex(fields, shifts, time):
    """Campo agregado: soma dos campos dos canais, cada um deslocado de shifts[k] [Hz]."""
    return sum(frequency_shift(c, s, time) for c, s in zip(fields, shifts))

def bandpass_response(freq, center, bandwidth, kind, order=1):
    """Transmissao de campo H(f) de um filtro optico centrado em `center` [Hz].

    kind = "rectangular": 1 dentro de [center - B/2, center + B/2), 0 fora (o `optical_filter`
    do professor); kind = "supergaussian": exp(-(ln 2/2)*(2*(f - center)/B)^(2*order)), com
    |H|^2 = 1/2 nas bordas (banda de 3 dB = B); order = 1 e a gaussiana.
    """
    x = 2.0 * (np.asarray(freq) - center) / bandwidth
    if kind == "rectangular":
        return ((x >= -1.0) & (x < 1.0)).astype(float)
    if kind == "supergaussian":
        return np.exp(-0.5 * np.log(2.0) * x ** (2 * order))
    raise ValueError(f"filtro desconhecido: {kind!r}")

def apply_optical_filter(field, response):
    """Filtra cada polarizacao do campo (n_pol, N) pela transmissao H(f) no eixo da FFT."""
    return np.fft.ifft(np.fft.fft(np.asarray(field), axis=-1) * response, axis=-1)
