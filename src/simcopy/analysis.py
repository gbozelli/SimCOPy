"""Grandezas de interesse de um sinal, para graficos e telas (docs/mapa_arquitetura.md, 5.12).

So arrays, nunca figuras: o mesmo calculo serve ao simcopy.plotting, aos notebooks e a uma
interface. A convencao de frequencia (Agrawal, portadora exp(-i*w0*t); decisao de 2026-10-22) fica
escrita aqui uma vez: o chirp e -(1/2pi)*dphi/dt e o eixo optico e TimeGrid.optical_frequency.
Depende so de core.
"""
from __future__ import annotations
import numpy as np
from .core.signals import Waveform

def optical_power(w: Waveform) -> np.ndarray:
    """Potencia optica instantanea [W], somando as polarizacoes."""
    return np.sum(np.abs(w.samples) ** 2, axis=0)

def chirp(w: Waveform, floor: float = 0.02) -> np.ndarray:
    """Desvio de frequencia optica [Hz] da polarizacao com mais potencia, -(1/2pi)*dphi/dt
    (Agrawal, Eq. 3.4.1 na convencao do projeto). Onde P < floor*P_max a fase nao e definida e o
    valor e nan."""
    k = int(np.argmax(np.mean(np.abs(w.samples) ** 2, axis=1)))
    field = w.samples[k]
    p = np.abs(field) ** 2
    dv = -np.gradient(np.unwrap(np.angle(field)), w.grid.dt) / (2 * np.pi)
    dv[p < floor * p.max()] = np.nan
    return dv

def _smooth(s: np.ndarray, df: float, resolution: float | None) -> np.ndarray:
    """Media movel de largura `resolution` [Hz] (como a banda de resolucao de um analisador)."""
    if not resolution or resolution <= df:
        return s
    k = int(round(resolution / df))
    return np.convolve(s, np.ones(k) / k, mode="same")

def optical_spectrum(w: Waveform, resolution: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(frequencia optica relativa [Hz], densidade [W/Hz]) em ordem crescente, somando as
    polarizacoes, com janela de Hann (a fase do laser nao fecha na janela). Positivo = azul.
    `resolution` [Hz]: media movel do espectro (None = resolucao da grade, fs/N)."""
    n = w.grid.n_samples
    win = np.hanning(n)
    s = np.sum(np.abs(np.fft.fft(w.samples * win, axis=-1)) ** 2, axis=0)
    s *= w.grid.dt / np.sum(win ** 2)
    f = w.grid.optical_frequency
    order = np.argsort(f)
    return f[order], _smooth(s[order], w.grid.sampling_rate / n, resolution)

def electrical_spectrum(w: Waveform, resolution: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(frequencia >= 0 [Hz], densidade unilateral [unidade^2/Hz]) da parte real, sem o nivel
    medio, com janela de Hann. `resolution` como em optical_spectrum."""
    x = w.samples[0].real
    win = np.hanning(x.size)
    s = 2 * np.abs(np.fft.rfft((x - x.mean()) * win)) ** 2 * w.grid.dt / np.sum(win ** 2)
    return np.fft.rfftfreq(x.size, w.grid.dt), _smooth(s, w.grid.sampling_rate / x.size, resolution)

def eye(x, samples_per_symbol: int, delay: int, first_symbol: int, last_symbol: int,
        span: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Diagrama de olho: (tempo [simbolos], tracos), um traco de `span` simbolos por simbolo
    medido, centrado no instante de decisao (que fica em t = 0)."""
    x = np.asarray(x)
    sps = int(samples_per_symbol)
    start = np.arange(first_symbol, last_symbol) * sps + delay - (span * sps) // 2
    idx = (start[:, None] + np.arange(span * sps)[None, :]) % x.size
    t = (np.arange(span * sps) - (span * sps) // 2) / sps
    return t, x[idx]
