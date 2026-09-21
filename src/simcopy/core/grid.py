"""Grade temporal compartilhada. No codigo antigo cada funcao reconstruia t e f
com convencoes ligeiramente diferentes (t[1]-t[0] aqui, t[2]-t[1] la)."""
from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
import numpy as np

@dataclass(frozen=True)
class TimeGrid:
    sampling_rate: float
    n_samples: int
    t0: float = 0.0

    def __post_init__(self):
        if self.sampling_rate <= 0: raise ValueError("sampling_rate deve ser positiva")
        if self.n_samples <= 0:     raise ValueError("n_samples deve ser positivo")

    @property
    def dt(self):       return 1.0 / self.sampling_rate
    @property
    def duration(self): return self.n_samples * self.dt

    @cached_property
    def time(self):  return self.t0 + np.arange(self.n_samples) * self.dt
    @cached_property
    def freq(self):  return np.fft.fftfreq(self.n_samples, self.dt)
    @cached_property
    def omega(self): return 2.0 * np.pi * self.freq
