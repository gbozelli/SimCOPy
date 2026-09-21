"""Tipos de sinal (secao 4.2.1 do documento).

Tres representacoes com dtypes e comprimentos diferentes, ligadas por uma
SignalGeometry unica. O campo `domain` existe como metadado para verificacao de
pre-condicao e escolha de visualizacao; nunca para escolher fisica dentro de uma
mesma funcao.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from enum import Enum
import numpy as np
from .grid import TimeGrid

class Domain(Enum):
    ELECTRICAL = "electrical"
    OPTICAL = "optical"

class SignalError(ValueError): pass

@dataclass(frozen=True)
class SignalGeometry:
    """O sizeof() da proposta original, promovido a objeto de primeira classe."""
    n_symbols: int
    bits_per_symbol: int
    samples_per_symbol: int
    n_pol: int

    def __post_init__(self):
        if self.n_pol not in (1, 2):
            raise SignalError(f"n_pol deve ser 1 ou 2, recebido {self.n_pol}")
        if self.samples_per_symbol < 1:
            raise SignalError("samples_per_symbol deve ser inteiro >= 1 ([D-4])")

    @property
    def n_bits(self):    return self.n_symbols * self.bits_per_symbol * self.n_pol
    @property
    def n_samples(self): return self.n_symbols * self.samples_per_symbol

@dataclass(frozen=True)
class Bits:
    data: np.ndarray
    geometry: SignalGeometry
    def __post_init__(self):
        exp = (self.geometry.n_pol, self.geometry.n_bits // self.geometry.n_pol)
        if self.data.shape != exp:
            raise SignalError(f"Bits: esperado {exp}, recebido {self.data.shape}")

@dataclass(frozen=True)
class Symbols:
    data: np.ndarray
    geometry: SignalGeometry
    def __post_init__(self):
        exp = (self.geometry.n_pol, self.geometry.n_symbols)
        if self.data.shape != exp:
            raise SignalError(f"Symbols: esperado {exp}, recebido {self.data.shape}")

@dataclass(frozen=True)
class Waveform:
    """Campo amostrado, sempre (n_pol, n_samples) complexo. Sinal eletrico real e
    guardado como complexo de parte imaginaria nula ([D-5])."""
    samples: np.ndarray
    grid: TimeGrid
    domain: Domain
    geometry: SignalGeometry
    f0: float = 0.0
    channel_id: int | None = None

    def __post_init__(self):
        if self.samples.ndim != 2 or self.samples.shape[0] not in (1, 2):
            raise SignalError(f"esperado (1|2, N), recebido {self.samples.shape}")
        if self.samples.shape[1] != self.grid.n_samples:
            raise SignalError(f"{self.samples.shape[1]} amostras nao batem com a "
                              f"grade ({self.grid.n_samples})")
        if not np.iscomplexobj(self.samples):
            raise SignalError("samples deve ser complexo ([D-5])")

    @property
    def n_pol(self): return self.samples.shape[0]

    @property
    def power(self):
        """Potencia media somando TODAS as polarizacoes (corrige C-5, que no
        codigo antigo somava a polarizacao X duas vezes)."""
        return float(np.sum(np.mean(np.abs(self.samples) ** 2, axis=1)))

    def with_samples(self, samples): return replace(self, samples=samples)

    def require(self, domain, who):
        """Pre-condicao. Unico uso legitimo de `domain`."""
        if self.domain is not domain:
            raise SignalError(f"{who} espera entrada {domain.value}, "
                              f"recebeu {self.domain.value}")
