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


@dataclass(frozen=True)
class DigitalSignal:
    """Quarta representacao: o sinal na taxa do DSP.

    E a entrada do DAC no transmissor e a saida do ADC no receptor (no arquivo
    VPI, SamplingRate = 2*SymbolRate e 12 bits). Fica FORA da grade global de
    simulacao: tem taxa propria, samples_per_symbol * symbol_rate, tipicamente
    2 amostras por simbolo contra as 16 da grade.

    Nao tem f0. Depois da deteccao coerente o sinal esta em banda base em relacao
    ao LO; o offset residual (LO_Offset, ruido de fase) e conteudo do sinal, que o
    DSP estima e remove, e nao um metadado de posicao na grade WDM.

    A quantizacao e um PROCESSO (feito pelo componente ADC), nao um metodo deste
    tipo. `resolution_bits` apenas registra se o sinal ja foi quantizado e com
    quantos bits; None significa amostras ideais, sem quantizacao.
    """
    samples: np.ndarray                   # (n_pol, n_symbols * samples_per_symbol)
    symbol_rate: float
    geometry: SignalGeometry              # geometry.samples_per_symbol = taxa do DSP
    resolution_bits: int | None = None
    channel_id: int | None = None

    def __post_init__(self):
        exp = (self.geometry.n_pol, self.geometry.n_samples)
        if self.samples.shape != exp:
            raise SignalError(f"DigitalSignal: esperado {exp}, recebido "
                              f"{self.samples.shape}")
        if not np.iscomplexobj(self.samples):
            raise SignalError("samples deve ser complexo (I + jQ)")
        if self.symbol_rate <= 0:
            raise SignalError("symbol_rate deve ser positiva")
        if self.resolution_bits is not None and self.resolution_bits < 1:
            raise SignalError("resolution_bits deve ser >= 1, ou None (sem quantizacao)")

    @property
    def samples_per_symbol(self) -> int:
        return self.geometry.samples_per_symbol

    @property
    def sampling_rate(self) -> float:
        """Taxa do DSP, derivada. Nao e a sampling_rate global da simulacao."""
        return self.samples_per_symbol * self.symbol_rate

    @property
    def n_pol(self) -> int:
        return self.samples.shape[0]

    def symbol_samples(self, offset: int = 0) -> np.ndarray:
        """Uma amostra por simbolo, a partir de `offset` (instante de decisao)."""
        return self.samples[:, offset::self.samples_per_symbol]
