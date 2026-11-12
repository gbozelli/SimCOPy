"""Controle de potencia de lancamento (arquitetura, anexo A.4: PowerControl no fim do build_tx)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, Waveform

@dataclass(frozen=True)
class PowerControl:
    """Escala o campo para a potencia media `power` [W], somando as polarizacoes."""
    power: float
    id: str = "tx.power"
    output_unit: ClassVar[str] = "sqrt(W)"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, self.id)
        if x.power <= 0.0:
            raise SignalError(f"{self.id}: potencia de entrada nula; nao ha o que escalar.")
        return x.with_samples(x.samples * np.sqrt(self.power / x.power))
