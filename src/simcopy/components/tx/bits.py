"""Fonte de bits do transmissor (docs/sistema_imdd.md, secao 2)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...core.context import SimContext, rng_from_seed
from ...core.signals import Bits, SignalGeometry

@dataclass(frozen=True)
class BitSource:
    """Bits equiprovaveis, um stream por canal, definido so pela semente do canal: os mesmos bits
    qualquer que seja a semente raiz, que fica para os ruidos (decisao de 2026-11-12; como o
    RandomNumberSeed do VPI). Fonte sem entrada: `x` e ignorado."""
    geometry: SignalGeometry
    seed: int
    id: str = "tx.bits"
    output_unit: ClassVar[str | None] = None

    def __call__(self, x, ctx: SimContext) -> Bits:
        g = self.geometry
        rng = rng_from_seed(self.seed, "payload")
        data = rng.integers(0, 2, (g.n_pol, g.n_bits // g.n_pol), dtype=np.uint8)
        return Bits(data=data, geometry=g)
