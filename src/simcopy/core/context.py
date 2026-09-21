"""Contexto: grade + aleatoriedade injetada. Substitui o np.random.seed global
que acoplava payload, sincronismo e ruido no mesmo stream (F-4)."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import numpy as np
from .grid import TimeGrid

def _stable_key(k):
    return k if isinstance(k, int) else int.from_bytes(
        hashlib.blake2b(k.encode(), digest_size=4).digest(), "big")

@dataclass(frozen=True)
class SimContext:
    grid: TimeGrid
    root_seed: int

    def rng_for(self, *keys):
        """Stream determinístico e independente por (canal, bloco, papel)."""
        return np.random.default_rng(np.random.SeedSequence(
            self.root_seed, spawn_key=tuple(_stable_key(k) for k in keys)))
