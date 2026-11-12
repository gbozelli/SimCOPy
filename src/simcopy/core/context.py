"""Contexto: grade + aleatoriedade injetada. Substitui o np.random.seed global
que acoplava payload, sincronismo e ruido no mesmo stream (F-4)."""

from __future__ import annotations
from dataclasses import dataclass
import hashlib
import numpy as np
from .grid import TimeGrid

def _stable_key(k: int | str) -> int:
    return k if isinstance(k, int) else int.from_bytes(
        hashlib.blake2b(k.encode(), digest_size=4).digest(), "big")

@dataclass(frozen=True)
class SimContext:
    grid: TimeGrid
    root_seed: int

    def rng_for(self, *keys: int | str) -> np.random.Generator:
        """Stream determinístico e independente por (canal, bloco, papel)."""
        return np.random.default_rng(np.random.SeedSequence(
            self.root_seed, spawn_key=tuple(_stable_key(k) for k in keys)))

def rng_from_seed(seed: int, *keys: int | str) -> np.random.Generator:
    """Stream que depende só de `seed` e das chaves, não da semente raiz.

    É o dos bits de cada canal: a semente do canal define os bits (como o RandomNumberSeed do
    VPI) e a semente raiz define os ruídos (decisão de 2026-11-12). Trocar root_seed muda os
    ruídos e mantém os bits.
    """
    return np.random.default_rng(np.random.SeedSequence(
        seed, spawn_key=tuple(_stable_key(k) for k in keys)))
