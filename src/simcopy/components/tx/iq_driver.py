"""Acionamento eletrico do transmissor coerente: bits -> simbolos -> I(t) + i*Q(t)
(docs/sistema_coerente.md, secao 2)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...core.context import SimContext
from ...core.signals import Bits, Domain, SignalError, Waveform
from ...kernels.mapping import BITS_PER_SYMBOL, map_bits
from ...kernels.pulse_shaping import (raised_cosine_response, root_raised_cosine_response,
                                      shape_pulses)

PULSES = {"raised_cosine": raised_cosine_response,
          "root_raised_cosine": root_raised_cosine_response}

@dataclass(frozen=True)
class NyquistDriver:
    """Bits -> sinal eletrico complexo I(t) + i*Q(t), adimensional, na grade da simulacao, uma
    linha por polarizacao.

    Mapeia os bits (kernels/mapping.py: BPSK, QPSK ou 16QAM, niveis em [-1, 1]) e formata o pulso
    (kernels/pulse_shaping.py) com fator `rolloff`:
    - "raised_cosine": cosseno levantado inteiro, como o `DAC_Nyquist` do professor. No centro de
      cada simbolo (amostra k*samples_per_symbol) o sinal vale o simbolo;
    - "root_raised_cosine": a raiz, como o transmissor do VPI. O sinal so vale o simbolo depois do
      filtro casado do receptor.
    A janela da simulacao tem que conter um numero inteiro de simbolos: o filtro e periodico.
    """
    symbol_rate: float
    modulation: str                   # "BPSK", "QPSK" ou "16QAM"
    rolloff: float
    channel_id: int | None = None
    id: str = "tx.driver"
    pulse: str = "raised_cosine"
    output_unit: ClassVar[str] = "adimensional"

    def __call__(self, bits: Bits, ctx: SimContext) -> Waveform:
        g = bits.geometry
        if self.modulation not in BITS_PER_SYMBOL:
            raise SignalError(f"{self.id}: formato {self.modulation!r} sem mapeamento; "
                              f"disponiveis: {sorted(BITS_PER_SYMBOL)}.")
        if self.pulse not in PULSES:
            raise SignalError(f"{self.id}: pulso {self.pulse!r}; disponiveis: {sorted(PULSES)}.")
        if g.bits_per_symbol != BITS_PER_SYMBOL[self.modulation]:
            raise SignalError(f"{self.id}: {self.modulation} pede "
                              f"{BITS_PER_SYMBOL[self.modulation]} bit(s) por simbolo; a geometria "
                              f"tem {g.bits_per_symbol}.")
        grid = ctx.grid
        if g.n_samples != grid.n_samples:
            raise SignalError(f"{self.id}: a janela tem {grid.n_samples} amostras e os "
                              f"{g.n_symbols} simbolos ocupam {g.n_samples}. O pulso de Nyquist "
                              f"pede um numero inteiro de simbolos na janela.")
        h = PULSES[self.pulse](grid.freq, self.symbol_rate, self.rolloff)
        x = np.array([shape_pulses(map_bits(row, self.modulation), g.samples_per_symbol, h)
                      for row in bits.data])
        return Waveform(samples=x, grid=grid, domain=Domain.ELECTRICAL, geometry=g,
                        channel_id=self.channel_id)
