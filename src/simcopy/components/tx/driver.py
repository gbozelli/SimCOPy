"""Acionamento eletrico do transmissor IMDD: bits -> NRZ -> filtro (docs/sistema_imdd.md, secao 2)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
from ...config.system import ElectricalFilterConfig
from ...core.context import SimContext
from ...core.signals import Bits, Domain, SignalError, Waveform
from ...kernels.electrical import nrz
from ...kernels.filters import apply_response, lowpass_response

@dataclass(frozen=True)
class NRZDriver:
    """Bits OOK -> sinal eletrico em [0, 1] na grade da simulacao (notebook 07, primeira celula).

    Com `filter`, passa-baixas de -3 dB em filter.bandwidth_factor * symbol_rate, aplicado na FFT
    (kernels/filters.py). Sem `filter`, NRZ ideal.
    """
    symbol_rate: float
    filter: ElectricalFilterConfig | None = None
    channel_id: int | None = None
    id: str = "tx.driver"
    output_unit: ClassVar[str] = "adimensional"

    def __call__(self, bits: Bits, ctx: SimContext) -> Waveform:
        g = bits.geometry
        if g.bits_per_symbol != 1 or g.n_pol != 1:
            raise SignalError(f"{self.id}: so OOK de uma polarizacao (1 bit por simbolo).")
        grid = ctx.grid
        x = nrz(bits.data[0], g.samples_per_symbol, grid.n_samples)
        if self.filter is not None:
            f = self.filter
            x = apply_response(x, lowpass_response(grid.freq, f.bandwidth_factor * self.symbol_rate,
                                                   f.filter_type, f.filter_order))
        return Waveform(samples=x.astype(complex)[None, :], grid=grid, domain=Domain.ELECTRICAL,
                        geometry=g, channel_id=self.channel_id)
