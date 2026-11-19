"""Acionamento eletrico do transmissor IMDD: bits -> NRZ -> filtro (docs/sistema_imdd.md, secao 2)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
from ...config.system import ElectricalFilterConfig
from ...core.context import SimContext
from ...core.signals import Bits, Domain, SignalError, Waveform
from ...kernels.electrical import nrz, pam_drive_levels
from ...kernels.filters import apply_response, lowpass_response

@dataclass(frozen=True)
class NRZDriver:
    """Bits -> sinal eletrico em [0, 1] na grade da simulacao (notebook 07, primeira celula). OOK
    (1 bit por simbolo: niveis 0 e 1) ou PAM4 (2 bits por simbolo: niveis 0, 1/3, 2/3 e 1, com
    codigo de Gray; kernels/electrical.py, `pam_drive_levels`).

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
        if g.bits_per_symbol not in (1, 2) or g.n_pol != 1:
            raise SignalError(f"{self.id}: so OOK ou PAM4 de uma polarizacao (1 ou 2 bits por "
                              f"simbolo); a geometria tem {g.bits_per_symbol} e n_pol = {g.n_pol}.")
        grid = ctx.grid
        x = nrz(pam_drive_levels(bits.data[0], g.bits_per_symbol), g.samples_per_symbol,
                grid.n_samples)
        if self.filter is not None:
            f = self.filter
            x = apply_response(x, lowpass_response(grid.freq, f.bandwidth_factor * self.symbol_rate,
                                                   f.filter_type, f.filter_order))
        return Waveform(samples=x.astype(complex)[None, :], grid=grid, domain=Domain.ELECTRICAL,
                        geometry=g, channel_id=self.channel_id)
