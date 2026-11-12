"""Fotodetector p-i-n como bloco (docs/fotodetectores.md, secao 5; notebook 07, funcao detect)."""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
from ...config.system import PhotodiodeConfig
from ...core.context import SimContext
from ...core.signals import Domain, Waveform
from ...kernels.filters import apply_response, lowpass_response
from ...kernels.photodetector import photocurrent, receiver_noise

@dataclass(frozen=True)
class Photodiode:
    """Campo optico (n_pol, N) -> corrente [A] num Waveform eletrico (1, N) com parte imaginaria
    nula (D-5): R*(|A_x|^2 + |A_y|^2) + I_d, mais ruido de disparo e termico, filtrado pelo
    passa-baixas de banda bandwidth_factor * symbol_rate (C-22). Ruido de ctx.rng_for(id, "noise").
    """
    photodiode: PhotodiodeConfig
    symbol_rate: float
    id: str = "rx.photodiode"
    output_unit: ClassVar[str] = "A"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, self.id)
        c = self.photodiode
        grid = x.grid
        i = photocurrent(x.samples, c.responsivity, c.dark_current)
        noise = receiver_noise(i, grid.sampling_rate, ctx.rng_for(self.id, "noise"),
                               thermal=c.thermal_noise_density, shot=c.shot_noise)
        h = lowpass_response(grid.freq, c.bandwidth_factor * self.symbol_rate,
                             c.filter_type, c.filter_order)
        out = apply_response(i + noise, h)
        return Waveform(samples=out.astype(complex)[None, :], grid=grid, domain=Domain.ELECTRICAL,
                        geometry=x.geometry, channel_id=x.channel_id)
