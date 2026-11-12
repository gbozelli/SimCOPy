"""Multiplexador e demultiplexador WDM como blocos (docs/mux_demux.md; kernels/wdm.py).

Quem leva cada canal para a sua frequencia e o mux (duvida 2, opcao (b), decisao de 2026-11-12):
os transmissores ficam em banda base do canal, com o desvio do canal em Waveform.f0, e o mux
desloca cada canal para esse desvio e soma. Depois do mux, f0 = 0: o envelope fica referenciado ao
centro da banda. O demux filtra o canal e o deixa na posicao dele (duvida 31, decisao de
2026-11-12): no IMDD o fotodetector ve |E|^2 e a posicao nao importa.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import ClassVar
import numpy as np
from ...config.system import MuxConfig, OpticalFilterConfig
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, SignalGeometry, Waveform
from ...kernels.wdm import apply_optical_filter, bandpass_response, frequency_shift

def _response(x: Waveform, center: float, f: OpticalFilterConfig) -> np.ndarray:
    return bandpass_response(x.grid.optical_frequency, center, f.bandwidth, f.filter_type,
                             f.filter_order or 1)

@dataclass(frozen=True)
class Multiplexer:
    """N campos (n_pol, N) -> campo agregado (2, N). Canal k vai para offsets[k] [Hz], positivo =
    azul, e confere com o f0 de cada entrada.

    Com o padrao de MuxConfig, soma ideal (professor). Com `mux.optical_filter`, cada canal passa
    por um filtro centrado nele antes da soma, e `mux.insertion_loss_db` atenua o campo agregado
    (VPI, WDM_MUX_N_1). A geometria da saida e a do primeiro canal: o campo agregado nao tem uma
    geometria propria (docs/mapa_arquitetura.md, 5.5).
    """
    offsets: tuple[float, ...]
    mux: MuxConfig = field(default_factory=MuxConfig)
    id: str = "mux"
    output_unit: ClassVar[str] = "sqrt(W)"

    def __call__(self, xs: list[Waveform], ctx: SimContext) -> Waveform:
        if len(xs) != len(self.offsets):
            raise SignalError(f"{self.id}: {len(xs)} entradas para {len(self.offsets)} desvios.")
        time = xs[0].grid.time
        total = np.zeros((2, xs[0].grid.n_samples), dtype=complex)
        for k, (x, offset) in enumerate(zip(xs, self.offsets)):
            x.require(Domain.OPTICAL, self.id)
            if abs(x.f0 - offset) > 1e-6 * max(1.0, abs(offset)):
                raise SignalError(f"{self.id}: entrada {k} chegou com f0 = {x.f0 / 1e9:.3f} GHz, e "
                                  f"o canal esta em {offset / 1e9:.3f} GHz na grade.")
            s = x.samples if x.n_pol == 2 else np.vstack([x.samples, np.zeros_like(x.samples)])
            s = frequency_shift(s, offset, time)
            if self.mux.optical_filter is not None:
                s = apply_optical_filter(s, _response(x, offset, self.mux.optical_filter))
            total += s
        total *= 10.0 ** (-self.mux.insertion_loss_db / 20.0)
        return Waveform(samples=total, grid=xs[0].grid, domain=Domain.OPTICAL,
                        geometry=xs[0].geometry, f0=0.0, channel_id=None)

@dataclass(frozen=True)
class Demultiplexer:
    """Filtro optico passa-banda centrado no canal (bandpass_response), aplicado na FFT em cada
    polarizacao. A saida volta a ter o canal (channel_id) e a geometria dele; o canal fica na sua
    posicao, entao f0 continua 0."""
    filter: OpticalFilterConfig
    center: float                     # desvio do canal [Hz], positivo = azul
    channel_id: int | None = None
    geometry: SignalGeometry | None = None
    id: str = "rx.demux"
    output_unit: ClassVar[str] = "sqrt(W)"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, self.id)
        return Waveform(samples=apply_optical_filter(x.samples, _response(x, self.center, self.filter)),
                        grid=x.grid, domain=Domain.OPTICAL, geometry=self.geometry or x.geometry,
                        f0=x.f0, channel_id=self.channel_id)
