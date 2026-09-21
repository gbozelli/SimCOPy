"""Resolve grandezas derivadas. O que esta aqui deixa de ser configuravel: e
calculado, registrado no resultado, e nunca aceito como entrada (secao 4.5.2).

No codigo antigo, {BaudRate, SpS, ts, fs} tinha um grau de liberdade mas cinco
funcoes recebiam SpS E ts como parametros independentes, sem nada verificar que
ts*SpS == 1/BaudRate.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..core.grid import TimeGrid
from ..core.signals import SignalGeometry
from ..core.units import (beta2_from_dispersion, gamma_from_n2, effective_length,
                          frequency_to_wavelength, db_per_m_to_np_per_m)
from .system import SystemConfig, FiberConfig
from .simulation import SimulationConfig

@dataclass(frozen=True)
class ResolvedChannel:
    index: int
    f0: float
    samples_per_symbol: int
    n_symbols: int
    geometry: SignalGeometry

@dataclass(frozen=True)
class ResolvedFiber:
    beta2: float
    gamma: float
    alpha_np: float
    effective_length: float
    dispersion_length: float | None

@dataclass(frozen=True)
class Resolved:
    grid: TimeGrid
    channels: tuple[ResolvedChannel, ...]
    fibers: tuple[ResolvedFiber, ...]
    n_samples: int
    window: float

def _next_pow2(n: int) -> int:
    return 1 << (int(n) - 1).bit_length()

def resolve_fiber(f: FiberConfig, center_frequency: float,
                  pulse_width: float | None = None) -> ResolvedFiber:
    lam = frequency_to_wavelength(center_frequency)
    b2 = beta2_from_dispersion(f.dispersion, lam)
    return ResolvedFiber(
        beta2=b2,
        gamma=gamma_from_n2(f.n2, f.effective_area, lam),
        alpha_np=db_per_m_to_np_per_m(f.attenuation),
        effective_length=effective_length(f.length, f.attenuation),
        dispersion_length=(pulse_width**2 / abs(b2)) if pulse_width else None,
    )

def resolve(system: SystemConfig, simulation: SimulationConfig) -> Resolved:
    fs = simulation.sampling_rate
    n = int(round(simulation.simulation_window * fs))
    if simulation.force_power_of_two:            # GreatestPrimeFactorLimit = 2
        n = _next_pow2(n)
    grid = TimeGrid(sampling_rate=fs, n_samples=n)

    offs = system.channel_offsets()
    chans = []
    for c, f0 in zip(system.channels, offs):
        sps_f = fs / c.symbol_rate
        sps = int(round(sps_f))                  # [D-4]: exigido inteiro
        n_sym = n // sps                         # truncado
        chans.append(ResolvedChannel(
            index=c.index, f0=f0, samples_per_symbol=sps, n_symbols=n_sym,
            geometry=SignalGeometry(n_symbols=n_sym,
                                    bits_per_symbol=c.modulation.bits_per_symbol,
                                    samples_per_symbol=sps,
                                    n_pol=c.pol.n_pol)))
    fibers = tuple(resolve_fiber(s.fiber, system.center_frequency)
                   for s in system.link.spans)
    return Resolved(grid=grid, channels=tuple(chans), fibers=fibers,
                    n_samples=n, window=grid.duration)
