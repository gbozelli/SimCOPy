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
                          frequency_to_wavelength, db_per_m_to_np_per_m,
                          C_LIGHT, H_PLANCK, Q_ELECTRON)
from .system import SystemConfig, FiberConfig, RateEquationLaserConfig
from .simulation import SimulationConfig

@dataclass(frozen=True)
class ResolvedChannel:
    index: int
    f0: float
    samples_per_symbol: int
    n_symbols: int
    geometry: SignalGeometry
    edge_guard_symbols: int           # simbolos descartados em cada borda da janela na medida

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

@dataclass(frozen=True)
class ResolvedRateEquationLaser:
    volume: float                     # V = L*W*H [m^3]
    group_velocity: float             # v_g [m/s]
    photon_lifetime: float            # tau_p [s]
    gain_coefficient: float           # Gamma*v_g*sigma_g [m^3/s]
    threshold_density: float          # N_th [1/m^3]
    threshold_current: float          # I_th [A]
    spontaneous_emission_rate: float  # R_sp [1/(m^3 s)]
    photon_energy: float              # h*f [J]
    differential_efficiency: float    # eta_d
    slope_efficiency: float           # dPe/dI por faceta [W/A]

def resolve_rate_equation_laser(cfg: RateEquationLaserConfig,
                                frequency: float) -> ResolvedRateEquationLaser:
    """Grandezas derivadas do laser (docs/laser.md, secoes 3.2 e 3.6).

    A frequencia optica vem do canal, para que lambda tenha um dono so.
    Agrawal: Eq. 3.3.4 (tau_p), 3.3.6 (limiar), 3.3.10 (eficiencia, eta_int = 1).
    """
    volume = cfg.cavity_length * cfg.cavity_width * cfg.cavity_height
    vg = C_LIGHT / cfg.group_index
    tau_p = 1.0 / (vg * (cfg.mirror_loss + cfg.internal_loss))
    gc = cfg.confinement * vg * cfg.differential_gain
    n_th = cfg.transparency_density + 1.0 / (gc * tau_p)
    i_th = Q_ELECTRON * n_th * volume / cfg.carrier_lifetime
    # Agrawal: R_sp = n_sp*G (texto apos a Eq. 3.3.3), avaliado no limiar, G = 1/tau_p,
    # como nos Ex. 3.10 e 3.14 (R_sp = 2/tau_p); dividido por V para a forma em
    # densidade. Decisao de 2026-10-16 (docs/laser.md, secao 3.5).
    r_sp = cfg.spontaneous_emission_factor / (tau_p * volume)
    hf = H_PLANCK * frequency
    eta_d = cfg.mirror_loss / (cfg.mirror_loss + cfg.internal_loss)
    return ResolvedRateEquationLaser(
        volume=volume, group_velocity=vg, photon_lifetime=tau_p, gain_coefficient=gc,
        threshold_density=n_th, threshold_current=i_th, spontaneous_emission_rate=r_sp,
        photon_energy=hf, differential_efficiency=eta_d,
        slope_efficiency=hf / (2.0 * Q_ELECTRON) * eta_d)

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

EDGE_SETTLING_TIME = 1e-9     # s, memoria do transmissor e dos filtros na borda da janela

def edge_guard_symbols(fibers, lengths, sampling_rate: float, symbol_rate: float) -> int:
    """Simbolos de cada borda da janela que a medida descarta (docs/sistema_imdd.md, secao 3).

    A FFT trata a janela como periodica, e o inicio e o fim do sinal nao se emendam (fase do
    laser, primeiro e ultimo bits). Na fibra, a componente na borda da banda de simulacao
    (omega = pi*fs) se atrasa pi*|beta2|*L*fs em relacao a portadora: e ate onde a emenda pode
    contaminar a janela (docs/mapa_arquitetura.md, 5.7). EDGE_SETTLING_TIME cobre a memoria do
    transmissor: o laser de modulacao direta parte do estacionario de I(0), sem a historia dos bits
    anteriores, e difere da janela longa por 0,6 ns (6 bits a 10 GBd, medido em 2026-11-12).
    """
    spread = np.pi * sum(abs(f.beta2) * L for f, L in zip(fibers, lengths)) * sampling_rate
    return int(np.ceil((spread + EDGE_SETTLING_TIME) * symbol_rate)) + 1

def resolve(system: SystemConfig, simulation: SimulationConfig) -> Resolved:
    fs = simulation.sampling_rate
    n = int(round(simulation.simulation_window * fs))
    if simulation.force_power_of_two:            # GreatestPrimeFactorLimit = 2
        n = _next_pow2(n)
    grid = TimeGrid(sampling_rate=fs, n_samples=n)

    # tem q fazer uns bons testes aqui p ver se os spans n vão conflitar
    fibers = tuple(resolve_fiber(s.fiber, system.center_frequency)
                   for s in system.link.spans)
    lengths = [s.fiber.length for s in system.link.spans]

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
                                    n_pol=c.pol.n_pol),
            edge_guard_symbols=edge_guard_symbols(fibers, lengths, fs, c.symbol_rate)))
    return Resolved(grid=grid, channels=tuple(chans), fibers=fibers,
                    n_samples=n, window=grid.duration)
