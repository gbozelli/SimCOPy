"""Lasers do transmissor (docs/laser.md, passo 8).

Dois blocos com a interface de todo bloco, __call__(x, ctx) -> Waveform:
- CWLaser: modelo B, potencia fixa e fase em passeio aleatorio (kernels/laser/cw.py).
- DirectModulatedLaser: modelo A, equacoes de taxa dirigidas pela corrente
  I(t) = I_th*(bias_ratio + modulation_ratio*x(t)) (Agrawal, Eq. 3.3.14, com I_b e I_m em
  multiplos de I_th, decisao de 2026-10-16).

Os dois devolvem o campo optico (2, N): o campo na polarizacao X e Y = 0 (decisao de
2026-10-08). O ruido vem de ctx.rng_for(canal, papel).

`frequency` e a frequencia absoluta do canal (h*nu). O campo sai em banda base do canal, e
`channel_offset` (o desvio do canal em relacao ao centro da banda) vai para Waveform.f0: quem leva
o canal para a sua posicao na grade e o mux (decisao de 2026-11-12, duvida 2).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ...config.system import CWLaserConfig, RateEquationLaserConfig
from ...config.resolve import ResolvedRateEquationLaser, resolve_rate_equation_laser
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, Waveform
from ...kernels.laser.conversions import facet_power
from ...kernels.laser.cw import cw_field
from ...kernels.laser.rate_equations import (RateEquationParams, derivatives, integrate,
                                             steady_state)

def rate_equation_params(cfg: RateEquationLaserConfig,
                         r: ResolvedRateEquationLaser) -> RateEquationParams:
    """Numeros que o kernel recebe: entradas da configuracao + grandezas do resolve."""
    return RateEquationParams(
        gain_coefficient=r.gain_coefficient, transparency_density=cfg.transparency_density,
        photon_lifetime=r.photon_lifetime, carrier_lifetime=cfg.carrier_lifetime,
        volume=r.volume, spontaneous_emission_rate=r.spontaneous_emission_rate,
        gain_compression=cfg.gain_compression,
        linewidth_enhancement=cfg.linewidth_enhancement)

def _rng(x: Waveform, ctx: SimContext, role: str, who: str) -> np.random.Generator:
    if x.channel_id is None:
        raise SignalError(f"{who} precisa de channel_id para sortear o ruido "
                          f"(ctx.rng_for(canal, '{role}')).")
    return ctx.rng_for(x.channel_id, role)

def _optical(x: Waveform, field_x: np.ndarray, channel_offset: float) -> Waveform:
    samples = np.zeros((2, field_x.size), dtype=complex)
    samples[0] = field_x
    return Waveform(samples=samples, grid=x.grid, domain=Domain.OPTICAL,
                    geometry=x.geometry, f0=channel_offset, channel_id=x.channel_id)

@dataclass(frozen=True)
class CWLaser:
    """Laser quase-monocromatico. A entrada so da a grade e o canal: em modulacao
    externa, o modulador combina esta portadora com o sinal eletrico."""
    laser: CWLaserConfig
    frequency: float                  # frequencia absoluta do canal [Hz]
    frequency_offset: float = 0.0     # desvio do campo em relacao ao canal [Hz]; no LO, lo_frequency_offset
    rng_role: str = "tx_laser_phase"
    channel_offset: float = 0.0       # desvio do canal em relacao ao centro da banda -> Waveform.f0

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        rng = _rng(x, ctx, self.rng_role, "CWLaser")
        field_x = cw_field(x.grid.n_samples, x.grid.dt, self.laser.power,
                           self.laser.linewidth, self.frequency_offset, rng)
        return _optical(x, field_x, self.channel_offset)

@dataclass(frozen=True)
class DirectModulatedLaser:
    """Laser de modulacao direta. O sinal eletrico x(t) (adimensional, parte real) vira
    a corrente I(t) = I_th*(bias_ratio + modulation_ratio*x(t)); as equacoes de taxa dao
    S(t) e phi(t), e a saida e sqrt(P_e)*exp(-i*phi_rel), com P_e por faceta (Eq. 3.3.8).

    phi_rel e a fase menos a do estacionario na polarizacao I_b = bias_ratio*I_th: a portadora
    fica na frequencia do canal e o chirp e a variacao em torno dela (decisao de 2026-10-22).
    O sinal: delta_nu = +(1/2pi)*dphi/dt e azul (Agrawal, Eq. 3.4.1), e na convencao do
    Agrawal adotada no projeto (portadora exp(-i*w0*t)) o envelope dessa frequencia gira como
    exp(-i*phi).

    Condicao inicial: estacionario em I(0). Com `noise`, forcas de Langevin e fase
    inicial uniforme, como no CW; sem `noise`, deterministico e com phi(0) = 0.
    """
    laser: RateEquationLaserConfig
    frequency: float                  # frequencia absoluta do canal [Hz]; da h*f
    bias_ratio: float                 # I_b/I_th
    modulation_ratio: float           # I_m/I_th
    noise: bool = True
    max_step: float = 0.1e-12         # subpasso maximo do integrador [s]
    channel_offset: float = 0.0       # desvio do canal em relacao ao centro da banda -> Waveform.f0

    def current(self, x: Waveform) -> np.ndarray:
        """Corrente de acionamento I(t) [A] para o sinal eletrico x."""
        i_th = resolve_rate_equation_laser(self.laser, self.frequency).threshold_current
        return i_th * (self.bias_ratio + self.modulation_ratio * x.samples[0].real)

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.ELECTRICAL, "DirectModulatedLaser")
        r = resolve_rate_equation_laser(self.laser, self.frequency)
        p = rate_equation_params(self.laser, r)
        current = self.current(x)
        n0, s0 = steady_state(current[0], p)
        rng, phi0 = None, 0.0
        if self.noise:
            rng = _rng(x, ctx, "tx_laser_noise", "DirectModulatedLaser")
            phi0 = rng.uniform(0.0, 2 * np.pi)
        _, s, phi = integrate(current, x.grid.dt, p, n0, s0, phi0=phi0,
                              max_step=self.max_step, rng=rng)
        pe = facet_power(s, r.volume, r.group_velocity, self.laser.mirror_loss, r.photon_energy)
        i_b = r.threshold_current * self.bias_ratio
        n_b, s_b = steady_state(i_b, p)
        w_b = derivatives(n_b, s_b, i_b, p)[2]           # dphi/dt no estacionario em I_b
        phi_rel = phi - w_b * x.grid.dt * np.arange(phi.size)
        return _optical(x, np.sqrt(pe) * np.exp(-1j * phi_rel), self.channel_offset)
