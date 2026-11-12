"""Fibra como bloco, no molde dos lasers (docs/proposta_bloco_fibra.md; docs/mapa_arquitetura.md, 4.5).

A configuracao e a frequencia central vem no construtor, como em CWLaser(laser, frequency); o
bloco calcula os proprios derivados com a mesma funcao do resolve do sistema (resolve_fiber), e
por isso o validador e o bloco veem o mesmo beta2 e o mesmo gamma.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
from typing import ClassVar
from ...config.resolve import ResolvedFiber, resolve_fiber
from ...config.system import FiberConfig
from ...core.context import SimContext
from ...core.signals import Domain, Waveform
from ...kernels.ssfm import FiberKernelParams, propagate

def fiber_kernel_params(cfg: FiberConfig, r: ResolvedFiber) -> FiberKernelParams:
    """Numeros que o kernel recebe: entradas da configuracao + grandezas do resolve. E o unico
    lugar em que a PMD da configuracao passa para o kernel."""
    b = cfg.birefringence
    return FiberKernelParams(
        length=cfg.length, beta2=r.beta2, gamma=r.gamma, alpha_np=r.alpha_np,
        pmd_coefficient=b.pmd_coefficient if b.enabled else 0.0,
        correlation_length=b.correlation_length)

@dataclass(frozen=True)
class Fiber:
    """Propagacao pelo SSFM (kernels/ssfm.py). A saida tem a mesma grade, geometria, canal e f0.

    Sorteio da PMD de ctx.rng_for(id, "pmd"), so com campo (2, N) e coeficiente > 0 (com rng e
    coeficiente zero o kernel embaralharia a polarizacao: proposta da fibra, ponto 2).
    O tempo e o do referencial que anda com o pulso: a saida nao e atrasada por L/v_g.
    """
    fiber: FiberConfig
    frequency: float                  # frequencia central da banda [Hz] (beta2, gamma)
    id: str = "fiber"
    max_step: float = 5e3             # passo maximo do SSFM [m] (SimulationConfig.ssfm_max_step)
    max_phase_deg: float = 0.5        # SimulationConfig.ssfm_max_phase_change_deg
    output_unit: ClassVar[str] = "sqrt(W)"

    @cached_property
    def params(self) -> FiberKernelParams:
        return fiber_kernel_params(self.fiber, resolve_fiber(self.fiber, self.frequency))

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, self.id)
        p = self.params
        rng = ctx.rng_for(self.id, "pmd") if x.n_pol == 2 and p.pmd_coefficient > 0 else None
        return x.with_samples(propagate(x.samples, x.grid.omega, p, rng,
                                        max_step=self.max_step, max_phase_deg=self.max_phase_deg))
