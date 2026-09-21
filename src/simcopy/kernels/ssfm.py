"""Split-step de Fourier simetrico, escalar e vetorial (Manakov).

Diferencas em relacao ao codigo antigo:
  C-8  simetrico (meio passo linear / passo nao linear / meio passo linear),
       erro local O(dz^3) em vez de O(dz^2)
  C-8  passo adaptativo por mudanca de fase nao linear, como o VPI
       (StepSelectionMethod = NonlinearPhaseChange, MaxPhaseChange = 0.5 deg)
  C-6  fator 8/9 de Manakov, confirmado pelo arquivo VPI
       (NonlinearAdjustmentFactors = (8./9.) ...)
  C-7  birrefringencia de passo grosso, desacoplada de dz
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from . import birefringence as bir

MANAKOV_FACTOR = 8.0 / 9.0

@dataclass(frozen=True)
class FiberKernelParams:
    length: float
    beta2: float
    gamma: float
    alpha_np: float
    pmd_coefficient: float = 0.0
    correlation_length: float = 50.0
    manakov: bool = True

def _linear_half(Ax, Ay, omega, beta2, alpha_np, dz, dtau=0.0):
    op = np.exp((-0.5 * alpha_np + 0.5j * beta2 * omega**2) * (dz / 2.0))
    Ax_f = np.fft.fft(Ax) * op
    if Ay is None:
        return np.fft.ifft(Ax_f), None
    Ay_f = np.fft.fft(Ay) * op
    Ax_f, Ay_f = bir.apply_dgd_freq(Ax_f, Ay_f, omega, dtau)
    return np.fft.ifft(Ax_f), np.fft.ifft(Ay_f)

def propagate(A, omega, p: FiberKernelParams, rng=None,
              max_step=5e3, max_phase_deg=0.5, min_step=1.0):
    """A: (n_pol, N) complexo. Devolve o campo apos p.length metros."""
    A = np.array(A, dtype=complex, copy=True)
    vec = A.shape[0] == 2
    Ax = A[0]
    Ay = A[1] if vec else None
    nl_factor = MANAKOV_FACTOR if p.manakov else 1.0
    max_phase = np.deg2rad(max_phase_deg)
    z = 0.0
    while z < p.length - 1e-9:
        inten = np.abs(Ax)**2 + (np.abs(Ay)**2 if vec else 0.0)
        pmax = float(np.max(inten)) * nl_factor * p.gamma
        dz = max_step if pmax <= 0 else min(max_step, max_phase / pmax)
        dz = max(min_step, min(dz, p.length - z))

        dtau = bir.dgd_for_step(p.pmd_coefficient, dz, p.correlation_length) if vec else 0.0
        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz, dtau)

        # Comprimento efetivo do passo. A intensidade aqui e a do ponto medio
        # (apos meio passo linear); integrar exp(-alpha*z) exatamente no passo da
        # 2*sinh(alpha*dz/2)/alpha em vez de dz. Com beta2 = 0 isso torna a SPM
        # exata; com dispersao, mantem 2a ordem com constante de erro menor.
        # Sem isso, passos de 5 km erravam a fase nao linear em ~5e-4.
        dz_nl = dz if p.alpha_np == 0 else 2.0 * np.sinh(p.alpha_np * dz / 2) / p.alpha_np
        inten = np.abs(Ax)**2 + (np.abs(Ay)**2 if vec else 0.0)
        nl = np.exp(1j * p.gamma * nl_factor * inten * dz_nl)
        Ax = Ax * nl
        if vec:
            Ay = Ay * nl

        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz)

        if vec and rng is not None and p.pmd_coefficient >= 0:
            if bir.sections_crossed(z, z + dz, p.correlation_length) > 0:
                Ax, Ay = bir.apply_rotation(Ax, Ay, bir.random_su2(rng))
        z += dz
    return np.vstack([Ax, Ay]) if vec else Ax[None, :]
