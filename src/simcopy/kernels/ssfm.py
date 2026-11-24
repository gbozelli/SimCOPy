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
    beta3: float = 0.0                # s^3/m (Agrawal, Eq. 2.4.4)

def _linear_half(Ax, Ay, omega, beta2, alpha_np, dz, dtau=0.0, beta3=0.0):
    phase = 0.5j * beta2 * omega**2
    if beta3 != 0.0:
        # Agrawal, Eq. 2.4.4: beta(w) tem o termo (beta3/6)*(dw)^3, com dw a frequencia optica
        # em relacao a portadora. Na convencao do projeto (portadora exp(-i*w0*t)) o eixo da
        # FFT e -dw: o termo cubico troca de sinal, o quadratico nao.
        phase = phase - (1j / 6.0) * beta3 * omega**3
    op = np.exp((-0.5 * alpha_np + phase) * (dz / 2.0))
    Ax_f = np.fft.fft(Ax) * op
    if Ay is None:
        return np.fft.ifft(Ax_f), None
    Ay_f = np.fft.fft(Ay) * op
    Ax_f, Ay_f = bir.apply_dgd_freq(Ax_f, Ay_f, omega, dtau)
    return np.fft.ifft(Ax_f), np.fft.ifft(Ay_f)

def propagate(A, omega, p: FiberKernelParams, rng=None,
              max_step=5e3, max_phase_deg=0.5, min_step=1.0, positions=None):
    """A: (n_pol, N) complexo. Devolve o campo apos p.length metros.

    positions: distancias [m] em ordem crescente, de 0 a p.length, em que o campo tambem e
    guardado (duvida 32, resposta de 2026-11-16). Com elas a funcao devolve (campo final, lista
    com o campo em cada posicao), no mesmo laco, sem propagar de novo. O passo e encurtado para
    cair em cada posicao, entao a sequencia de passos (e o sorteio da PMD) nao e a mesma de uma
    chamada sem `positions`.
    """
    A = np.array(A, dtype=complex, copy=True)
    vec = A.shape[0] == 2
    Ax = A[0]
    Ay = A[1] if vec else None
    pending, snapshots = [], []
    if positions is not None:
        pending = [float(v) for v in positions]
        if any(b < a for a, b in zip(pending, pending[1:])) or \
                (pending and (pending[0] < 0.0 or pending[-1] > p.length + 1e-9)):
            raise ValueError("positions: use distancias crescentes entre 0 e o comprimento da fibra")
        while pending and pending[0] <= 1e-9:          # posicao 0: o campo de entrada
            snapshots.append(A.copy())
            pending.pop(0)
    nl_factor = MANAKOV_FACTOR if p.manakov else 1.0
    max_phase = np.deg2rad(max_phase_deg)
    z = 0.0
    while z < p.length - 1e-9:
        inten = np.abs(Ax)**2 + (np.abs(Ay)**2 if vec else 0.0)
        pmax = float(np.max(inten)) * nl_factor * p.gamma
        dz = max_step if pmax <= 0 else min(max_step, max_phase / pmax)
        # O ultimo passo nunca pode ultrapassar o fim da fibra: min_step so vale
        # enquanto couber. Antes, max(min_step, ...) fazia a fibra ficar mais
        # longa do que a pedida quando L nao era multiplo do passo.
        dz = min(max(dz, min_step), p.length - z)
        if pending:
            dz = min(dz, pending[0] - z)

        dtau = bir.dgd_for_step(p.pmd_coefficient, dz, p.correlation_length) if vec else 0.0
        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz, dtau, p.beta3)

        # Comprimento efetivo do passo. A intensidade aqui e a do ponto medio
        # (apos meio passo linear), integrar exp(-alpha*z) exatamente no passo da
        # 2*sinh(alpha*dz/2)/alpha em vez de dz. Com beta2 = 0 isso torna a SPM
        # exata; com dispersao, mantem 2a ordem com constante de erro menor.
        # Sem isso, passos de 5 km erravam a fase nao linear em ~5e-4.
        dz_nl = dz if p.alpha_np == 0 else 2.0 * np.sinh(p.alpha_np * dz / 2) / p.alpha_np
        inten = np.abs(Ax)**2 + (np.abs(Ay)**2 if vec else 0.0)
        nl = np.exp(1j * p.gamma * nl_factor * inten * dz_nl)
        Ax = Ax * nl
        if vec:
            Ay = Ay * nl

        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz, 0.0, p.beta3)

        if vec and rng is not None and p.pmd_coefficient >= 0:
            if bir.sections_crossed(z, z + dz, p.correlation_length) > 0:
                Ax, Ay = bir.apply_rotation(Ax, Ay, bir.random_su2(rng))
        z += dz
        while pending and z >= pending[0] - 1e-9:
            snapshots.append(np.vstack([Ax, Ay]) if vec else Ax[None, :].copy())
            pending.pop(0)
    out = np.vstack([Ax, Ay]) if vec else Ax[None, :]
    return (out, snapshots) if positions is not None else out
