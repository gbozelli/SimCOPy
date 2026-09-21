"""Birefringencia e PMD, modelo de passo grosso (secao 10.4).

Corrige dois defeitos do codigo antigo: o sorteio usava 2 angulos quando SU(2)
tem 3 (o estado ficava preso no plano S3=0 partindo de X linear), e a rotacao era
re-sorteada a cada passo DeltaL do SSFM, tornando a fisica dependente da
discretizacao numerica.

O sorteio uniforme segue o metodo de quaternios de BIFROST (Banner et al.,
Phys. Rev. Applied, arXiv:2510.01212, eq. 17): quatro gaussianas normalizadas
para um ponto na 3-esfera.
"""
from __future__ import annotations
import numpy as np

_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)

def random_su2(rng: np.random.Generator) -> np.ndarray:
    """Matriz de Jones uniformemente distribuida em SU(2) (medida de Haar)."""
    g = rng.standard_normal(4)
    g /= np.linalg.norm(g)
    ct, a = g[0], g[1:]
    na = np.linalg.norm(a)
    if na == 0.0:
        return np.eye(2, dtype=complex)
    st = np.sqrt(max(0.0, 1.0 - ct * ct))
    a = a / na
    return (np.eye(2, dtype=complex) * ct
            - 1j * st * (a[0] * _SX + a[1] * _SY + a[2] * _SZ))

def sections_crossed(z0: float, z1: float, correlation_length: float) -> int:
    """Numero de fronteiras de secao cruzadas no intervalo [z0, z1)."""
    if correlation_length <= 0:
        return 1
    n = int(np.floor(z1 / correlation_length)) - int(np.floor(z0 / correlation_length))
    return max(0, n) + (1 if z0 == 0.0 else 0)

def dgd_for_step(pmd_coefficient: float, step: float,
                 correlation_length: float) -> float:
    """DGD de um passo do SSFM.

    Dentro de uma secao de birrefringencia (comprimento de correlacao L_c) a
    orientacao nao muda, entao o DGD cresce LINEARMENTE com a distancia. Entre
    secoes a orientacao e sorteada de novo, e o DGD acumula em sqrt(L).

      passo <  L_c : dtau = D_PMD * passo / sqrt(L_c)   (linear dentro da secao)
      passo >= L_c : dtau = D_PMD * sqrt(passo)         (passeio aleatorio)

    Versao anterior usava sempre D_PMD*sqrt(passo): com passo de 5 m e L_c de
    50 m o DGD saia 3,18x maior, ou seja, dependia do passo numerico.
    """
    return pmd_coefficient * step / np.sqrt(max(step, correlation_length))

def apply_dgd_freq(Ax_f, Ay_f, omega, dtau):
    """Elemento de DGD de 1a ordem, aplicado no dominio da frequencia.
    Roda dentro do passo linear do SSFM, entao nao custa FFT extra."""
    if dtau == 0.0:
        return Ax_f, Ay_f
    ph = np.exp(-0.5j * omega * dtau)
    return Ax_f * ph, Ay_f * np.conj(ph)

def apply_rotation(Ax, Ay, U):
    return U[0, 0] * Ax + U[0, 1] * Ay, U[1, 0] * Ax + U[1, 1] * Ay
