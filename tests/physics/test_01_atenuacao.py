"""1. ATENUACAO — P(z) = P(0) * exp(-alpha * z)

O que se testa: so a perda. Sem dispersao (beta2 = 0) e sem nao linearidade
(gamma = 0), a potencia media tem que cair exatamente como exp(-alpha*z), com
alpha em Np/m (alpha[dB/m] / 4,343).

Por que a tolerancia e tao apertada (1e-12): nesse caso o SSFM aplica so um
fator multiplicativo por passo; qualquer desvio acima do arredondamento de
ponto flutuante e erro de conta, nao aproximacao numerica.
"""
import numpy as np
import pytest
from tests._approx import approx
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import db_per_m_to_np_per_m

@pytest.mark.parametrize("L_km", [1, 25, 80, 150])
def test_potencia_cai_como_exponencial(L_km):
    A = np.ones((1, 256), dtype=complex)
    omega = 2 * np.pi * np.fft.fftfreq(256, 1e-12)
    alpha = db_per_m_to_np_per_m(0.2e-3)                   # 0,2 dB/km
    out = propagate(A, omega, FiberKernelParams(L_km * 1e3, 0.0, 0.0, alpha))
    P = np.mean(np.abs(out) ** 2)
    assert P == approx(np.exp(-alpha * L_km * 1e3), rel=1e-12)

def test_perda_em_db_e_linear_na_distancia():
    """Consequencia direta: 0,2 dB/km vezes 125 km = 25 dB."""
    A = np.ones((1, 64), dtype=complex)
    omega = 2 * np.pi * np.fft.fftfreq(64, 1e-12)
    out = propagate(A, omega, FiberKernelParams(125e3, 0.0, 0.0,
                                                db_per_m_to_np_per_m(0.2e-3)))
    assert -10 * np.log10(np.mean(np.abs(out) ** 2)) == approx(25.0, abs=1e-9)

@pytest.mark.parametrize("L,passo", [(1234.0, 500.0), (99.0, 500.0), (150e3, 7e3)])
def test_comprimento_nao_multiplo_do_passo(L, passo):
    """O ultimo passo nao pode ultrapassar o fim da fibra. Com L = 1234 m e passo
    de 500 m, propagar 1500 m daria 0,05 dB a mais de perda."""
    A = np.ones((1, 64), dtype=complex)
    omega = 2 * np.pi * np.fft.fftfreq(64, 1e-12)
    alpha = db_per_m_to_np_per_m(0.2e-3)
    out = propagate(A, omega, FiberKernelParams(L, 0.0, 0.0, alpha),
                    max_step=passo, min_step=passo)
    assert np.mean(np.abs(out) ** 2) == approx(np.exp(-alpha * L), rel=1e-12)
