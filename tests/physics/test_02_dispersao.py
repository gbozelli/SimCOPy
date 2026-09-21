"""2. DISPERSAO CROMATICA — alargamento de pulso gaussiano (Agrawal, NLFO, sec. 3.2)

Sem chirp (C = 0):
    T1 / T0 = sqrt(1 + (z / L_D)^2),     L_D = T0^2 / |beta2|
e como a FWHM de uma gaussiana e 2*sqrt(ln 2)*T0, a mesma razao vale para a FWHM.

Com chirp (C != 0):
    T1 / T0 = sqrt((1 + C*beta2*z/T0^2)^2 + (beta2*z/T0^2)^2)

Por que o caso com chirp importa: o caso sem chirp depende so de |beta2|, entao
um ERRO DE SINAL na dispersao passaria despercebido. Com C*beta2 < 0 o pulso
primeiro COMPRIME e depois alarga; com C*beta2 > 0 so alarga. So o sinal certo de
beta2 reproduz os dois. E o teste que pega convencao de dispersao trocada, a
mesma armadilha da comparacao com o VPI (secao 14.3 do documento).
"""
import numpy as np
import pytest
from tests._approx import approx
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import beta2_from_dispersion
from ._helpers import gaussian, fwhm

N, FS = 1 << 14, 1e12
T = (np.arange(N) - N // 2) / FS
OMEGA = 2 * np.pi * np.fft.fftfreq(N, 1 / FS)
T0 = 20e-12
B2 = beta2_from_dispersion(16e-6, 1550e-9)      # < 0: dispersao anomala
LD = T0**2 / abs(B2)

def _fwhm_out(C, z):
    A = gaussian(T, T0, C=C)[None, :]
    out = propagate(A, OMEGA, FiberKernelParams(z, B2, 0.0, 0.0), max_step=z / 200)
    return fwhm(T, np.abs(out[0]) ** 2)

@pytest.mark.parametrize("z_LD", [0.5, 1.0, 2.0, 4.0])
def test_fwhm_sem_chirp(z_LD):
    z = z_LD * LD
    esperado = fwhm(T, np.abs(gaussian(T, T0)) ** 2) * np.sqrt(1 + (z / LD) ** 2)
    assert _fwhm_out(0.0, z) == approx(esperado, rel=2e-3)

@pytest.mark.parametrize("C", [-2.0, +2.0])
@pytest.mark.parametrize("z_LD", [0.25, 0.5, 1.0])
def test_fwhm_com_chirp_verifica_sinal_de_beta2(C, z_LD):
    z = z_LD * LD
    fator = np.sqrt((1 + C * B2 * z / T0**2) ** 2 + (B2 * z / T0**2) ** 2)
    esperado = fwhm(T, np.abs(gaussian(T, T0)) ** 2) * fator
    assert _fwhm_out(C, z) == approx(esperado, rel=2e-3)

def test_compressao_inicial_quando_C_beta2_negativo():
    """Com beta2 < 0, C > 0 da C*beta2 < 0: o pulso tem que ficar MAIS ESTREITO."""
    z = 0.25 * LD
    assert _fwhm_out(+2.0, z) < fwhm(T, np.abs(gaussian(T, T0)) ** 2)
