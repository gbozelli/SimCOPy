"""3. NAO LINEARIDADE — automodulacao de fase, SPM (Agrawal, NLFO, sec. 4.1)

Isolando a nao linearidade (beta2 = 0), a equacao tem solucao fechada:

    |A(L,T)|^2 = |A(0,T)|^2 * exp(-alpha*L)        (a forma do pulso NAO muda)
    phi_NL(L,T) = gamma * |A(0,T)|^2 * L_eff        (so a FASE muda)
    L_eff = (1 - exp(-alpha*L)) / alpha

Consequencias testadas aqui, uma por teste:
  a) a intensidade normalizada na saida e identica a da entrada;
  b) a fase nao linear tem a MESMA FORMA do pulso: maxima onde a amplitude e
     maxima, e proporcional a |A(0,T)|^2 em todo instante;
  c) a defasagem maxima vale phi_max = gamma * P0 * L_eff;
  d) em DP (Manakov), phi_max = (8/9) * gamma * P0 * L_eff, com P0 = |Ax|^2+|Ay|^2
     de pico. O 8/9 e o NonlinearAdjustmentFactors do arquivo VPI.

Por que isso responde "como testar nao linearidade": a SPM pura tem solucao
analitica, e ela envolve exatamente gamma e L_eff.
"""
import numpy as np
import pytest
from tests._approx import approx
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import db_per_m_to_np_per_m, effective_length
from ._helpers import gaussian

N, FS = 4096, 1e12
T = (np.arange(N) - N // 2) / FS
OMEGA = 2 * np.pi * np.fft.fftfreq(N, 1 / FS)
T0, L = 20e-12, 50e3
GAMMA = 1.3e-3
ALPHA = db_per_m_to_np_per_m(0.2e-3)
LEFF = effective_length(L, 0.2e-3)
P0 = 2.0 / (GAMMA * LEFF)                     # escolhido para phi_max = 2 rad

def _spm(A, manakov):
    p = FiberKernelParams(L, 0.0, GAMMA, ALPHA, manakov=manakov)
    return propagate(A, OMEGA, p, rng=None)

def _fase_nl(entrada, saida, limiar=1e-6):
    """Fase adicionada pela fibra, so onde ha sinal. Nas caudas do pulso a
    amplitude e ruido de arredondamento das FFTs e a fase e aleatoria; incluir
    essa regiao no unwrap acumula saltos espurios."""
    I0 = np.abs(entrada) ** 2
    m = I0 > limiar * I0.max()
    phi = np.zeros_like(I0)
    phi[m] = np.unwrap(np.angle(saida[m] * np.conj(entrada[m])))
    return phi

def test_a_forma_do_pulso_nao_muda():
    A = gaussian(T, T0, P0)[None, :]
    out = _spm(A, manakov=False)[0]
    I0, I1 = np.abs(A[0]) ** 2, np.abs(out) ** 2
    np.testing.assert_allclose(I1 / I1.max(), I0 / I0.max(), atol=1e-9)

def test_b_fase_tem_a_mesma_forma_do_pulso():
    A = gaussian(T, T0, P0)[None, :]
    phi = _fase_nl(A[0], _spm(A, manakov=False)[0])
    I0 = np.abs(A[0]) ** 2
    assert np.argmax(phi) == np.argmax(I0)                 # maxima no mesmo instante
    janela = I0 > 1e-3 * I0.max()                          # onde ha sinal
    np.testing.assert_allclose(phi[janela] / phi.max(), I0[janela] / I0.max(),
                               atol=2e-4)

def test_c_defasagem_maxima_escalar():
    A = gaussian(T, T0, P0)[None, :]
    phi = _fase_nl(A[0], _spm(A, manakov=False)[0])
    assert phi.max() == approx(GAMMA * P0 * LEFF, rel=1e-9)

def test_d_defasagem_maxima_manakov_dp():
    A = np.vstack([gaussian(T, T0, P0 / 2), gaussian(T, T0, P0 / 2)])
    out = _spm(A, manakov=True)
    phi = _fase_nl(A[0], out[0])
    # 8/9 escrito LITERALMENTE, de proposito. Se o gabarito fosse a constante do
    # proprio codigo (MANAKOV_FACTOR), mudar o codigo mudaria o gabarito junto e o
    # teste nunca falharia. O valor vem da fisica (Manakov) e do arquivo VPI.
    assert phi.max() == approx((8.0 / 9.0) * GAMMA * P0 * LEFF, rel=1e-9)

@pytest.mark.parametrize("L_km", [10, 50, 125])
def test_e_escala_com_L_eff_e_nao_com_L(L_km):
    """A fase acumulada satura em L_eff (~21,7 km para 0,2 dB/km), nao cresce com L."""
    Lm = L_km * 1e3
    A = gaussian(T, T0, 1e-2)[None, :]
    out = propagate(A, OMEGA, FiberKernelParams(Lm, 0.0, GAMMA, ALPHA, manakov=False))
    phi = _fase_nl(A[0], out[0])
    assert phi.max() == approx(GAMMA * 1e-2 * effective_length(Lm, 0.2e-3),
                                      rel=1e-9)
