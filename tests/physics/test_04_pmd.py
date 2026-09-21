"""4. PMD — atraso diferencial de grupo (DGD) estatistico

Referencia: Agrawal, Fiber-Optic Communication Systems, cap. 2 (PMD):

    sigma_T^2 = 2 * (dBeta1 * l_c)^2 * [exp(-z/l_c) + z/l_c - 1]

  z >> l_c :  sigma_T ~= D_p * sqrt(z),   D_p = dBeta1 * sqrt(2 * l_c)
  z << l_c :  sigma_T ~= dBeta1 * z       (fibra de birrefringencia uniforme)

O que se testa:
  a) regime longo: DGD rms = D_p * sqrt(z). E o regime do cenario 400ZR
     (80 a 150 km contra l_c = 50 m) e o que precisa bater;
  b) o DGD NAO depende do passo do SSFM. Foi este teste que pegou um erro: com
     passo de 5 m e l_c de 50 m o DGD saia 3,18x maior;
  c) a distribuicao e Maxwelliana: <tau> / tau_rms = sqrt(8/(3*pi)) = 0,921;
  d) UMA rotacao sorteada cobre a esfera de Poincare inteira, os dois hemisferios.
     E o defeito original do rascunho: com 2 angulos em vez de 3, partindo de
     polarizacao linear em X o estado ficava preso no plano S3 = 0. Precisa de
     teste proprio porque, acumulando muitas secoes, ate a rotacao errada acaba
     espalhando o estado, e os testes a-c nao percebem.

O que NAO se testa, e por que: no regime z ~ l_c o nosso modelo (secoes
discretas com rotacao abrupta) e o do Agrawal (birrefringencia com correlacao
exponencial) diferem por ate sqrt(2). Sao modelos diferentes do mesmo fenomeno,
e so o regime longo e relevante aqui. A figura de validacao mostra a curva toda.

Tolerancias: 200 realizacoes dao erro estatistico de ~3 % no valor rms;
usamos 10 % (cerca de 3 desvios).
"""
import numpy as np
import pytest
from tests._approx import approx
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from ._helpers import jones_dgd

DF, NS = 1.25e9, 64
OMEGA = 2 * np.pi * np.fft.fftfreq(NS, 1 / (NS * DF))
D_P = 0.1e-12 / np.sqrt(1e3)                # 0,1 ps/sqrt(km), valor do arquivo VPI
L_C = 50.0

def _dgd(L, n_real, max_step):
    def run(A, r):
        return propagate(A, OMEGA, FiberKernelParams(L, 0.0, 0.0, 0.0, D_P, L_C),
                         rng=np.random.default_rng(10_000 + r), max_step=max_step)
    return jones_dgd(run, n_real, DF, NS)

@pytest.mark.parametrize("L_km", [2, 10, 40])
def test_a_regime_longo_segue_D_p_raiz_de_z(L_km):
    tau = _dgd(L_km * 1e3, 200, max_step=1e3)
    assert np.sqrt(np.mean(tau**2)) == approx(D_P * np.sqrt(L_km * 1e3), rel=0.10)

@pytest.mark.parametrize("passo", [1000.0, 50.0, 5.0])
def test_b_dgd_nao_depende_do_passo_numerico(passo):
    tau = _dgd(2e3, 200, max_step=passo)
    assert np.sqrt(np.mean(tau**2)) == approx(D_P * np.sqrt(2e3), rel=0.10)

def test_c_distribuicao_maxwelliana():
    tau = _dgd(10e3, 400, max_step=1e3)
    razao = np.mean(tau) / np.sqrt(np.mean(tau**2))
    assert razao == approx(np.sqrt(8 / (3 * np.pi)), abs=0.03)


def test_d_uma_rotacao_cobre_os_dois_hemisferios_da_esfera():
    from simcopy.kernels.birefringence import random_su2
    rng = np.random.default_rng(3)
    S = np.empty((20000, 3))
    for i in range(len(S)):
        U = random_su2(rng)
        ax, ay = U[0, 0], U[1, 0]                     # entrada linear em X
        S[i] = [abs(ax)**2 - abs(ay)**2, 2*np.real(ax*np.conj(ay)), 2*np.imag(ax*np.conj(ay))]
    # Haar-uniforme na esfera => cada parametro de Stokes uniforme em [-1, 1]
    for k in range(3):
        h = np.histogram(S[:, k], bins=8, range=(-1, 1))[0] / (len(S) / 8)
        assert np.abs(h - 1).max() < 0.06, f"S{k+1} nao uniforme: {np.round(h, 3)}"
    # e, explicitamente, os dois hemisferios circulares:
    assert 0.47 < np.mean(S[:, 2] > 0) < 0.53
