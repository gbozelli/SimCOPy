"""5. DISPERSAO + NAO LINEARIDADE JUNTAS — soliton fundamental (Agrawal, NLFO, sec. 5.2)

Os testes 01 a 04 isolam um efeito por vez. Este e o unico em que dispersao e
nao linearidade atuam juntas, que e o regime real do enlace.

Com beta2 < 0 e potencia de pico P0 = |beta2| / (gamma * T0^2) (ordem N = 1),
o pulso A(0,T) = sqrt(P0) * sech(T/T0) se propaga SEM mudar de forma: a
dispersao e a SPM se cancelam exatamente.

O que se testa:
  a) a forma |A(z,T)|^2 e preservada por varios periodos de soliton;
  b) o SSFM e de SEGUNDA ORDEM: dividir o passo por 2 divide o erro por ~4.
     Um split-step assimetrico (primeira ordem) dividiria so por ~2. E o teste
     direto da correcao C-8.
"""
import numpy as np
import pytest
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import beta2_from_dispersion

N, T0 = 2048, 10e-12
FS = N / (80 * T0)                        # janela de +-40 T0
T = (np.arange(N) - N // 2) / FS
OMEGA = 2 * np.pi * np.fft.fftfreq(N, 1 / FS)
B2 = beta2_from_dispersion(16e-6, 1550e-9)
GAMMA = 1.3e-3
LD = T0**2 / abs(B2)
P0 = abs(B2) / (GAMMA * T0**2)            # N = 1
Z0 = np.pi / 2 * LD                       # periodo do soliton
A0 = (np.sqrt(P0) / np.cosh(T / T0)).astype(complex)[None, :]

def _erro(z, passo):
    """Erro maximo de intensidade, normalizado, com passo FIXO."""
    out = propagate(A0, OMEGA, FiberKernelParams(z, B2, GAMMA, 0.0, manakov=False),
                    max_step=passo, max_phase_deg=1e9, min_step=passo)
    I0, I1 = np.abs(A0[0])**2, np.abs(out[0])**2
    return np.max(np.abs(I1 - I0)) / I0.max()

def test_a_soliton_preserva_a_forma():
    out = propagate(A0, OMEGA, FiberKernelParams(4 * Z0, B2, GAMMA, 0.0, manakov=False))
    I0, I1 = np.abs(A0[0])**2, np.abs(out[0])**2
    assert np.max(np.abs(I1 - I0)) / I0.max() < 1e-3

def test_b_ssfm_e_de_segunda_ordem():
    z = 2 * Z0
    e1, e2 = _erro(z, LD / 10), _erro(z, LD / 20)
    ordem = np.log2(e1 / e2)
    assert ordem > 1.7, f"ordem de convergencia {ordem:.2f}; esperado ~2 (1 = assimetrico)"
