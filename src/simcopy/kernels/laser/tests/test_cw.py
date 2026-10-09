"""Laser CW quase-monocromático (docs/laser.md, seção 2.4).

Os testes estatísticos usam seeds fixas (resultado reproduzível) e tolerâncias tiradas
do erro do estimador, escrito em cada teste.
"""
import numpy as np
from simcopy.kernels.laser.cw import cw_field
from simcopy.testing import approx

FS = 1e12          # 1 THz -> dt = 1 ps
DT = 1 / FS

def test_cw_field_potencia_constante():
    # Laser de 1 mW e 500 kHz (valores do VPI): |E|^2 = P em toda amostra. O ruído de
    # fase não muda a potência (Agrawal, Eq. 3.3.33: E = sqrt(P)*exp(i*phi)).
    e = cw_field(4096, DT, 1e-3, 500e3, 0.0, np.random.default_rng(1))
    assert np.max(np.abs(np.abs(e) ** 2 / 1e-3 - 1)) < 1e-12

def test_cw_field_variancia_dos_incrementos_de_fase():
    # Var(dphi) = 2*pi*Delta_nu*dt = 2*pi*500e3*1e-12 = 3,1415927e-6 rad^2 (docs/laser.md,
    # seção 2.2). Com 2e5 incrementos, o desvio relativo do estimador é sqrt(2/2e5) = 0,32 %;
    # tolerância de 1,6 % (5 desvios).
    e = cw_field(200001, DT, 1e-3, 500e3, 0.0, np.random.default_rng(2))
    dphi = np.diff(np.unwrap(np.angle(e)))
    assert np.var(dphi) == approx(3.141592653589793e-06, rel=0.016)

def test_cw_field_autocorrelacao_lorentziana():
    # Agrawal, Eq. 3.3.34 e 3.3.37: |<E*(t)E(t+tau)>|/P = exp(-pi*Delta_nu*|tau|).
    # Delta_nu = 1 GHz: tau = 100 ps -> exp(-0,31416) = 0,73040; 300 ps -> exp(-0,94248) = 0,38966.
    # Média temporal sobre 1 us (tempo de coerência 0,32 ns): erro ~ sqrt(0,32e-9/1e-6) ~ 2 %;
    # tolerância absoluta 0,05.
    e = cw_field(2 ** 20, DT, 1e-3, 1e9, 0.0, np.random.default_rng(3))
    for atraso, esperado in ((100, 0.7304026910486456), (300, 0.3896611373753468)):
        gamma = np.abs(np.mean(np.conj(e[:-atraso]) * e[atraso:])) / 1e-3
        assert abs(gamma - esperado) < 0.05

def test_cw_field_sem_largura_de_linha_e_um_tom():
    # Delta_nu = 0: tom puro em Delta_f. Grade de 1024 amostras a 1,024 THz -> bins de 1 GHz;
    # Delta_f = +10 GHz e -10 GHz caem cada um no seu bin (ordem do NumPy, fftfreq).
    freq = np.fft.fftfreq(1024, 1 / 1.024e12)
    for desvio in (10e9, -10e9):
        e = cw_field(1024, 1 / 1.024e12, 1e-3, 0.0, desvio, np.random.default_rng(4))
        espectro = np.abs(np.fft.fft(e)) ** 2
        k = np.argmin(np.abs(freq - desvio))
        assert espectro[k] / espectro.sum() > 1 - 1e-12

def test_cw_field_fase_inicial_uniforme():
    # Decisão de 2026-10-08: fase inicial uniforme em [0, 2*pi). Com Delta_nu = Delta_f = 0
    # o campo é constante e igual a sqrt(P)*exp(i*phi0). Sobre 4000 seeds, a média de
    # exp(i*phi0) tende a 0 (desvio 1/sqrt(2*4000) = 0,011 por componente; tolerância 0,05)
    # e cada quadrante recebe ~1/4 (desvio 0,007; tolerância 0,03).
    fases = []
    for seed in range(4000):
        e = cw_field(8, DT, 1e-3, 0.0, 0.0, np.random.default_rng(seed))
        assert np.all(e == e[0])
        fases.append(np.angle(e[0]) % (2 * np.pi))
    fases = np.array(fases)
    assert abs(np.mean(np.exp(1j * fases))) < 0.05
    quadrantes = np.histogram(fases, bins=4, range=(0, 2 * np.pi))[0] / len(fases)
    assert np.all(np.abs(quadrantes - 0.25) < 0.03)

def test_cw_field_determinismo():
    # Mesma seed -> mesmo campo; seeds diferentes -> campos diferentes (regra 5 do projeto).
    a = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(7))
    b = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(7))
    c = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(8))
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
