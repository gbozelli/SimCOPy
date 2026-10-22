"""Laser CW quase-monocromatico (docs/laser.md, secao 2.4)."""

import numpy as np

from simcopy.kernels.laser.cw import cw_field
from simcopy.testing import approx

FS = 1e12
DT = 1 / FS


def test_cw_field_potencia_constante():
    # Laser de 1 mW e 500 kHz: |E|^2 = P em toda amostra (Agrawal, Eq. 3.3.33).
    field = cw_field(4096, DT, 1e-3, 500e3, 0.0, np.random.default_rng(1))
    assert np.max(np.abs(np.abs(field) ** 2 / 1e-3 - 1)) < 1e-12


def test_cw_field_variancia_dos_incrementos_de_fase():
    # Var(dphi) = 2*pi*Delta_nu*dt = pi*1e-6 rad^2. Com 2e5 incrementos,
    # o desvio relativo do estimador e sqrt(2/2e5) = 0,32 %; usamos 5 desvios.
    field = cw_field(200001, DT, 1e-3, 500e3, 0.0, np.random.default_rng(2))
    dphi = np.diff(np.unwrap(np.angle(field)))
    assert np.var(dphi) == approx(np.pi * 1e-6, rel=0.016)


def test_cw_field_autocorrelacao_lorentziana():
    # |<E*(t)E(t+tau)>|/P = exp(-pi*Delta_nu*|tau|), Agrawal 3.3.34/3.3.37.
    field = cw_field(2**20, DT, 1e-3, 1e9, 0.0, np.random.default_rng(3))
    for delay, expected in ((100, 0.7304026910486456), (300, 0.3896611373753468)):
        gamma = np.abs(np.mean(np.conj(field[:-delay]) * field[delay:])) / 1e-3
        assert abs(gamma - expected) < 0.05


def test_cw_field_sem_largura_de_linha_e_um_tom():
    # Grade de 1024 amostras a 1,024 THz: bins de 1 GHz, sem leakage espectral.
    frequency = np.fft.fftfreq(1024, 1 / 1.024e12)
    for offset in (10e9, -10e9):
        field = cw_field(1024, 1 / 1.024e12, 1e-3, 0.0, offset, np.random.default_rng(4))
        spectrum = np.abs(np.fft.fft(field)) ** 2
        expected_bin = np.argmin(np.abs(frequency - offset))
        assert spectrum[expected_bin] / spectrum.sum() > 1 - 1e-12


def test_cw_field_fase_inicial_uniforme():
    # Com Delta_nu = Delta_f = 0, E/sqrt(P) = exp(i*phi0). Para 4000 sementes,
    # a media circular tende a zero e cada quadrante recebe aproximadamente 1/4.
    phases = []
    for seed in range(4000):
        field = cw_field(8, DT, 1e-3, 0.0, 0.0, np.random.default_rng(seed))
        assert np.all(field == field[0])
        phases.append(np.angle(field[0]) % (2 * np.pi))
    phases = np.array(phases)
    assert abs(np.mean(np.exp(1j * phases))) < 0.05
    quadrants = np.histogram(phases, bins=4, range=(0, 2 * np.pi))[0] / len(phases)
    assert np.all(np.abs(quadrants - 0.25) < 0.03)


def test_cw_field_determinismo():
    # A mesma seed reproduz o campo; seeds diferentes produzem realizacoes diferentes.
    a = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(7))
    b = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(7))
    c = cw_field(1000, DT, 1e-3, 500e3, 0.0, np.random.default_rng(8))
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
