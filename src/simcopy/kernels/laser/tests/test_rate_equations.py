"""Equacoes de taxa: ganho, derivadas e estado estacionario (docs/laser.md, passo 3).

Dois conjuntos de parametros, como decidido em 2026-10-08:
- professor (notebook 2_Laser_sources, formulacao em densidade);
- Agrawal, exercicios 3.7, 3.8 e 3.10 (formulacao em numero, convertida com V arbitrario).
Os valores esperados vem das formas fechadas do Agrawal (Eq. 3.3.4-3.3.10), calculados
a parte e escritos aqui como numeros.
"""
import dataclasses
import pytest
from simcopy.kernels.laser.conversions import facet_power
from simcopy.kernels.laser.rate_equations import (RateEquationParams, gain, derivatives,
                                                  steady_state)
from simcopy.testing import approx

# Professor: L*W*H = 300*3*0,2 um^3 = 1,8e-16 m^3; v_g = c/3,4;
# tau_p = 1/(v_g*(2410 + 1510)) = 2,8931580 ps; g_c = 1*v_g*4e-20 = 3,5269701e-12 m^3/s;
# R_sp = 1*1,5/(2,1e-9*1,8e-16)/1e40 = 3,9682540e-16 1/(m^3 s)
PROFESSOR = RateEquationParams(
    gain_coefficient=3.526970094117647e-12, transparency_density=1.25e24,
    photon_lifetime=2.8931579685554005e-12, carrier_lifetime=2.1e-9, volume=1.8e-16,
    spontaneous_emission_rate=3.968253968253968e-16, linewidth_enhancement=3.0)
N_TH_PROFESSOR = 1.348e24          # N_tr + (alpha_mir + alpha_int)/(Gamma*sigma_g) = 1,25e24 + 9,8e22
I_TH_PROFESSOR = 0.01851200659398857   # q*N_th*V/tau_c = 18,512 mA

def agrawal(volume, r_sp_number=0.0):
    """Agrawal Ex. 3.7: G_N = 6e3 1/s, N_0 = 1e8, tau_c = 2 ns, tau_p = 1,2593408 ps
    (L = 250 um, alpha_int = 40 1/cm, n = 3,3, n_g = 3,4). Em densidade: g_c = G_N*V,
    N_tr = N_0/V, R_sp/V."""
    return RateEquationParams(
        gain_coefficient=6e3 * volume, transparency_density=1e8 / volume,
        photon_lifetime=1.2593408321736882e-12, carrier_lifetime=2e-9, volume=volume,
        spontaneous_emission_rate=r_sp_number / volume)
I_TH_AGRAWAL = 0.018612835881208975   # Ex. 3.8: q*N_th/tau_c, N_th = N_0 + 1/(G_N*tau_p)

def test_gain():
    # G = g_c*(N - N_tr)*(1 - eps*S) = 3,5269701e-12*(1,3e24 - 1,25e24) = 1,7634850e11 1/s
    assert gain(1.3e24, 1e20, PROFESSOR) == approx(176348504705.88168, rel=1e-12)
    # com compressao eps = 1e-22 m^3 e S = 1e20: fator (1 - 1e-2)
    comprimido = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    assert gain(1.3e24, 1e20, comprimido) == approx(174585019658.82288, rel=1e-12)

def test_derivatives():
    # Ponto fora do equilibrio: N = 1,3e24, S = 1e20, I = 30 mA, eps = 0.
    #   dS/dt = (G - 1/tau_p)*S + R_sp = (1,763485e11 - 3,456432e11)*1e20 = -1,6929456e31
    #   dN/dt = I/(qV) - G*S - N/tau_c = 1,0402601e33 - 1,763485e31 - 6,1904762e32 = 4,0356904e32
    #   dphi/dt = (1/2)*alpha_H*(G - 1/tau_p) = 1,5*(-1,6929456e11) = -2,5394185e11 rad/s
    dn, ds, dphi = derivatives(1.3e24, 1e20, 30e-3, PROFESSOR)
    assert dn == approx(4.0356904289192e32, rel=1e-9)
    assert ds == approx(-1.6929456451764774e31, rel=1e-9)
    assert dphi == approx(-253941846776.47162, rel=1e-9)

def test_derivatives_fase_sem_compressao():
    # Agrawal Eq. 3.3.16: a fase usa G_N*(N - N_0), sem o fator (1 - eps*S)
    comprimido = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    assert derivatives(1.3e24, 1e20, 30e-3, comprimido)[2] == approx(-253941846776.47162, rel=1e-9)

@pytest.mark.parametrize("corrente, s_esperado", [
    (20e-3, 1.4927609132902528e20),    # tau_p*(I - I_th)/(qV), Eq. 3.3.7 em densidade
    (30e-3, 1.1524800754727474e21),
])
def test_steady_state_acima_do_limiar(corrente, s_esperado):
    n, s = steady_state(corrente, PROFESSOR)
    assert s == approx(s_esperado, rel=1e-12)
    assert n == approx(N_TH_PROFESSOR, rel=1e-12)   # o ganho trava N no limiar

def test_steady_state_abaixo_do_limiar():
    # I = 10 mA < I_th: N = tau_c*I/(qV) = 7,2817606e23 (abaixo ate da transparencia).
    # S = R_sp*tau_p/(1 - G*tau_p), G = g_c*(N - N_tr) < 0  ->  1,8152203e-28 m^-3
    n, s = steady_state(10e-3, PROFESSOR)
    assert n == approx(7.281760586870891e23, rel=1e-12)
    assert s == approx(1.815220279146239e-28, rel=1e-9)

def test_steady_state_sem_emissao_espontanea_abaixo_do_limiar():
    # R_sp = 0: abaixo do limiar S = 0 exatamente (Agrawal, texto antes da Eq. 3.3.6)
    sem_rsp = dataclasses.replace(PROFESSOR, spontaneous_emission_rate=0.0)
    n, s = steady_state(10e-3, sem_rsp)
    assert s == 0.0
    assert n == approx(7.281760586870891e23, rel=1e-12)

def test_steady_state_agrawal_exercicio_3_8():
    # I = 2*I_th: P = (tau_p/q)*(I - I_th) = 1,4630037e5 fotons; N = N_th = 2,3234437e8
    p = agrawal(volume=1e-16)
    n, s = steady_state(2 * I_TH_AGRAWAL, p)
    assert s * p.volume == approx(146300.37494201775, rel=1e-12)
    assert n * p.volume == approx(232344368.10802943, rel=1e-12)

def test_steady_state_nao_depende_do_volume_escolhido():
    # A formulacao em densidade e a do Agrawal dividida por V: o numero de fotons nao
    # pode depender do V usado na conversao.
    for volume in (1e-16, 1.8e-16, 1.0):
        p = agrawal(volume)
        assert steady_state(2 * I_TH_AGRAWAL, p)[1] * volume == approx(146300.37494201775,
                                                                       rel=1e-12)

def test_steady_state_agrawal_exercicio_3_10_com_emissao_espontanea():
    # R_sp = 2/tau_p = 1,5881324e12 fotons/s. Com eps = 0 o estacionario e a raiz positiva
    # de (G_N*tau_c/tau_p)*P^2 + (1/tau_p - G_N*tau_c*(I/q + R_sp) + G_N*N_0)*P - R_sp = 0
    # -> P = 1,4630351e5 (pouco acima de 1,4630037e5 sem R_sp) e N = 2,3234256e8
    p = agrawal(volume=1e-16, r_sp_number=1588132417296.3528)
    n, s = steady_state(2 * I_TH_AGRAWAL, p)
    assert s * p.volume == approx(146303.5141263201, rel=1e-10)
    assert n * p.volume == approx(232342558.93250957, rel=1e-10)

def test_steady_state_zera_as_derivadas_com_compressao():
    # Com compressao de ganho nao ha forma fechada simples; o estado devolvido tem de
    # anular dN/dt e dS/dt, e a compressao exige mais portadores e da menos fotons.
    comprimido = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    n, s = steady_state(30e-3, comprimido)
    dn, ds, _ = derivatives(n, s, 30e-3, comprimido)
    assert abs(dn) * comprimido.carrier_lifetime / n < 1e-12
    assert abs(ds) * comprimido.photon_lifetime / s < 1e-12
    assert n > N_TH_PROFESSOR
    assert s < 1.1524800754727474e21

def test_curva_potencia_corrente():
    # P_e x I acima do limiar: inclinacao (h*f/2q)*alpha_mir/(alpha_mir + alpha_int)
    # = 0,39994642 * 2410/3920 = 0,24588703 W/A por faceta (Eq. 3.3.10, eta_int = 1)
    def pe(corrente):
        s = steady_state(corrente, PROFESSOR)[1]
        return facet_power(s, 1.8e-16, 88174252.35294119, 2410.0, 1.2815779723541474e-19)
    assert pe(20e-3) == approx(0.0003658782784898758, rel=1e-9)
    assert pe(30e-3) == approx(0.0028247485732891897, rel=1e-9)
    assert (pe(30e-3) - pe(20e-3)) / 10e-3 == approx(0.2458870294799314, rel=1e-9)


# ---- integrador no tempo: Euler explicito com subpassos (docs/laser.md, passo 4)
import numpy as np
from simcopy.kernels.laser.rate_equations import integrate

S_30MA = 1.1524800754727474e21     # estacionario em 30 mA (acima)
S_31MA = 1.2528004738871203e21     # tau_p*(31 mA - I_th)/(qV)

def test_integrate_estacionario_e_ponto_fixo():
    # Partindo do estacionario com corrente constante, N e S nao mudam: no ponto fixo as
    # derivadas sao nulas, e o Euler preserva o ponto fixo para qualquer passo.
    n0, s0 = steady_state(30e-3, PROFESSOR)
    n, s, phi = integrate(np.full(2000, 30e-3), 1e-12, PROFESSOR, n0, s0)
    assert n[-1] == approx(N_TH_PROFESSOR, rel=1e-9)
    assert s[-1] == approx(S_30MA, rel=1e-9)
    # dphi/dt = (alpha_H/2)*(g_c*(N - N_tr) - 1/tau_p): fisicamente ~1e-36 rad/s, mas a
    # subtracao de dois numeros de ~3,46e11 deixa ~1e-3 rad/s de arredondamento -> ~1e-12
    # rad em 2 ns. Limite de 1e-9 rad: muito abaixo de qualquer chirp real.
    assert abs(phi[-1]) < 1e-9

def test_integrate_degrau_converge_para_o_estacionario():
    # Degrau de 30 para 31 mA: depois de 8 ns (~20 tempos de amortecimento 1/Gamma_R) o
    # laser esta no estacionario de 31 mA.
    n0, s0 = steady_state(30e-3, PROFESSOR)
    n, s, _ = integrate(np.full(8000, 31e-3), 1e-12, PROFESSOR, n0, s0)
    assert s[-1] == approx(S_31MA, rel=1e-6)
    assert n[-1] == approx(N_TH_PROFESSOR, rel=1e-6)

def test_integrate_oscilacao_de_relaxacao():
    # Agrawal Eq. 3.3.20-3.3.21 em densidade, no estacionario de 31 mA (eps = 0):
    #   Gamma_N = 1/tau_c + g_c*S = 4,8095e8 + 4,4186e9; Gamma_P = R_sp/S ~ 0
    #   Omega_R = sqrt(g_c*S/tau_p - (Gamma_P - Gamma_N)^2/4) = 3,9003400e10 rad/s (6,2076 GHz)
    #   Gamma_R = (Gamma_P + Gamma_N)/2 = 2,4473901e9 1/s
    # Euler (1a ordem, h = 0,1 ps): erro relativo ~ Omega_R*h/2 = 0,2 % na frequencia e
    # ~ Omega_R^2*h/(2*Gamma_R) = 3 % no amortecimento. Tolerancias: 1 % e 5 %.
    n0, s0 = steady_state(30e-3, PROFESSOR)
    dt = 1e-12
    _, s, _ = integrate(np.full(8000, 31e-3), dt, PROFESSOR, n0, s0)
    x = s - S_31MA
    t = np.arange(x.size) * dt
    i = np.where(np.diff(np.sign(x)) != 0)[0][:9]
    cruzamentos = t[i] - x[i] * dt / (x[i + 1] - x[i])
    omega = 2 * np.pi / (2 * np.mean(np.diff(cruzamentos)))
    assert omega == approx(39003400152.31235, rel=0.01)
    picos = [k for k in range(1, x.size - 1) if x[k] > x[k - 1] and x[k] > x[k + 1]][:2]
    gamma = np.log(x[picos[0]] / x[picos[1]]) / ((picos[1] - picos[0]) * dt)
    assert gamma == approx(2447390140.743383, rel=0.05)

def test_integrate_subpassos_equivalem_a_grade_fina():
    # Rampa de 20 para 30 mA em 500 ps. Grade de 1 ps com 10 subpassos de 0,1 ps (I(t)
    # interpolada linearmente) tem de dar o mesmo que a grade de 0,1 ps sem subpassos,
    # olhando uma amostra a cada 10.
    n0, s0 = steady_state(20e-3, PROFESSOR)
    grossa = np.linspace(20e-3, 30e-3, 501)          # passo de 1 ps
    fina = np.linspace(20e-3, 30e-3, 5001)           # passo de 0,1 ps
    n_g, s_g, phi_g = integrate(grossa, 1e-12, PROFESSOR, n0, s0, max_step=0.1e-12)
    n_f, s_f, phi_f = integrate(fina, 0.1e-12, PROFESSOR, n0, s0, max_step=0.1e-12)
    assert np.max(np.abs(s_g / s_f[::10] - 1)) < 1e-9
    assert np.max(np.abs(n_g / n_f[::10] - 1)) < 1e-9
    assert np.max(np.abs(phi_g - phi_f[::10])) < 1e-9
