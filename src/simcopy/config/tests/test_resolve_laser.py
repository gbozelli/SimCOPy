"""Grandezas derivadas do laser de equacoes de taxa (docs/laser.md, secao 4).

Parametros do notebook do professor (2_Laser_sources), em 1550 nm. Valores esperados
das formulas do Agrawal (Eq. 3.3.4, 3.3.6, 3.3.10), calculados a parte.
"""
from simcopy.config.system import RateEquationLaserConfig
from simcopy.config.resolve import resolve_rate_equation_laser
from simcopy.testing import approx

PROFESSOR = RateEquationLaserConfig(
    cavity_length=300e-6, cavity_width=3e-6, cavity_height=0.2e-6, group_index=3.4,
    mirror_loss=2410.0, internal_loss=1510.0, differential_gain=4e-20,
    transparency_density=1.25e24, confinement=1.0, carrier_lifetime=2.1e-9,
    spontaneous_emission_factor=1.5, gain_compression=1e-22, linewidth_enhancement=3.0)
F_1550 = 299792458.0 / 1550e-9      # c/lambda = 193,41449 THz

def test_resolve_rate_equation_laser():
    r = resolve_rate_equation_laser(PROFESSOR, F_1550)
    assert r.volume == approx(1.8e-16, rel=1e-12)                       # 300*3*0,2 um^3
    assert r.group_velocity == approx(88174252.35294119, rel=1e-12)      # c/3,4
    # tau_p = 1/(v_g*(alpha_mir + alpha_int)) = 1/(88174252,35*3920) = 2,8931580 ps (Eq. 3.3.4)
    assert r.photon_lifetime == approx(2.8931579685554005e-12, rel=1e-12)
    assert r.gain_coefficient == approx(3.526970094117647e-12, rel=1e-12)  # Gamma*v_g*sigma_g
    # N_th = N_tr + 1/(g_c*tau_p) = 1,25e24 + 3920/4e-20 = 1,348e24 m^-3
    assert r.threshold_density == approx(1.348e24, rel=1e-12)
    # I_th = q*N_th*V/tau_c = 1,6021766e-19*1,348e24*1,8e-16/2,1e-9 = 18,512007 mA (Eq. 3.3.6)
    assert r.threshold_current == approx(0.01851200659398857, rel=1e-12)
    # R_sp = Gamma*n_sp/(tau_c*V)/1e40 (formula do notebook, decisao de 2026-10-08)
    assert r.spontaneous_emission_rate == approx(3.968253968253968e-16, rel=1e-12)
    assert r.photon_energy == approx(1.2815779723541474e-19, rel=1e-12)   # h*c/1550 nm
    assert r.differential_efficiency == approx(0.6147959183673469, rel=1e-12)  # 2410/3920
    # dPe/dI = (h*f/2q)*eta_d = 0,39994642*0,61479592 = 0,24588703 W/A (Eq. 3.3.10)
    assert r.slope_efficiency == approx(0.2458870294799314, rel=1e-12)

def test_resolve_rate_equation_laser_limiar_cresce_com_perda():
    # Mais perda interna -> tau_p menor -> mais portadores para atingir o limiar
    import dataclasses
    mais_perda = dataclasses.replace(PROFESSOR, internal_loss=3020.0)
    r = resolve_rate_equation_laser(mais_perda, F_1550)
    # N_th = 1,25e24 + (2410 + 3020)/4e-20 = 1,38575e24
    assert r.threshold_density == approx(1.38575e24, rel=1e-12)
