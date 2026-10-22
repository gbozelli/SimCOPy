"""Conversoes numero <-> densidade e fotons -> potencia (docs/laser.md, secao 3.6)."""

from simcopy.kernels.laser.conversions import (
    density_to_number,
    facet_power,
    gain_coefficient_from_number_form,
    number_to_density,
)
from simcopy.testing import approx


def test_number_to_density():
    # N_0 do professor em numero: 1,25e24 m^-3 * 1,8e-16 m^3 = 2,25e8 eletrons.
    assert number_to_density(2.25e8, 1.8e-16) == approx(1.25e24, rel=1e-12)


def test_density_to_number():
    assert density_to_number(1.25e24, 1.8e-16) == approx(2.25e8, rel=1e-12)
    assert number_to_density(density_to_number(3.7e20, 2e-16), 2e-16) == approx(
        3.7e20, rel=1e-12
    )


def test_gain_coefficient_from_number_form():
    # Agrawal Ex. 3.7: G_N = 6e3 1/s; com V = 1e-16 m^3, g_c = 6e-13 m^3/s.
    assert gain_coefficient_from_number_form(6e3, 1e-16) == approx(6e-13, rel=1e-12)


def test_facet_power():
    # Agrawal Ex. 3.8, I = 2*I_th: P = 1,4630037e5 fotons. Pela Eq. 3.3.8,
    # Pe = (1/2)*v_g*alpha_mir*h*f*P = 4,13773 mW por faceta.
    photon_density = 146300.37494201775 / 1e-16
    power = facet_power(
        photon_density,
        1e-16,
        88174252.35294119,
        5005.647198115303,
        1.2815779723541474e-19,
    )
    assert power == approx(0.00413773335687485, rel=1e-9)
