"""Conversoes entre a formulacao em numero (Agrawal) e em densidade (notebook do
professor), e de densidade de fotons para potencia (docs/laser.md, secao 3.6).

As duas formulacoes sao a mesma equacao dividida por V: S = P/V, N' = N/V,
Gamma*v_g*sigma_g = G_N*V e N_tr = N_0/V (Agrawal, texto apos a Eq. 3.3.3).
"""
from __future__ import annotations

def number_to_density(x, volume):
    """Numero de particulas (eletrons ou fotons) -> densidade [1/m^3]."""
    return x / volume

def density_to_number(x, volume):
    """Densidade [1/m^3] -> numero de particulas."""
    return x * volume

def gain_coefficient_from_number_form(gain_number_form, volume):
    """G_N do Agrawal [1/s] -> Gamma*v_g*sigma_g [m^3/s] da formulacao em densidade."""
    return gain_number_form * volume

def facet_power(photon_density, volume, group_velocity, mirror_loss, photon_energy):
    """Potencia emitida por faceta [W]: Pe = (1/2)*v_g*alpha_mir*h*f*P, com P = S*V.

    Agrawal, Eq. 3.3.8 (laser Fabry-Perot com facetas iguais).
    """
    return 0.5 * group_velocity * mirror_loss * photon_energy * photon_density * volume
