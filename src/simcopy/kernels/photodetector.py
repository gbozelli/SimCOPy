"""Fotodetector p-i-n: deteccao quadratica, ruido e banda eletrica (docs/fotodetectores.md).

Agrawal, Sistemas de Comunicacao por Fibra Optica, secoes 4.1 e 4.4:
    I_p = R_d*P_in,  R_d = eta*q/(h*nu)                                   (Eq. 4.1.1, 4.1.3)
    disparo: densidade bilateral q*I, sigma_s^2 = 2*q*(I_p + I_d)*delta_f    (Eq. 4.4.2-4.4.5)
    termico: densidade bilateral 2*k_B*T/R_L, sigma_T^2 = (4*k_B*T/R_L)*F_n*delta_f (Eq. 4.4.7-4.4.9)
    delta_f = integral de 0 a infinito de |H(f)|^2 df (banda equivalente de ruido) (Eq. 4.4.4)
Modelo do codigo do professor (notebook 4_EDFAsAndReceivers, secao 2):
    i_pd(t) = h_pd(t) * [R*(|A_x|^2 + |A_y|^2) + n_shot(t) + n_thermal(t)]
O ruido e somado antes do filtro do receptor, como no professor e no Agrawal (Eq. 4.4.4).
"""
from __future__ import annotations
import numpy as np
from ..core.units import Q_ELECTRON, K_BOLTZMANN
# Os filtros moraram aqui ate 2026-11-12; continuam importaveis deste modulo.
from .filters import lowpass_response, noise_bandwidth, apply_response  # noqa: F401

def photocurrent(field, responsivity, dark_current=0.0):
    """Corrente media R*(|A_x|^2 + |A_y|^2) + I_d [A] de um campo (n_pol, N) (Eq. 4.1.1)."""
    return responsivity * np.sum(np.abs(np.asarray(field)) ** 2, axis=0) + dark_current

def thermal_density(temperature, load_resistance, noise_figure=1.0):
    """Densidade unilateral de corrente do ruido termico, sqrt(4*k_B*T*F_n/R_L) [A/sqrt(Hz)],
    de modo que sigma_T^2 = densidade^2 * delta_f (Eq. 4.4.9)."""
    return np.sqrt(4.0 * K_BOLTZMANN * temperature * noise_figure / load_resistance)

def receiver_noise(current, sampling_rate, rng, thermal=0.0, shot=True):
    """Ruido de disparo + termico [A], amostra a amostra, antes do filtro.

    Na grade de taxa f_s o ruido branco ocupa a banda [-f_s/2, f_s/2]: a variancia por amostra e
    a densidade bilateral vezes f_s. Disparo: q*I(t)*f_s, com a corrente instantanea (como o
    professor); termico: thermal^2/2 * f_s, com `thermal` a densidade unilateral [A/sqrt(Hz)].
    """
    current = np.asarray(current, dtype=float)
    var = 0.5 * thermal ** 2 * sampling_rate * np.ones_like(current)
    if shot:
        var = var + Q_ELECTRON * np.maximum(current, 0.0) * sampling_rate
    return np.sqrt(var) * rng.standard_normal(current.shape)
