"""Equacoes de taxa de um laser monomodo, em densidade (docs/laser.md, secao 3).

    dS/dt   = (G - 1/tau_p)*S + R_sp
    dN/dt   = I/(q*V) - G*S - N/tau_c
    dphi/dt = (1/2)*alpha_H*[g_c*(N - N_tr) - 1/tau_p]
    G       = g_c*(N - N_tr)*(1 - eps*S),   g_c = Gamma*v_g*sigma_g

E a forma do Agrawal (Eq. 3.3.1, 3.3.2, 3.3.15, 3.3.16) dividida pelo volume V,
como no notebook do professor. S e N sao densidades [1/m^3]; S e o "fluxo de
fotons" da aula. A fase usa o ganho sem compressao, como na Eq. 3.3.16.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from ...core.units import Q_ELECTRON

@dataclass(frozen=True)
class RateEquationParams:
    gain_coefficient: float           # g_c = Gamma*v_g*sigma_g [m^3/s]
    transparency_density: float       # N_tr [1/m^3]
    photon_lifetime: float            # tau_p [s]
    carrier_lifetime: float           # tau_c [s]
    volume: float                     # V [m^3]
    spontaneous_emission_rate: float  # R_sp [1/(m^3 s)]
    gain_compression: float = 0.0     # eps [m^3]
    linewidth_enhancement: float = 0.0  # alpha_H

def gain(n, s, p: RateEquationParams):
    """Taxa liquida de emissao estimulada G [1/s] (Agrawal, Eq. 3.3.15)."""
    return p.gain_coefficient * (n - p.transparency_density) * (1.0 - p.gain_compression * s)

def derivatives(n, s, current, p: RateEquationParams):
    """(dN/dt, dS/dt, dphi/dt) no estado (N, S) com corrente `current` [A]."""
    g = gain(n, s, p)
    dn = current / (Q_ELECTRON * p.volume) - g * s - n / p.carrier_lifetime
    ds = (g - 1.0 / p.photon_lifetime) * s + p.spontaneous_emission_rate
    dphi = 0.5 * p.linewidth_enhancement * (
        p.gain_coefficient * (n - p.transparency_density) - 1.0 / p.photon_lifetime)
    return dn, ds, dphi

def steady_state(current, p: RateEquationParams):
    """(N, S) com dN/dt = dS/dt = 0 para corrente constante [A].

    De dS/dt = 0 vem G*S = S/tau_p - R_sp; levando em dN/dt = 0,
    N(S) = tau_c*(I/(qV) + R_sp - S/tau_p), linear em S. Resta uma equacao em S,
    resolvida por bissecao/Brent no intervalo em que o ganho vai de positivo a
    nao positivo. Com R_sp = 0 e abaixo do limiar, S = 0 e N = tau_c*I/(qV).
    """
    tau_p, r_sp, eps = p.photon_lifetime, p.spontaneous_emission_rate, p.gain_compression
    injection = current / (Q_ELECTRON * p.volume)

    def carriers(s):
        return p.carrier_lifetime * (injection + r_sp - s / tau_p)

    # Limite superior: S em que N(S) = 0 ou em que (1 - eps*S) = 0, o que vier antes.
    # Nos dois o ganho e <= 0, entao dS/dt < 0 ali.
    s_upper = tau_p * (injection + r_sp)
    if eps > 0.0:
        s_upper = min(s_upper, 1.0 / eps)

    if r_sp > 0.0:
        def balance(s):
            return (gain(carriers(s), s, p) - 1.0 / tau_p) * s + r_sp
    else:
        if gain(carriers(0.0), 0.0, p) <= 1.0 / tau_p:
            return carriers(0.0), 0.0
        def balance(s):
            return gain(carriers(s), s, p) - 1.0 / tau_p

    s = brentq(balance, 0.0, s_upper, xtol=1e-300, rtol=1e-15, maxiter=2000)
    return carriers(s), s

def integrate(current, dt, p: RateEquationParams, n0, s0, phi0=0.0, max_step=0.1e-12,
              rng=None):
    """Integra N, S e phi no tempo por Euler explicito (docs/laser.md, secoes 3.4 e 3.5).

    `current` sao as amostras de I(t) [A] na grade de passo `dt`. Cada intervalo da
    grade e dividido em k = ceil(dt/max_step) subpassos de h = dt/k, com I(t)
    interpolada linearmente entre as amostras. Devolve (N, S, phi) nas amostras da
    grade, comecando em (n0, s0, phi0), com o tempo no ultimo eixo.

    Com `rng`, soma as forcas de Langevin (Agrawal, Eq. 3.3.25-3.3.28) por
    Euler-Maruyama: a cada subpasso S ganha sqrt(2*D_SS*h)*xi e phi ganha
    sqrt(2*D_phiphi*h)*xi', xi e xi' normais independentes, com D_SS = R_sp*S e
    D_phiphi = R_sp/(4*S) em densidade. D_NN e os cruzados sao desprezados, como no
    Agrawal. Com poucos fotons a forca gaussiana deixa de valer: S fica limitado a
    >= 0 e D_phiphi usa S >= R_sp*tau_p, isto e, n_sp fotons (proposta, secao 3.5).

    n0, s0 e phi0 podem ser arrays de mesma forma (M,): M realizacoes independentes,
    integradas juntas, para medir estatisticas de conjunto. A saida fica (M, T).
    """
    current = np.asarray(current, dtype=float)
    k = max(1, int(np.ceil(dt / max_step * (1 - 1e-12))))
    h = dt / k
    n, s, phi = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (n0, s0, phi0)))
    n_out = np.empty(n.shape + (current.size,))
    s_out = np.empty(n.shape + (current.size,))
    phi_out = np.empty(n.shape + (current.size,))
    n_out[..., 0], s_out[..., 0], phi_out[..., 0] = n, s, phi
    r_sp = p.spontaneous_emission_rate
    noisy = rng is not None and r_sp > 0.0
    s_floor = r_sp * p.photon_lifetime
    for j in range(current.size - 1):
        i0, i1 = current[j], current[j + 1]
        if noisy:
            xi = rng.standard_normal((k, 2) + n.shape)
        for m in range(k):
            i_m = i0 + (i1 - i0) * m / k
            dn, ds, dphi = derivatives(n, s, i_m, p)
            n_new, s_new, phi_new = n + dn * h, s + ds * h, phi + dphi * h
            if noisy:
                s_new = np.maximum(s_new + np.sqrt(2.0 * r_sp * s * h) * xi[m, 0], 0.0)
                phi_new = phi_new + np.sqrt(r_sp * h / (2.0 * np.maximum(s, s_floor))) * xi[m, 1]
            n, s, phi = n_new, s_new, phi_new
        n_out[..., j + 1], s_out[..., j + 1], phi_out[..., j + 1] = n, s, phi
    return n_out, s_out, phi_out
