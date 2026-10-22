"""Resposta de pequeno sinal do laser de equacoes de taxa (docs/laser.md, passo 5).

Agrawal, Eq. 3.3.17-3.3.23, com a linearizacao completa do ganho com compressao
(decisao de 2026-10-15). O livro descarta termos de ordem eps*P, validos so para
eps*P << 1 (Eq. 3.3.15); aqui eles sao mantidos. Com eps = 0 as expressoes sao
identicas as do livro.

Em torno do estacionario (N_b, S_b), com G(N, S) = g_c*(N - N_tr)*(1 - eps*S):

    d(dS)/dt = -Gamma_P*dS + a_SN*dN
    d(dN)/dt = dI/(qV) - a_NS*dS - Gamma_N*dN

    a_SN    = S_b*dG/dN = S_b*g_c*(1 - eps*S_b)
    a_NS    = G_b + S_b*dG/dS
    Gamma_P = R_sp/S_b - S_b*dG/dS          (= R_sp/S + eps*G*S do Agrawal, a menos de eps*S)
    Gamma_N = 1/tau_c + S_b*dG/dN           (= 1/tau_c + G_N*P do Agrawal, a menos de eps*S)

Com dI ~ exp(i*omega*t): dS/dI = (a_SN/qV) / (Omega_R^2 + Gamma_R^2 - omega^2 + 2i*omega*Gamma_R),
Gamma_R = (Gamma_P + Gamma_N)/2 e Omega_R^2 + Gamma_R^2 = Gamma_P*Gamma_N + a_SN*a_NS
(Eq. 3.3.19, 3.3.20). Se Omega_R^2 < 0 a resposta e superamortecida: nao ha oscilacao
de relaxacao, mas a funcao de transferencia continua valendo.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ...core.units import Q_ELECTRON
from .rate_equations import RateEquationParams, gain, steady_state

@dataclass(frozen=True)
class SmallSignal:
    carrier_density: float      # N_b [1/m^3]
    photon_density: float       # S_b [1/m^3]
    a_sn: float                 # [1/s]
    a_ns: float                 # [1/s]
    gamma_p: float              # [1/s]
    gamma_n: float              # [1/s]

    @property
    def gamma_r(self) -> float:
        """Taxa de amortecimento das oscilacoes de relaxacao [1/s] (Eq. 3.3.20)."""
        return 0.5 * (self.gamma_p + self.gamma_n)

    @property
    def omega_r_squared(self) -> float:
        """Omega_R^2 [rad^2/s^2]; negativo quando a resposta e superamortecida."""
        return self.a_sn * self.a_ns - 0.25 * (self.gamma_p - self.gamma_n) ** 2

    @property
    def omega_r(self) -> float:
        """Frequencia angular das oscilacoes de relaxacao [rad/s]; nan se superamortecido."""
        w2 = self.omega_r_squared
        return float(np.sqrt(w2)) if w2 > 0 else float("nan")

def linearize(current, p: RateEquationParams) -> SmallSignal:
    """Coeficientes da linearizacao em torno do estacionario em `current` [A]."""
    n, s = steady_state(current, p)
    g = gain(n, s, p)
    dg_dn = p.gain_coefficient * (1.0 - p.gain_compression * s)
    dg_ds = -p.gain_coefficient * (n - p.transparency_density) * p.gain_compression
    return SmallSignal(
        carrier_density=n, photon_density=s,
        a_sn=s * dg_dn,
        a_ns=g + s * dg_ds,
        gamma_p=p.spontaneous_emission_rate / s - s * dg_ds,
        gamma_n=1.0 / p.carrier_lifetime + s * dg_dn)

def modulation_response(frequency, current, p: RateEquationParams):
    """dS/dI complexo [1/(m^3 A)] na frequencia de modulacao `frequency` [Hz] (Eq. 3.3.19).

    Convencao dI ~ exp(+i*omega*t): a fase negativa e atraso da resposta.
    """
    ss = linearize(current, p)
    w = 2.0 * np.pi * np.asarray(frequency, dtype=float)
    a = ss.gamma_p * ss.gamma_n + ss.a_sn * ss.a_ns          # Omega_R^2 + Gamma_R^2
    return (ss.a_sn / (Q_ELECTRON * p.volume)) / (a - w ** 2 + 2j * w * ss.gamma_r)

def normalized_modulation_response(frequency, current, p: RateEquationParams):
    """H(omega) = dS(omega)/dS(0), com H(0) = 1 (Eq. 3.3.22)."""
    return modulation_response(frequency, current, p) / modulation_response(0.0, current, p)

def bandwidth_3db(current, p: RateEquationParams) -> float:
    """Frequencia [Hz] em que |H| cai a 1/2 (Eq. 3.3.23).

    Com a = Omega_R^2 + Gamma_R^2: omega_3dB^2 = (a - 2*Gamma_R^2) + sqrt((a - 2*Gamma_R^2)^2 + 3*a^2),
    que e a Eq. 3.3.23 reescrita; vale tambem no regime superamortecido.
    """
    ss = linearize(current, p)
    a = ss.gamma_p * ss.gamma_n + ss.a_sn * ss.a_ns
    b = a - 2.0 * ss.gamma_r ** 2
    return float(np.sqrt(b + np.sqrt(b ** 2 + 3.0 * a ** 2)) / (2.0 * np.pi))
