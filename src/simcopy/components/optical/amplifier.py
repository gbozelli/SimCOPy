"""Amplificador optico (EDFA) como bloco (docs/amplificadores.md, secao 6).

Decisoes de 2026-10-22 aplicadas (docs/amplificadores.md, secao 8): n_sp pela Eq. 7.2.15 invertida
(A1), S_ASE em cada polarizacao (A2), campo (1, N) promovido a (2, N) (A3), output_power por canal
(potencia total = n_channels x output_power) e ganho calculado com a potencia total de entrada,
sinal mais ASE acumulada (A4). h*nu na frequencia central e ASE na banda inteira da grade (A5 e
A6, recomendacoes ainda sem resposta).

Saturacao (decisao de 2026-11-12). Agrawal, Exercicio 7.2: o coeficiente de ganho satura como
g = g0/(1 + P/Ps). Integrando dP/dz = g0*P/(1 + P/Ps) no comprimento do amplificador:
    ln(G) + (G - 1)*P_in/Ps = ln(G0),
com G0 o ganho de pequeno sinal (`gain_db`) e P_in a potencia media de entrada: o EDFA responde a
potencia media, nao a pulsos isolados (Agrawal, secao 7.2.2, texto apos a Fig. 7.4).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from scipy.optimize import brentq
from ...config.system import AmplifierConfig, ConfigError
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, Waveform
from ...core.units import db_to_lin
from ...kernels.amplifier import amplify, ase_psd, nsp_from_noise_figure

@dataclass(frozen=True)
class Amplifier:
    amplifier: AmplifierConfig
    frequency: float                  # frequencia central da banda [Hz] (h*nu da ASE)
    n_channels: int = 1               # canais no campo: alvo total = n_channels x output_power
    id: str = "amplifier"
    output_unit: ClassVar[str] = "sqrt(W)"

    def __post_init__(self):
        a = self.amplifier
        if a.mode == "fixed_gain" and a.gain_db is None:
            raise ConfigError(f"{self.id}: mode = 'fixed_gain' exige gain_db.")
        if a.mode == "target_output_power" and a.output_power is None:
            raise ConfigError(f"{self.id}: mode = 'target_output_power' exige output_power "
                              f"(W por canal) ou output_power_dbm.")
        if a.saturation_power is not None:
            if a.saturation_power <= 0:
                raise ConfigError(f"{self.id}: saturation_power deve ser positiva.")
            if a.mode != "fixed_gain":
                raise ConfigError(f"{self.id}: saturation_power so vale no modo 'fixed_gain', em "
                                  f"que gain_db e o ganho de pequeno sinal. No controlado por "
                                  f"potencia o ganho e o que o controle impoe.")

    def gain(self, x: Waveform) -> float:
        """Ganho de potencia linear para a entrada x (saturado, se houver saturation_power)."""
        a = self.amplifier
        if a.mode != "fixed_gain":
            return self.n_channels * a.output_power / x.power
        g0 = db_to_lin(a.gain_db)
        if a.saturation_power is None or g0 <= 1.0:
            return g0
        ratio = x.power / a.saturation_power
        return brentq(lambda g: np.log(g) + (g - 1.0) * ratio - np.log(g0), 1.0, g0,
                      xtol=1e-12, rtol=1e-12)

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, self.id)
        g = self.gain(x)
        if g < 1.0:
            raise SignalError(f"{self.id}: o alvo de {self.n_channels * self.amplifier.output_power:.3e} W "
                              f"exigiria ganho {g:.3f} < 1 (entrada de {x.power:.3e} W).")
        f_n = db_to_lin(self.amplifier.noise_figure_db)
        psd = ase_psd(nsp_from_noise_figure(f_n, g), g, self.frequency) if g > 1.0 else 0.0
        out = amplify(x.samples, g, psd, x.grid.sampling_rate, ctx.rng_for(self.id, "ase"))
        return x.with_samples(out)
