"""Modulador externo do transmissor IMDD: laser CW + MZM em push-pull (docs/moduladores.md).

Decisoes de 2026-10-22 aplicadas (docs/moduladores.md, secao 8):
- M4: extincao finita pelo desbalanco dos bracos, a1 + a2 = 10^(-IL/20) e
  (a1 - a2)/(a1 + a2) = 10^(-ER/20);
- M5: S21 Bessel de ordem s21_order na tensao de cada braco; sem s21_bandwidth, banda infinita.
A renormalizacao da saida (M9) fica com o PowerControl do fim do transmissor, quando o canal tem
launch_power.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
from ...config.system import CWLaserConfig, MZMConfig
from ...core.context import SimContext
from ...core.signals import Domain, Waveform
from ...kernels.filters import apply_response, lowpass_response
from ...kernels.modulator.mzm import mzm_field, push_pull
from .laser import CWLaser

def arm_amplitudes(insertion_loss_db: float, extinction_ratio_db: float) -> tuple[float, float]:
    """(a1, a2) com a1 + a2 = 10^(-IL/20) e (a1 - a2)/(a1 + a2) = 10^(-ER/20) (decisao M4).
    IL = 0 e ER infinita dao a1 = a2 = 1/2, o MZM ideal da Eq. 10.1.3."""
    total = 10.0 ** (-insertion_loss_db / 20.0)
    diff = total * 10.0 ** (-extinction_ratio_db / 20.0)
    return 0.5 * (total + diff), 0.5 * (total - diff)

@dataclass(frozen=True)
class ExternalModulator:
    """Sinal eletrico x(t) em [0, 1] -> campo optico (2, N), polarizacao X.

    Tensao de cada braco v = drive_amplitude*V_pi*(2x - 1), em push-pull em torno de V_b
    (`mzm.vbias`). Com V_b = V_pi/2 (quadratura) e drive_amplitude = 0,25 a transmissao vai de 0
    (x = 0) a 1 (x = 1): e a montagem do notebook 07, que usava 0,5*V_pi*(x - 0,5).
    A portadora e o bloco CWLaser (ruido de fase de ctx.rng_for(canal, "tx_laser_phase")).
    """
    laser: CWLaserConfig
    mzm: MZMConfig
    frequency: float                  # frequencia absoluta do canal [Hz]
    drive_amplitude: float = 0.25     # pico por braco / V_pi
    id: str = "tx.modulator"
    channel_offset: float = 0.0       # desvio do canal em relacao ao centro da banda -> Waveform.f0
    output_unit: ClassVar[str] = "sqrt(W)"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.ELECTRICAL, self.id)
        carrier = CWLaser(self.laser, self.frequency, channel_offset=self.channel_offset)(x, ctx)
        m = self.mzm
        v = self.drive_amplitude * m.vpi * (2.0 * x.samples[0].real - 1.0)
        if m.s21_bandwidth is not None:
            v = apply_response(v, lowpass_response(x.grid.freq, m.s21_bandwidth, "bessel",
                                                   m.s21_order))
        v1, v2 = push_pull(v, m.vbias)
        a1, a2 = arm_amplitudes(m.insertion_loss_db, m.extinction_ratio_db)
        return carrier.with_samples(mzm_field(carrier.samples, v1, v2, m.vpi, a1, a2))
