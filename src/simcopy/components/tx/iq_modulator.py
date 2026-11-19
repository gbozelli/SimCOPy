"""Modulador em quadratura do transmissor coerente (docs/sistema_coerente.md, secao 2;
docs/moduladores.md).

Agrawal, secao 10.1.2 e Fig. 10.3: um MZM em push-pull em cada braco de um interferometro externo,
com defasagem de pi/2 entre os bracos. Cada MZM fica polarizado no nulo (V_b = V_pi), em que o
campo e sin(pi*v/V_pi): bipolar, como pedem o BPSK e cada quadratura do QAM (duvida 10; o professor
usa Vpi = Vbias no `DPQAM_transmitter`). O bloco supoe o nulo e nao le `mzm.vbias`: quem garante
vbias = vpi e a regra de validacao `coherent_mzm_at_null` (duvida 10, opcao a, 2026-11-16).

Decisoes de 2026-10-22 aplicadas (docs/moduladores.md, secao 8): extincao finita pelo desbalanco
dos bracos (M4), S21 Bessel na tensao (M5), pre-distorcao em arcsin (M6), skew por fase linear na
FFT (M7), desbalanco de ganho do braco Q. A renormalizacao da saida (M9) fica com o PowerControl
do fim do transmissor, quando o canal tem launch_power, como no IMDD.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...config.system import IQImpairments, MZMConfig
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, Waveform
from ...kernels.filters import apply_response, lowpass_response
from ...kernels.modulator.mzm import arcsin_predistortion, iq_field
from .modulator import arm_amplitudes

@dataclass(frozen=True)
class IQModulator:
    """(portadora optica, sinal eletrico I + i*Q) -> campo optico modulado.

    Entradas, nesta ordem: a saida do laser CW e a saida do acionamento (niveis em +-1 no centro
    dos simbolos). Tensao de cada MZM, para o sinal s de uma quadratura:
    - com `predistortion`: v = (V_pi/pi)*arcsin(drive_depth*s), e o campo fica drive_depth*s. Onde
      |drive_depth*s| passa de 1 (picos entre os simbolos) o sinal e ceifado em +-1;
    - sem: v = (V_pi/2)*drive_depth*s, e o campo e sin[(pi/2)*drive_depth*s].
    Nos dois casos drive_depth = 1 leva os niveis externos de -1 a 1 (o MZM inteiro).

    O campo de saida e (1/2)*[t_I + g_Q*exp(i*theta)*t_Q]*portadora, com theta = pi/2 mais o
    desbalanco de fase (kernels/modulator/mzm.py, `iq_field`).

    Dupla polarizacao: com duas linhas no sinal eletrico, a portadora e dividida meio a meio
    (decisao M8 de 2026-10-22, como o VPI e o professor) e cada metade passa pelo seu modulador
    IQ. A saida tem X na primeira linha e Y na segunda.
    """
    mzm: MZMConfig
    iq: IQImpairments = IQImpairments()
    drive_depth: float = 1.0
    predistortion: bool = True
    id: str = "tx.modulator"
    output_unit: ClassVar[str] = "sqrt(W)"

    def arm_voltage(self, s: np.ndarray) -> np.ndarray:
        """Tensao por braco [V] para o sinal real s de uma quadratura."""
        if self.predistortion:
            limit = 1.0 / self.drive_depth
            return arcsin_predistortion(np.clip(s, -limit, limit), self.mzm.vpi, self.drive_depth)
        return 0.5 * self.mzm.vpi * self.drive_depth * s

    def __call__(self, xs, ctx: SimContext) -> Waveform:
        if not isinstance(xs, (list, tuple)) or len(xs) != 2:
            raise SignalError(f"{self.id}: espera duas entradas, a portadora e o sinal eletrico.")
        carrier, drive = xs
        carrier.require(Domain.OPTICAL, self.id)
        drive.require(Domain.ELECTRICAL, self.id)
        fields = [self._modulate(carrier.samples[0] / np.sqrt(drive.n_pol), s, drive.grid)
                  for s in drive.samples]
        if drive.n_pol == 1:
            fields.append(np.zeros_like(fields[0]))
        return carrier.with_samples(np.array(fields))

    def _modulate(self, carrier: np.ndarray, s: np.ndarray, grid) -> np.ndarray:
        """Um modulador IQ: portadora (N,) e sinal eletrico I + i*Q (N,) -> campo (N,)."""
        m = self.mzm
        v_i, v_q = self.arm_voltage(s.real), self.arm_voltage(s.imag)
        if self.iq.skew != 0.0:       # atraso de I em relacao a Q: fase linear na FFT (M7)
            v_i = apply_response(v_i, np.exp(-1j * grid.omega * self.iq.skew))
        if m.s21_bandwidth is not None:
            h = lowpass_response(grid.freq, m.s21_bandwidth, "bessel", m.s21_order)
            v_i, v_q = apply_response(v_i, h), apply_response(v_q, h)
        a1, a2 = arm_amplitudes(m.insertion_loss_db, m.extinction_ratio_db)
        theta = 0.5 * np.pi + np.deg2rad(self.iq.phase_imbalance_deg)
        q_gain = 10.0 ** (self.iq.gain_imbalance_db / 20.0)
        return iq_field(carrier, v_i, v_q, m.vpi, m.vpi, theta, a1, a2, q_gain)
