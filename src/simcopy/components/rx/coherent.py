"""Receptor coerente de uma polarizacao: oscilador local, hibrida de 90 graus e deteccao balanceada
(docs/sistema_coerente.md, secao 3; kernels/coherent.py).

O oscilador local livre (laser proprio, com largura de linha e desvio de frequencia) e o bloco
CWLaser de components/tx/laser.py, com o papel de semente "lo_laser_phase". Aqui fica o oscilador
local travado no laser do transmissor, o homodino ideal da duvida 54.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...config.system import PhotodiodeConfig
from ...core.context import SimContext
from ...core.signals import Domain, SignalError, Waveform
from ...kernels.coherent import hybrid_90
from ...kernels.filters import apply_response, lowpass_response
from ...kernels.photodetector import photocurrent, receiver_noise
from ...kernels.wdm import frequency_shift

@dataclass(frozen=True)
class LockedLocalOscillator:
    """Oscilador local derivado do laser do transmissor: o mesmo campo, com a potencia `power`.

    A entrada e a saida do laser CW do transmissor. O ruido de fase e identico ao do sinal e se
    cancela na mistura: e o homodino ideal (Agrawal, secao 4.5.2, com phi_s = phi_LO), sem
    recuperacao de portadora. `frequency_offset` [Hz, positivo = azul] desloca o campo e deixa um
    desvio de frequencia puro, sem ruido de fase relativo.
    """
    power: float                      # W
    frequency_offset: float = 0.0
    id: str = "rx.lo"
    output_unit: ClassVar[str] = "sqrt(W)"

    def __call__(self, tx_laser: Waveform, ctx: SimContext) -> Waveform:
        tx_laser.require(Domain.OPTICAL, self.id)
        if tx_laser.power <= 0.0:
            raise SignalError(f"{self.id}: o laser do transmissor chegou sem potencia.")
        s = tx_laser.samples * np.sqrt(self.power / tx_laser.power)
        if self.frequency_offset != 0.0:
            s = frequency_shift(s, self.frequency_offset, tx_laser.grid.time)
        return tx_laser.with_samples(s)

@dataclass(frozen=True)
class CoherentFrontEnd:
    """(sinal optico, oscilador local) -> corrente complexa I + i*Q [A], Waveform eletrico (1, N).

    Hibrida de 90 graus e quatro fotodiodos, como no codigo do professor (`DPQAM_receiver`): cada
    fotodiodo tem o seu ruido de disparo e termico (ctx.rng_for(id, "noise0") a "noise3"), as
    correntes sao subtraidas duas a duas e a diferenca passa pelo passa-baixas do receptor, de
    banda photodiode.bandwidth_factor * symbol_rate. Sem ruido, I + i*Q = R*A_s*conj(A_LO).

    n_pol = 1, receptor de uma polarizacao: mistura a componente X do sinal com o oscilador local
    inteiro. A potencia do sinal em Y nao e detectada.
    n_pol = 2, diversidade de polarizacao (Agrawal, secao 8.7.3 e Fig. 8.39): o sinal e separado em
    X e Y, o oscilador local e dividido meio a meio (professor: `X_fraction = 0.5`), e cada
    polarizacao tem a sua hibrida e os seus quatro fotodiodos. A saida tem uma linha por
    polarizacao, na base do receptor: depois da fibra, cada linha mistura as duas polarizacoes
    transmitidas, e quem as separa e o DSP.
    """
    photodiode: PhotodiodeConfig
    symbol_rate: float
    id: str = "rx.frontend"
    n_pol: int = 1
    output_unit: ClassVar[str] = "A"

    def __call__(self, xs, ctx: SimContext) -> Waveform:
        if not isinstance(xs, (list, tuple)) or len(xs) != 2:
            raise SignalError(f"{self.id}: espera duas entradas, o sinal e o oscilador local.")
        signal, lo = xs
        signal.require(Domain.OPTICAL, self.id)
        lo.require(Domain.OPTICAL, self.id)
        c, grid = self.photodiode, signal.grid
        if self.n_pol == 2 and signal.n_pol != 2:
            raise SignalError(f"{self.id}: diversidade de polarizacao pede campo (2, N).")
        h = lowpass_response(grid.freq, c.bandwidth_factor * self.symbol_rate,
                             c.filter_type, c.filter_order)
        lo_field = lo.samples[0] / np.sqrt(self.n_pol)
        out = []
        for p in range(self.n_pol):
            currents = []
            for k, port in enumerate(hybrid_90(signal.samples[p], lo_field)):
                i = photocurrent(port[None, :], c.responsivity, c.dark_current)
                # com uma polarizacao, os papeis de semente sao os da primeira versao
                role = f"noise{k}" if self.n_pol == 1 else f"noise{'xy'[p]}{k}"
                currents.append(i + receiver_noise(i, grid.sampling_rate,
                                                   ctx.rng_for(self.id, role),
                                                   thermal=c.thermal_noise_density,
                                                   shot=c.shot_noise))
            in_phase = apply_response(currents[0] - currents[1], h)
            quadrature = apply_response(currents[2] - currents[3], h)
            out.append(in_phase + 1j * quadrature)
        return Waveform(samples=np.array(out), grid=grid, domain=Domain.ELECTRICAL,
                        geometry=signal.geometry, channel_id=signal.channel_id)
