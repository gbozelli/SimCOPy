"""Blocos do DSP do receptor coerente que nao precisam dos simbolos transmitidos
(docs/sistema_coerente.md, secao 11; kernels/dsp.py). A parte que usa os simbolos conhecidos
(instante, polarizacao, fase, equalizador) fica na medida (kernels/demodulation.py).

O desvio de frequencia tambem e estimado na medida, com os simbolos: o bloco
FrequencyOffsetCompensator daqui e opcional (rx.dsp.frequency_offset) e fragil para 16QAM.

Os tres blocos trabalham na grade da simulacao, sobre a corrente complexa I + i*Q de cada
polarizacao. No modo com simbolos conhecidos nao ha ADC: a reamostragem para 2 amostras por
simbolo acontece na medida, sem quantizacao (duvida 40). O ADC e o receptor cego ficam em
components/rx/adc.py e vem depois destes blocos (rx.dsp.mode = "blind").
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from ...core.context import SimContext
from ...core.signals import Domain, Waveform
from ...kernels.dsp import dispersion_compensation_response, estimate_frequency_offset
from ...kernels.pulse_shaping import root_raised_cosine_response

@dataclass(frozen=True)
class DispersionCompensator:
    """Desfaz a dispersao acumulada do enlace (Agrawal, secao 8.7.3), com os valores verdadeiros
    das fibras (decisao D-7 do pacote: sem estimativa cega).

    beta2_length: soma de beta2*L das fibras [s^2]. frequency_shift: frequencia do oscilador local
    em relacao ao centro da banda de simulacao [Hz] (desvio do canal mais o desvio do oscilador
    local); com ela a compensacao tira tambem o atraso do canal em relacao ao centro da banda.
    """
    beta2_length: float
    frequency_shift: float = 0.0
    id: str = "rx.cd"
    beta3_length: float = 0.0         # soma de beta3*L das fibras [s^3]
    output_unit: ClassVar[str] = "A"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.ELECTRICAL, self.id)
        h = dispersion_compensation_response(x.grid.freq, self.beta2_length, self.frequency_shift,
                                             self.beta3_length)
        return x.with_samples(np.fft.ifft(np.fft.fft(x.samples, axis=1) * h, axis=1))

@dataclass(frozen=True)
class FrequencyOffsetCompensator:
    """Estima o desvio de frequencia entre o sinal e o oscilador local pela raia da quarta
    potencia (kernels/dsp.py) e o desfaz. O desvio e procurado ate symbol_rate/16 de cada lado."""
    symbol_rate: float
    id: str = "rx.frequency"
    output_unit: ClassVar[str] = "A"

    def offset(self, x: Waveform) -> float:
        """Desvio estimado [Hz]: o sinal gira como exp(+i*2*pi*desvio*t)."""
        return estimate_frequency_offset(x.samples, x.grid.sampling_rate, self.symbol_rate / 16.0)

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.ELECTRICAL, self.id)
        return x.with_samples(x.samples * np.exp(-2j * np.pi * self.offset(x) * x.grid.time))

@dataclass(frozen=True)
class MatchedFilter:
    """Filtro casado ao pulso do transmissor: raiz de cosseno levantado de fator `rolloff`. O
    produto com o filtro do transmissor e o cosseno levantado, que nao tem interferencia entre
    simbolos."""
    symbol_rate: float
    rolloff: float
    id: str = "rx.filter"
    output_unit: ClassVar[str] = "A"

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.ELECTRICAL, self.id)
        h = root_raised_cosine_response(x.grid.freq, self.symbol_rate, self.rolloff)
        return x.with_samples(np.fft.ifft(np.fft.fft(x.samples, axis=1) * h, axis=1))
