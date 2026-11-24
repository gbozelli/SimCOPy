"""Conversor analogico-digital e receptor cego como blocos (docs/sistema_coerente.md, secao 12;
kernels/adc.py, kernels/dsp.py).

Os dois trabalham com o sinal na taxa do DSP (`DigitalSignal`), fora da grade da simulacao: o ADC
entrega 2 amostras por simbolo e o receptor cego devolve 1 amostra por simbolo, os simbolos.
Nenhum dos dois usa os dados transmitidos nem sorteia nada.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import ClassVar
from ...config.system import ADCConfig
from ...core.context import SimContext
from ...core.signals import DigitalSignal, Domain, SignalError, Waveform
from ...kernels.adc import sample_and_quantize
from ...kernels.dsp import blind_receiver

@dataclass(frozen=True)
class ADC:
    """Corrente complexa I + i*Q na grade da simulacao -> `DigitalSignal` a adc.samples_per_symbol
    amostras por simbolo, com a parte real e a imaginaria de cada polarizacao quantizadas em
    adc.resolution_bits bits (0 = sem quantizar).

    A amostragem pega uma amostra a cada (amostras por simbolo da grade)/(amostras por simbolo do
    ADC), a partir de `sampling_offset`. Nao ha recuperacao de relogio: o relogio do ADC e o da
    simulacao, e a fase de amostragem, qualquer que seja, fica com o equalizador fracionario.
    """
    adc: ADCConfig
    symbol_rate: float
    sampling_offset: int = 0
    id: str = "rx.adc"
    output_unit: ClassVar[str] = "A"

    def __call__(self, x: Waveform, ctx: SimContext) -> DigitalSignal:
        x.require(Domain.ELECTRICAL, self.id)
        g = x.geometry
        if g.samples_per_symbol % self.adc.samples_per_symbol:
            raise SignalError(f"{self.id}: {g.samples_per_symbol} amostras por simbolo na grade nao "
                              f"sao multiplo das {self.adc.samples_per_symbol} do ADC.")
        decimation = g.samples_per_symbol // self.adc.samples_per_symbol
        samples = sample_and_quantize(x.samples[:g.n_pol], decimation, self.adc.resolution_bits,
                                      self.adc.full_scale_rms, self.sampling_offset % decimation)
        return DigitalSignal(samples=samples, symbol_rate=self.symbol_rate,
                             geometry=replace(g, samples_per_symbol=self.adc.samples_per_symbol),
                             resolution_bits=self.adc.resolution_bits or None,
                             channel_id=x.channel_id)

@dataclass(frozen=True)
class BlindReceiver:
    """`DigitalSignal` a 2 amostras por simbolo -> `DigitalSignal` com os simbolos (1 amostra por
    simbolo), pela cadeia cega de kernels/dsp.py (`blind_receiver`): espaco de Stokes, equalizador
    adaptativo, desvio de frequencia, busca cega de fase e ajuste dirigido por decisao.

    A saida ainda nao esta casada com os dados: cada linha pode ser qualquer das duas polarizacoes
    transmitidas, com um atraso de alguns simbolos e um giro de multiplo de 2*pi/S. Quem resolve
    isso e a medida (kernels/demodulation.py, `measure_blind`), como o `SynchDSP` do VPI.
    """
    modulation: str
    taps: int = 7
    step_size: float = 1e-3
    cma_iterations: int = 4
    rde_iterations: int = 0
    dd_iterations: int = 4
    test_angles: int = 64
    phase_window: int = 21
    stokes_demux: bool = True
    id: str = "rx.dsp"
    output_unit: ClassVar[str] = "adimensional"

    def frequency_offset(self, x: DigitalSignal) -> float:
        """Desvio de frequencia estimado [Hz] (refaz a cadeia)."""
        return self._run(x)[1] * x.symbol_rate

    def _run(self, x: DigitalSignal):
        if x.samples_per_symbol != 2:
            raise SignalError(f"{self.id}: o equalizador fracionario pede 2 amostras por simbolo; "
                              f"chegaram {x.samples_per_symbol} (rx.adc.samples_per_symbol).")
        return blind_receiver(x.samples, self.modulation, taps=self.taps, step=self.step_size,
                              cma_passes=self.cma_iterations, rde_passes=self.rde_iterations,
                              dd_iterations=self.dd_iterations, test_angles=self.test_angles,
                              window=self.phase_window, stokes=self.stokes_demux)

    def __call__(self, x: DigitalSignal, ctx: SimContext) -> DigitalSignal:
        symbols, _ = self._run(x)
        return DigitalSignal(samples=symbols, symbol_rate=x.symbol_rate,
                             geometry=replace(x.geometry, samples_per_symbol=1),
                             resolution_bits=x.resolution_bits, channel_id=x.channel_id)
