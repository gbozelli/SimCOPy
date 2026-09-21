"""Como simular, separado de o que simular (PR-2)."""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class SimulationConfig:
    sampling_rate: float                  # unica no sistema inteiro
    simulation_window: float              # s
    root_seed: int = 0
    ssfm_max_step: float = 5e3            # m
    ssfm_max_phase_change_deg: float = 0.5
    force_power_of_two: bool = True
    broadening_margin: float = 1.2        # [D-6]: constante configurada
    shot_noise: bool = True
    thermal_noise: bool = True

    @classmethod
    def from_symbols(cls, symbol_rate: float, samples_per_symbol: int,
                     n_symbols: int, **kw) -> "SimulationConfig":
        """Atalho no estilo do VPI: SampleRate = SymbolRate*SamplesPerSymbol,
        TimeWindow = NumberOfSymbols/SymbolRate."""
        return cls(sampling_rate=symbol_rate * samples_per_symbol,
                   simulation_window=n_symbols / symbol_rate, **kw)
