"""Configuracao do sistema, com uniao discriminada em `kind`.

Secao 4.6.1: o bloqueio mais barato e mais confiavel e tornar o estado invalido
irrepresentavel. Escrever `vpi` num canal de modulacao direta nao e erro de
validacao, e erro de tipo: o campo nao existe la.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Literal, Union

class ConfigError(ValueError): pass

class Pol(Enum):
    SP = 1
    DP = 2
    @property
    def n_pol(self) -> int: return self.value

# --------------------------------------------------------------- modulacao ---
@dataclass(frozen=True)
class Modulation:
    name: str
    bits_per_symbol: int
    gray: bool = True

OOK   = Modulation("OOK", 1)
PAM4  = Modulation("PAM4", 2)
QPSK  = Modulation("QPSK", 2)
QAM16 = Modulation("16QAM", 4)
QAM32 = Modulation("32QAM", 5)
QAM64 = Modulation("64QAM", 6)
IMDD_FORMATS = {"OOK": OOK, "PAM4": PAM4}
DCS_FORMATS = {"QPSK": QPSK, "16QAM": QAM16, "32QAM": QAM32, "64QAM": QAM64}

# ---------------------------------------------------------------- devices ----
@dataclass(frozen=True)
class LaserConfig:
    power: float = 1e-3           # W
    linewidth: float = 500e3      # Hz

@dataclass(frozen=True)
class MZMConfig:
    vpi: float = 5.0                     # V
    vbias: float = 2.5                   # V
    insertion_loss_db: float = 6.0
    extinction_ratio_db: float = 30.0
    s21_bandwidth: float = 40e9          # Hz
    s21_order: int = 3

@dataclass(frozen=True)
class IQImpairments:
    """Impairments deliberados do arquivo de referencia. Sao o termo dominante da
    curva entre 8 e 15 dBm (piso de 19,12 dB) - ver secao 14.14 e C-25."""
    phase_imbalance_deg: float = 0.0
    skew: float = 0.0                    # s, atraso de I em relacao a Q

@dataclass(frozen=True)
class PhotodiodeConfig:
    responsivity: float = 1.0                    # A/W
    thermal_noise_density: float = 10e-12        # A/sqrt(Hz), parametro direto
    dark_current: float = 0.0
    shot_noise: bool = True
    bandwidth_factor: float = 1.15               # banda = fator * symbol_rate (C-22)
    filter_type: Literal["bessel", "butterworth"] = "bessel"   # fase linear (C-21)
    filter_order: int = 4

@dataclass(frozen=True)
class ADCConfig:
    samples_per_symbol: int = 2
    resolution_bits: int = 12

# ---------------------------------------------------- variantes de canal ------
@dataclass(frozen=True)
class DirectModulation:
    scheme: Literal["direct"] = "direct"
    laser: LaserConfig = field(default_factory=LaserConfig)
    chirp_alpha: float = 0.0

@dataclass(frozen=True)
class ExternalModulation:
    scheme: Literal["external"] = "external"
    laser: LaserConfig = field(default_factory=LaserConfig)
    mzm: MZMConfig = field(default_factory=MZMConfig)

@dataclass(frozen=True)
class IMDDTransmitter:
    frontend: Union[DirectModulation, ExternalModulation] = field(
        default_factory=DirectModulation)
    rolloff: float = 0.2

@dataclass(frozen=True)
class IMDDReceiver:
    photodiode: PhotodiodeConfig = field(default_factory=PhotodiodeConfig)
    adc: ADCConfig = field(default_factory=ADCConfig)

@dataclass(frozen=True)
class CoherentTransmitter:
    laser: LaserConfig = field(default_factory=LaserConfig)
    mzm: MZMConfig = field(default_factory=MZMConfig)
    iq: IQImpairments = field(default_factory=IQImpairments)
    rolloff: float = 0.18

@dataclass(frozen=True)
class CoherentReceiver:
    lo: LaserConfig = field(default_factory=LaserConfig)
    lo_frequency_offset: float = 500e6           # [D-20]: fica aqui, nao no comum
    photodiode: PhotodiodeConfig = field(default_factory=PhotodiodeConfig)
    adc: ADCConfig = field(default_factory=ADCConfig)
    optical_bandwidth: float = 70e9

@dataclass(frozen=True)
class IMDDChannel:
    index: int
    symbol_rate: float
    seed: int
    modulation: Modulation = PAM4
    launch_power: float | None = None
    tx: IMDDTransmitter = field(default_factory=IMDDTransmitter)
    rx: IMDDReceiver = field(default_factory=IMDDReceiver)
    kind: Literal["imdd"] = "imdd"
    # pol nao existe: IMDD e sempre SP ([D-15])
    @property
    def pol(self) -> Pol: return Pol.SP

    def __post_init__(self):
        if self.modulation.name not in IMDD_FORMATS:
            raise ConfigError(
                f"Canal {self.index}: modulacao '{self.modulation.name}' nao se aplica.\n"
                f"  Causa: channel.kind = IMDD.\n"
                f"  Formatos disponiveis: {sorted(IMDD_FORMATS)}.\n"
                f"  Se voce quer {self.modulation.name}, mude kind para 'dcs'.")

@dataclass(frozen=True)
class CoherentChannel:
    index: int
    symbol_rate: float
    seed: int
    modulation: Modulation = QAM16
    pol: Pol = Pol.DP
    launch_power: float | None = None
    tx: CoherentTransmitter = field(default_factory=CoherentTransmitter)
    rx: CoherentReceiver = field(default_factory=CoherentReceiver)
    kind: Literal["dcs"] = "dcs"

    def __post_init__(self):
        if self.modulation.name not in DCS_FORMATS:
            raise ConfigError(
                f"Canal {self.index}: modulacao '{self.modulation.name}' nao se aplica.\n"
                f"  Causa: channel.kind = DCS.\n"
                f"  Formatos disponiveis: {sorted(DCS_FORMATS)}.\n"
                f"  Se voce quer {self.modulation.name}, mude kind para 'imdd'.")

ChannelConfig = Union[IMDDChannel, CoherentChannel]

# ------------------------------------------------------------------ enlace ---
@dataclass(frozen=True)
class BirefringenceConfig:
    """[D-17] e [D-18] fechadas com os valores do arquivo VPI."""
    pmd_coefficient: float = 0.1e-12 / 31.62     # s/sqrt(m) = 0,1 ps/sqrt(km)
    correlation_length: float = 50.0             # m
    enabled: bool = True

@dataclass(frozen=True)
class FiberConfig:
    length: float
    attenuation: float = 0.2e-3          # dB/m
    dispersion: float = 16e-6            # s/m^2
    n2: float = 2.6e-20                  # m^2/W
    effective_area: float = 80e-12       # m^2
    birefringence: BirefringenceConfig = field(default_factory=BirefringenceConfig)

@dataclass(frozen=True)
class AmplifierConfig:
    mode: Literal["fixed_gain", "target_output_power"] = "target_output_power"
    gain_db: float | None = None
    output_power: float | None = None    # W, total do campo agregado
    noise_figure_db: float = 5.0

@dataclass(frozen=True)
class SpanConfig:
    fiber: FiberConfig
    amplifier: AmplifierConfig | None = None
    line_compensation: bool = False

@dataclass(frozen=True)
class LinkConfig:
    spans: tuple[SpanConfig, ...]
    booster: AmplifierConfig | None = None     # antes da fibra, define o lancamento
    @property
    def n_spans(self): return len(self.spans)

@dataclass(frozen=True)
class SystemConfig:
    channels: tuple[ChannelConfig, ...]
    link: LinkConfig
    center_frequency: float = 193.7e12
    channel_spacing: float = 75e9

    @property
    def n_channels(self): return len(self.channels)

    def channel_offsets(self) -> list[float]:
        """f0 fica no nivel comum: o multiplexador precisa da grade inteira."""
        k = [i - (self.n_channels - 1) / 2 for i in range(self.n_channels)]
        return [x * self.channel_spacing for x in k]

    @property
    def link_is_vectorial(self) -> bool:
        """O enlace nao escolhe SP/DP: ele descobre a partir dos canais (4.4.1)."""
        return any(c.pol is Pol.DP for c in self.channels)
