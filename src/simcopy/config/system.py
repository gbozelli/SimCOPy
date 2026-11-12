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

@dataclass(frozen=True)
class CWLaserConfig:
    """Laser quase-monocromatico (modelo B, docs/laser.md, secao 2)."""
    power: float = 1e-3           # W
    linewidth: float = 500e3      # Hz

@dataclass(frozen=True)
class RateEquationLaserConfig:
    """Laser monomodo descrito por equacoes de taxa em densidade (docs/laser.md, 3.6).

    So entradas fisicas; tempo de vida de fotons, limiar, R_sp e eficiencia sao
    derivados em `config/resolve.py` (decisao de 2026-10-08). Sem valores padrao:
    nenhum destes parametros tem valor tipico unico na literatura.
    """
    cavity_length: float              # m
    cavity_width: float               # m
    cavity_height: float              # m
    group_index: float                # n_gr
    mirror_loss: float                # alpha_mir [1/m]
    internal_loss: float              # alpha_int [1/m]
    differential_gain: float          # sigma_g [m^2]
    transparency_density: float       # N_tr [1/m^3]
    confinement: float                # Gamma
    carrier_lifetime: float           # tau_c [s]
    spontaneous_emission_factor: float  # n_sp
    gain_compression: float           # epsilon [m^3]
    linewidth_enhancement: float      # alpha_H (beta_c no Agrawal)

@dataclass(frozen=True)
class MZMConfig:
    vpi: float = 5.0                     # V
    vbias: float = 2.5                   # V
    insertion_loss_db: float = 6.0
    extinction_ratio_db: float = 30.0
    s21_bandwidth: float | None = None   # Hz; None = banda infinita, como no VPI (decisao M5, 2026-10-22)
    s21_order: int = 3                   # Bessel (decisao M5)

@dataclass(frozen=True)
class IQImpairments:
    """Impairments deliberados do arquivo de referencia. Sao o termo dominante da
    curva entre 8 e 15 dBm (piso de 19,12 dB) - ver secao 14.14 e C-25."""
    phase_imbalance_deg: float = 0.0
    skew: float = 0.0                    # s, atraso de I em relacao a Q
    gain_imbalance_db: float = 0.0       # ganho do braco Q em relacao ao I (decisao 2026-10-22)

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

@dataclass(frozen=True)
class ElectricalFilterConfig:
    """Passa-baixas eletrico do acionamento do transmissor IMDD (aula: "nao se usa Nyquist, e sim
    um filtro de Bessel"). Banda de -3 dB = bandwidth_factor * symbol_rate, aplicada na FFT com a
    mesma funcao do filtro do fotodetector (kernels/filters.py). Ordem e banda sem padrao: nenhum
    valor tem base unica na literatura (duvida 50)."""
    filter_order: int
    bandwidth_factor: float
    filter_type: Literal["bessel", "butterworth"] = "bessel"

# ---------------------------------------------------- variantes de canal ------
@dataclass(frozen=True)
class DirectModulation:
    """Laser de modulacao direta: I(t) = I_th*(bias_ratio + modulation_ratio*x(t)), com x(t) o
    sinal eletrico em [0, 1] (docs/laser.md, passo 8). Sem padroes (duvida 5a)."""
    laser: RateEquationLaserConfig
    bias_ratio: float                    # I_b/I_th
    modulation_ratio: float              # I_m/I_th
    scheme: Literal["direct"] = "direct"

@dataclass(frozen=True)
class ExternalModulation:
    """Laser CW + MZM em push-pull na quadratura (docs/moduladores.md). O sinal eletrico x(t) em
    [0, 1] vira a tensao de cada braco v = drive_amplitude*V_pi*(2x - 1); com V_b = V_pi/2,
    drive_amplitude = 0,25 leva o OOK de nulo a pico (notebook 07; duvida 11)."""
    scheme: Literal["external"] = "external"
    laser: CWLaserConfig = field(default_factory=CWLaserConfig)
    mzm: MZMConfig = field(default_factory=MZMConfig)
    drive_amplitude: float = 0.25        # pico por braco / V_pi

@dataclass(frozen=True)
class IMDDTransmitter:
    # Padrao: modulacao externa. A direta exige os parametros do laser de equacoes de taxa,
    # que nao tem padrao (duvida 5b, opcao B).
    frontend: Union[DirectModulation, ExternalModulation] = field(
        default_factory=ExternalModulation)
    driver: ElectricalFilterConfig | None = None   # None = NRZ sem filtro
    rolloff: float = 0.2                 # so a regra de Nyquist usa (duvida 50)

@dataclass(frozen=True)
class OpticalFilterConfig:
    """Filtro optico do demux (kernels/wdm.py, bandpass_response; docs/mux_demux.md). Banda de
    3 dB `bandwidth` [Hz] centrada no canal. "rectangular" e a mascara do professor;
    "supergaussian" exige a ordem (1 = gaussiana). Sem padroes: tipo e banda estao na duvida 30."""
    bandwidth: float
    filter_type: Literal["rectangular", "supergaussian"]
    filter_order: int | None = None

    def __post_init__(self):
        if self.filter_type == "supergaussian" and self.filter_order is None:
            raise ConfigError("optical_filter: filter_type = 'supergaussian' exige filter_order.")

@dataclass(frozen=True)
class IMDDReceiver:
    photodiode: PhotodiodeConfig = field(default_factory=PhotodiodeConfig)
    adc: ADCConfig = field(default_factory=ADCConfig)
    optical_filter: OpticalFilterConfig | None = None   # demux; obrigatorio com mais de um canal

@dataclass(frozen=True)
class CoherentTransmitter:
    laser: CWLaserConfig = field(default_factory=CWLaserConfig)
    mzm: MZMConfig = field(default_factory=MZMConfig)
    iq: IQImpairments = field(default_factory=IQImpairments)
    rolloff: float = 0.18

@dataclass(frozen=True)
class CoherentReceiver:
    lo: CWLaserConfig = field(default_factory=CWLaserConfig)
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
    output_power: float | None = None    # W por canal (decisao 2026-10-22); total = N canais x isso
    noise_figure_db: float = 5.0
    # Saturacao do ganho (decisao 2026-11-12): g = g0/(1 + P/Ps), Agrawal, Exercicio 7.2.
    # None = sem saturacao. So no modo fixed_gain, em que gain_db e o ganho de pequeno sinal.
    saturation_power: float | None = None   # W

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
class MuxConfig:
    """Multiplexador WDM (duvida 29, decisao 2026-11-12: as duas formas). Padrao: soma ideal, sem
    filtro nem perda, como o codigo do professor. Com `optical_filter`, cada entrada passa por um
    filtro centrado no seu canal antes da soma, e `insertion_loss_db` e a perda do mux, como o
    WDM_MUX_N_1 do VPI (gaussiano de ordem 3, 70 GHz, 3 dB)."""
    insertion_loss_db: float = 0.0
    optical_filter: OpticalFilterConfig | None = None

@dataclass(frozen=True)
class SystemConfig:
    channels: tuple[ChannelConfig, ...]
    link: LinkConfig
    center_frequency: float = 193.7e12
    channel_spacing: float = 75e9
    mux: MuxConfig = field(default_factory=MuxConfig)

    @property
    def n_channels(self): return len(self.channels)

    @property
    def center_wavelength(self) -> float:
        """Derivado, nunca armazenado: lambda = c / f. Existe um unico dono
        (center_frequency), entao os dois nao podem divergir."""
        from ..core.units import C_LIGHT
        return C_LIGHT / self.center_frequency

    def channel_offsets(self) -> list[float]:
        """f0 fica no nivel comum: o multiplexador precisa da grade inteira."""
        k = [i - (self.n_channels - 1) / 2 for i in range(self.n_channels)]
        return [x * self.channel_spacing for x in k]

    @property
    def link_is_vectorial(self) -> bool:
        """O enlace nao escolhe SP/DP: ele descobre a partir dos canais (4.4.1)."""
        return any(c.pol is Pol.DP for c in self.channels)
