# Documentação da Simulação VPItransmissionMaker — WDM 400G ZR (4 Canais, Foco Ch2)


## Resumo Executivo da Topologia

O esquemático simula um enlace óptico WDM de curta distância (Optical Interconnect) com as seguintes características principais:

- **4 canais WDM**, cada um com **16-QAM** (BitsPerSymbol = 4), polarização dupla (Dual-Polarization), taxa de símbolo de **59,84 GBd** → bit rate de **~239,4 Gb/s por canal** (2 pol × 4 bits/símbolo × 59,84 GBd).
- Espaçamento entre canais: **75 GHz**.
- Frequência central (canal de referência, Canal 3 na numeração interna / "TxCarrierFrequency"): **193,7 THz**.
- Potência de lançamento: **13 dBm** (`PowerdBm`).
- Enlace de fibra: **95 km** (`Length = 95e3` m), atenuação 0,2 dB/km, dispersão 16 ps/nm/km.
- Amplificação: um amplificador óptico controlado por potência, compensando exatamente a perda da fibra (Gain = Length × 0,2e-3 dB), NF = 5 dB.
- Recepção coerente completa (Rx_DualPol) com oscilador local, ADC de 12 bits e cadeia de DSP de 6 estágios (CD compensation → equalização MIMO → demux de polarização (Stokes) → recuperação de clock → correção de offset de frequência → recuperação de fase - Blind Phase Search).
- Medição de desempenho via **BER_4D** (contagem de erros/constelação) para o Canal 2 (Ch_X2/Ch_Y2).

---

# Documentação da Simulação VPItransmissionMaker

## WDM 400G ZR — 4 Canais 16-QAM (Foco no Canal 2)

**Arquivo analisado:** `WDM_of_400G_ZR_4_Channels_Ch2-_versao_editavel.vtmu`

### Metadados do Projeto

| Campo | Valor |
|---|---|
| Descrição | 400G ZR 16-QAM Optical Interconnect transmission system. |
| Palavras-chave | 16QAM, Short reach, WDM, 400G, Optical Interconnect |
| Produto | VPI_TMM_OPTICALSYSTEMS.pro |
| Versão do produto | 11.5 |
| Criado em | 2018-08-02T08:10:08.3880282Z |

> Este é um esquema do VPIphotonics/VPItransmissionMaker que simula um sistema WDM óptico de 4 canais em formato 16-QAM com polarização dupla (400G ZR), com foco na recepção coerente do **Canal 2**. Abaixo estão listados **todos os parâmetros** configurados em cada componente do esquemático, extraídos diretamente do arquivo `.vtmu`.

---
## 1. Parâmetros Globais da Simulação (INTERFACE)

Esses parâmetros são definidos no nível do projeto e referenciados (por nome) pelos demais componentes.

| Parâmetro | Valor | Unidade |
|---|---|---|
| PowerdBm [Physical] | 13 |  |
| Length [Physical] | 95e3 |  |
| TxCarrierFrequency [Physical] | 193.7e12 | Hz |
| SymbolRate [System] | 59.84e9 |  |
| NumberOfSymbols [System] | 2^18 |  |
| BitsPerSymbol [System] | 4 |  |
| BitRateDefault [System] | SymbolRate*BitsPerSymbol | bit/s |
| S21_Bandwidth [MZM] | 40e9 | Hz |
| TimeWindow [Simulation Control] | NumberOfSymbols/SymbolRate | s |
| SamplesPerSymbol [Simulation Control] | 16 |  |
| SampleRateDefault [Simulation Control] | SymbolRate*SamplesPerSymbol | Hz |
| GreatestPrimeFactorLimit [Simulation Control] | 2 |  |
| BoundaryConditions [Simulation Control] | Periodic |  |
| InBandNoiseBins [Simulation Control] | OFF |  |
| LogicalInformation [Simulation Control] | ON |  |
| SampleModeBandwidth [Simulation Control] | SymbolRate*SamplesPerSymbol | Hz |
| SampleModeCenterFrequency [Simulation Control] | 193.7e12 | Hz |
| TrackingMode [Design Rules] | None |  |
| TrackingStepSize [Design Rules] | 1e4 | m |
| OpticalBandwidth [Design Rules] | 4*BitRateDefault | Hz |
| ElectricalBandwidth [Design Rules] | 0.75*BitRateDefault | Hz |


---
## 2. Transmissores Ópticos — Tx_mQAM_DP (4 canais WDM)

O esquema contém **4 instâncias idênticas** do módulo `Tx_mQAM_DP` (transmissor coerente 16-QAM com polarização dupla), uma para cada canal WDM. Elas compartilham a mesma configuração-base; as únicas diferenças entre elas estão na tabela de "Diferenças entre canais" logo abaixo.

### 2.1 Diferenças entre os 4 canais Tx

| Parâmetro | Valores por canal (Tx1..Tx4) | Unidade |
|---|---|---|
| ChannelLabelX | Ch1: Ch_X3 \| Ch2: Ch_X4 \| Ch3: Ch_X2 \| Ch4: Ch_X1 |  |
| ChannelLabelY | Ch1: Ch_Y3 \| Ch2: Ch_Y4 \| Ch3: Ch_Y2 \| Ch4: Ch_Y1 |  |
| RandomNumberSeedX | Ch1: 2 \| Ch2: 0 \| Ch3: 4 \| Ch4: 6 |  |
| RandomNumberSeedY | Ch1: 3 \| Ch2: 1 \| Ch3: 5 \| Ch4: 7 |  |
| EmissionFrequency | Ch1: TxCarrierFrequency-75e9 \| Ch2: TxCarrierFrequency-150e9 \| Ch3: TxCarrierFrequency \| Ch4: TxCarrierFrequency+75e9 | Hz |

### 2.2 Parâmetros comuns aos 4 transmissores (Tx_mQAM_DP)

| Parâmetro | Valor | Unidade |
|---|---|---|
| BaudRate | SymbolRate | Bd |
| ModulationFormat | mQAM |  |
| GrayMapping | No |  |
| FormatDetails |  |  |
| ModulationFormatFilename |  |  |
| BitsPerSymbol | 4 |  |
| PreSpaces | 0 |  |
| PostSpaces | 0 |  |
| PRBS_Type | PRBS |  |
| MarkProbability | 0.5 |  |
| PRBS_Order | 7 |  |
| DiffPreCoding | No |  |
| MarkNumber | 7 |  |
| CodeWord |  |  |
| InputFilenameX |  |  |
| InputFilenameY |  |  |
| OutputFilenameX |  |  |
| OutputFilenameY |  |  |
| ControlFlagReset | Continue |  |
| ControlFlagWrite | Overwrite |  |
| PreDistortion | Yes |  |
| PulsedOutput | No |  |
| Gain | 0 | dB |
| Renormalize | On |  |
| SkewX | 0.75e-12 #adds delay of I w.r.t. Q driver signal for X branch | s |
| SkewY | 0.75e-12 #adds delay of I w.r.t. Q driver signal for Y branch | s |
| FilterActive | Advanced |  |
| LPF_CutoffFrequency | BaudRate | Hz |
| LPF_FilterOrder | 4 |  |
| RiseTime | 1.0/8.0/BaudRate | s |
| FilterType | LowPass |  |
| A |  |  |
| B | 1.0 |  |
| TransferFunction | RaisedCosine |  |
| Bandwidth | 0.75*BaudRate/2 | Hz |
| PassBandRipple | 5e-1 | dB |
| FilterOrderDefinition | FilterOrder |  |
| FilterOrder | 4 |  |
| GaussianOrder | 1 |  |
| StopBandwidth | 2*BaudRate | Hz |
| StopbandAttenuation | 40 | dB |
| FilterDataFilename | data.dat |  |
| DataFormat | Frequency Magnitude |  |
| LogScale | Automatic |  |
| FrequencyShift | 0 | Hz |
| FrequencyScaling | 1.0 |  |
| PowerScaling | 1.0 |  |
| MinimumLoss | 0 | dB |
| FrequencyAxisScaling | Logarithmic |  |
| MinimumPhase | Off |  |
| RollOff | 0.18 | a.u |
| SymbolRate | BaudRate | Bd |
| NyquistResponse | squareRootRaisedCosine |  |
| InterpolationOrder | Cubic |  |
| ExtrapolationType | LastPoint |  |
| SetTransmission | -200 | dB |
| SetPhase | 0 | deg |
| ConserveMemory | ON |  |
| SampleRate | SampleRateDefault | Hz |
| AveragePower | 1e-3 | W |
| Linewidth | 500e3 | Hz |
| LaserRandomNumberSeed | 0 |  |
| QuantizedEmissionFreq | Yes |  |
| SimulationBandDefinition | EmissionFrequency |  |
| SimulationBandCenterFreq | 193.1e12 | Hz |
| S21_Bandwidth | S21_Bandwidth | Hz |
| S21_FilterOrder | 3 |  |
| Vpi | 5 | V |
| InsertionLoss | 6 | dB |
| ExtinctionRatio | 30 | dB |
| PhaseShiftX | -90+5 #Pi/2, add IQ phase imbalance here for X branch | deg |
| PhaseShiftY | -90+5 #Pi/2, add IQ phase imbalance here for X branch | deg |
| SplitterExtinctionRatioX | 30 | dB |
| SplitterExtinctionRatioY | 30 | dB |
| CombinerExtinctionRatioX | 30 | dB |
| CombinerExtinctionRatioY | 30 | dB |
| PowerControl | Off |  |
| OutputPower_dBm | 0 |  |
| OutputDataType | Blocks |  |
| Active | On |  |
| randX | !"expr {int(rand()*10000000)}" |  |
| randY | !"expr {int(rand()*10000000)}" |  |


---
## 3. Multiplexador WDM — WDM_MUX_N_1

Combina os 4 canais transmitidos em uma única fibra (grade de canais).

| Parâmetro | Valor | Unidade |
|---|---|---|
| InsertionLoss | 3 | dB |
| FilterType | BandPass |  |
| TransferFunction | Gaussian |  |
| FilterOrderDefinition | FilterOrder |  |
| FilterOrder | 3 |  |
| GaussianOrder | 3 |  |
| Bandwidth | 70e9 | Hz |
| PassBandRipple | 500e-3 | dB |
| StopBandwidth | 8*SymbolRate | Hz |
| StopbandAttenuation | 40 | dB |
| DataFilename | filter.dat |  |
| DataFormat | Frequency Magnitude Phase |  |
| LogScale | On |  |
| PowerScaling | 1.0 |  |
| FrequencyScaling | 1.0 |  |
| ScalingReferenceFrequency | 193.1e12 | Hz |
| TruncateFrequencyMin | 50e12 | Hz |
| TruncateFrequencyMax | 350e12 | Hz |
| ReferenceFrequency | 193.1e12 | Hz |
| DelayAtReferenceFrequency | 0 | s |
| MinimumLoss | 0 | dB |
| MinimumPhase | Off |  |
| FrequencyAxisScaling | Linear |  |
| ChannelNumber | 4 |  |
| ChannelSpacing | 75e9 |  |
| CenterFrequencyOfBand | TxCarrierFrequency-75e9/2 |  |
| NoiseDynamic | 3.0 | dB |
| NoiseResolution | SymbolRate/4.0 | Hz |
| NoiseThreshold | -100.0 | dB |
| ActiveFilterBandwidth | 2*SampleRateDefault | Hz |
| ConserveMemory | On |  |
| FiltersEnabled | On |  |
| LimitToActiveFilterBandwidth | SampledBands |  |
| DigitalFilterOrderDefinition | Parameter |  |
| DigitalFilterOrder | 64 |  |
| DigitalFilterAccuracyGoal | 1.e-4 |  |

---
## 4. Fibra Óptica de Transmissão — UniversalFiberFwd

Modelo de fibra óptica (não-linear, dispersão, PMD, Raman) que representa o enlace de transmissão.

| Parâmetro | Valor | Unidade |
|---|---|---|
| NumberOfFiberSpans | 1 |  |
| Length | Length | m |
| GroupRefractiveIndex | 1.47 |  |
| AttenuationDescription | AttenuationParameter |  |
| Attenuation | 0.2e-3 | dB/m |
| AttFilename |  |  |
| ReferenceFrequency | TxCarrierFrequency | Hz |
| DispersionDescription | DispersionParameters |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0.08e3 | s/m^3 |
| DispersionFilename |  |  |
| PMDCoefficient | 0.1e-12/31.62 | s/sqrt(m) |
| CorrelationLength | 50.0 | m |
| RamanScattering | No |  |
| RamanFilename | $BNROOT/data/UFM/RamanGain.dat |  |
| SpontaneousRamanScattering | No |  |
| RamanAdjustmentFactors | 0.5 |  |
| Temperature | 300.0 | K |
| NonlinearDescription | NonlinearIndexParameter |  |
| NonLinearIndex | 2.6e-20 | m^2/W |
| NonLinearFilename |  |  |
| CoreAreaDescription | CoreAreaParameter |  |
| CoreArea | 80.0e-12 | m^2 |
| CoreAreaFilename |  |  |
| OverlapIntegralFilename |  |  |
| RamanFraction | 0.17 |  |
| SPM_EC | Yes |  |
| XPM_EC | Yes |  |
| XPM_MC | Yes |  |
| IntraBandRaman | Yes |  |
| InterBandRaman | Yes |  |
| FWM | No |  |
| HigherOrderNLEffects | No |  |
| NonlinearAdjustmentFactors | (8./9.) 1.0 (2./3.) 1.0 (1./2.) 1.0 1.0 |  |
| EventLossDescription | EventLossParameter |  |
| EventLoss | 0.0 | dB |
| EventLossFilename |  |  |
| FreqResolutionNB | 100e9 | Hz |
| FreqResolutionSFB | 100e9 | Hz |
| GridReferenceFrequency | 193.1e12 | Hz |
| SpectralDiscretizerDescription | FixedFrequencies |  |
| FWMThreshold | -30.0 | (dBm) |
| FWMPowerSeed | 1e-10 | W |
| NoiseCenterFrequency | 193.1e12 | Hz |
| NoiseBandwidth | 10e12 | Hz |
| AccuracyGoal | 1e-3 | dB |
| MaximumIterations | 1000 |  |
| IterationAccuracyFactor | 10.0 |  |
| InitialStepSize | 1e3 | m |
| MinimumStepSize | 1 | m |
| StabilizationFactor | 0.5 |  |
| StepSelectionMethod | NonlinearPhaseChange |  |
| SplitStepType | Symmetric |  |
| MaxStepWidth | 5.0e3 | m |
| MinStepWidth | 0.0 | m |
| InitialStepWidth | 1e3 | m |
| TargetLocalError | 2.0e-3 |  |
| MaxPhaseChange | 0.5 | deg |
| MeanStepWidth | 1.0e3 | m |
| WidthDeviation | 10 | m |
| PolarizationAnalysis | VectorPMD |  |
| FieldAnalysis | Yes |  |
| SkipFieldAnalysis | 0 |  |
| NoiseBinInteractions | Yes |  |
| DistortionInteractions | Yes |  |
| NoisePolarization | Unpolarized |  |
| DistortionPolarization | Unpolarized |  |
| WidthRandomNumberSeed | 78 |  |
| ScatteringRandomNumberSeed | 79 |  |
| FreezeWidth | 0 |  |
| FreezeScattering | 0 |  |
| CorrelationScatteringFactor | 0.0 |  |
| BirefringenceProfile | None |  |
| BirefringenceFileDataType | RotationMatrices |  |
| BirefringenceInputFilename |  |  |
| BirefringenceOutputFilename |  |  |
| LikelihoodRatiosFilename |  |  |
| ConserveMemory | No |  |
| NoiseRandomNumberSeed | 0 |  |
| Active | On |  |
| VisualizationMode | None |  |
| Title |  |  |
| SaveToFile | Off |  |
| LogFilename |  |  |
| PlotFrequencyRange | (100e12, 1000e12) | Hz |
| PowerPlotPoints | 101 |  |
| SeparateWindows | Off |  |
| Signals | Off |  |
| NoiseBins | Off |  |
| Distortions | Off |  |
| CrosstalkDistortions | On |  |
| RayleighDistortions | On |  |
| BrillouinDistortions | On |  |
| FWMDistortions | On |  |
| GenericDistortions | On |  |
| StokesVector | Off |  |

---
## 5. Amplificador Óptico — AmpSysOpt

Amplificador óptico (modelo de sistema) usado para compensar a atenuação da fibra.

| Parâmetro | Valor | Unidade |
|---|---|---|
| AmplifierType | PowerControlled |  |
| GainShapeDescription | GainShapeParameters |  |
| Gain | !" expr {Length}*0.2e-3" | dB |
| GainTilt | 0 | dB/Hz |
| GainTiltRefFreq | 193.1e12 | Hz |
| GainShapeFilename |  |  |
| OutputPowerMax | 1e6 | W |
| OutputPower | 10^((PowerdBm-30)/10) | W |
| GainMax | 100.0 | dB |
| SaturationPower | 10e-3 | W |
| SaturationPowerType | OutputPower |  |
| NoiseDescription | NoiseFigureParameters |  |
| NoiseFigure | 5 | dB |
| NoisePower | 0.0 | W |
| NoiseTilt | 0.0 | dB/Hz |
| NoiseTiltRefFreq | 193.1e12 | Hz |
| NoiseFigureFilename |  |  |
| NoisePowerFilename |  |  |
| NoiseBandwidth | 0 | Hz |
| NoiseCenterFrequency | TxCarrierFrequency | Hz |
| PolarizFilter | None |  |
| FrequencyResolution | 100.0e9 | Hz |
| GridReferenceFrequency | 193.1e12 | Hz |
| NoiseDynamic | 3.0 | dB |
| NoiseResolution | BitRateDefault/4 | Hz |
| NoiseThreshold | -100 | dB |
| DigitalFilterOrder | 64 |  |
| SampledSignals | On |  |
| ParameterizedSignals | On |  |
| Distortions | On |  |
| NoiseBins | On |  |
| RandomNumberSeed | 0 |  |
| Active | On |  |

---
## 6. Filtro Óptico de Seleção de Canal — FilterOpt

Filtro óptico passa-banda usado no receptor para selecionar/isolar o canal de interesse (Canal 2) antes da detecção coerente.

| Parâmetro | Valor | Unidade |
|---|---|---|
| FilterType | BandPass |  |
| DataFilename | filter.dat |  |
| NumberOfStepAndRepeat | 1 |  |
| ReA |  |  |
| ImA |  |  |
| ReB | 1.0 |  |
| ImB | 0.0 |  |
| TransferFunction | Gaussian |  |
| CenterFrequency | TxCarrierFrequency | Hz |
| Bandwidth | 70e9 | Hz |
| PassBandRipple | 0.5 | dB |
| FilterOrderDefinition | FilterOrder |  |
| FilterOrder | 3 |  |
| GaussianOrder | 3 |  |
| StopBandwidth | 8*BitRateDefault | Hz |
| StopbandAttenuation | 80 | dB |
| FreeSpectralRange | 100e9 | Hz |
| PowerScaling | 1.0 |  |
| FrequencyShift | 0.0 | Hz |
| FrequencyScaling | 1.0 |  |
| ScalingReferenceFrequency | 193.1e12 | Hz |
| InterpolationOrder | Cubic |  |
| ExtrapolationType | LastPoint |  |
| SetTransmission | -200 | dB |
| SetPhase | 0.0 | deg |
| TruncateFrequencyMin | 50e12 | Hz |
| TruncateFrequencyMax | 350e12 | Hz |
| DataValidThreshold | -200 | dB |
| DataFormat |  |  |
| LogScale | On |  |
| ReferenceFrequency | 193.1e12 | Hz |
| DelayAtReferenceFrequency | 0.0 | s |
| MinimumLoss | 0.0 | dB |
| DigitalFilterOrderDefinition | Automatic |  |
| DigitalFilterAccuracyGoal | 1.0e-4 |  |
| DigitalFilterHighestOrder | 2048 |  |
| DigitalFilterOrder | 64 |  |
| ImpulseResponseLength | 8*DigitalFilterOrder |  |
| ImpulseResponseAutoShift | On |  |
| ImpulseResponseShift | 0 |  |
| WindowFunction | Tukey |  |
| GaussianParameter | 0.0 : 0.2 : 4.0 |  |
| KaiserParameter | 0.0 : 0.05 : 1.0 |  |
| TukeyParameter | 0.0, 0.01, 0.05, 0.10, 0.15 |  |
| OptimizationFrequencies |  | Hz |
| OptimizationMagnitudeWeights | 1.0 | 1/dB |
| OptimizationDelayWeights | 1e12 | 1/s |
| TransferFunctionRegularization | None |  |
| StitchingWidth | 0.1 |  |
| IgnoreDataColumn | None |  |
| MinimumPhase | Off |  |
| ComplementaryPower | Off |  |
| FrequencyAxisScaling | Linear |  |
| NoiseDynamic | 0.25 | dB |
| NoiseResolution | SymbolRate/4. | Hz |
| NoiseThreshold | -100 | dB |
| ActiveFilterBandwidth | 2*SampleRateDefault | Hz |
| ActiveFilterBandwidthCenterFrequency | TxCarrierFrequency | Hz |
| LimitToActiveFilterBandwidth | SampledBands |  |
| ConserveMemory | On |  |
| InternalVisualizersList | Magnitude, Phase, Delay, ImpulseResponse |  |
| MagnitudeTitle | Magnitude |  |
| PhaseTitle | Phase |  |
| DelayTitle | Delay |  |
| ImpulseResponseTitle | Impulse Response |  |
| VisualizationMode | None |  |
| Title |  |  |
| SaveToFile | Off |  |
| LogFilename |  |  |
| Magnitude | Off |  |
| Phase | Off |  |
| Delay | Off |  |
| ImpulseResponse | Off |  |
| Active | On |  |

---
## 7. Receptor Óptico Coerente — Rx_DualPol

Front-end óptico coerente de dupla polarização (oscilador local, híbrida 90°, fotodetectores balanceados, filtro elétrico e ADC).

| Parâmetro | Valor | Unidade |
|---|---|---|
| ModelType | Ideal |  |
| ExtinctionRatioX | 50 | dB |
| ExtinctionRatioY | 50 | dB |
| InsertionLossSignal_X | 0.0 | dB |
| InsertionLossLocalOscillator_X | 0.0 | dB |
| PhaseImbalance_SignalLocalOscillator_X | 0.0 | deg |
| PhaseImbalance_IQ_X | 0.0 | deg |
| InsertionLossImbalanceI_X | 0.0 | dB |
| InsertionLossImbalanceQ_X | 0.0 | dB |
| InsertionLossSignal_Y | 0.0 | dB |
| InsertionLossLocalOscillator_Y | 0.0 | dB |
| PhaseImbalance_SignalLocalOscillator_Y | 0.0 | deg |
| PhaseImbalance_IQ_Y | 0.0 | deg |
| InsertionLossImbalanceI_Y | 0.0 | dB |
| InsertionLossImbalanceQ_Y | 0.0 | dB |
| EmissionFrequency | TxCarrierFrequency | Hz |
| LO_Offset | 500e6 | Hz |
| SampleRate | SampleRateDefault | Hz |
| AveragePower | 1.0e-3 | W |
| Linewidth | 500e3 | Hz |
| LO_RandomNumberSeed | 0 |  |
| QuantizedEmissionFreq | Yes |  |
| SimulationBandDefinition | FreeParameter |  |
| SimulationBandCenterFreq | EmissionFrequency | Hz |
| Responsivity | 1.0 | A/W |
| DarkCurrent | 0.0 | A |
| ThermalNoise | 10.0e-12 | A/Hz^(1/2) |
| ShotNoise | On |  |
| PD_RandomNumberSeed | 123 |  |
| LPF_Active | On |  |
| TransferFunction | Bessel |  |
| Bandwidth | 1.15*SymbolRate | Hz |
| PassBandRipple | 5e-1 | dB |
| FilterOrderDefinition | FilterOrder |  |
| FilterOrder | 4 |  |
| GaussianOrder | 1 |  |
| StopBandwidth | 2*BitRateDefault | Hz |
| StopbandAttenuation | 40 | dB |
| FrequencyAxisScaling | Logarithmic |  |
| MinimumPhase | Off |  |
| DigitalFilterOrder | 64 |  |
| ConserveMemory | On |  |
| ClockActive | Off |  |
| DelayCalculation | Individual |  |
| CompensationType | Individual |  |
| Resynchronize | Off |  |
| ChannelLabel |  |  |
| ChannelLabelX |  |  |
| ChannelLabelY |  |  |
| ChannelIndex | -1 |  |
| SamplingRate | 2*SymbolRate | Hz |
| StartTime | (0.0/SamplingRate) | s |
| TransferFunctionADC | Bipolar |  |
| TransferFunctionFilenameADC |  |  |
| Resolution | 12 | bit |
| FSRdefinition | Relative |  |
| Range | 1.0 |  |
| Vmin | 0.0 |  |
| Vmax | 1.0 |  |
| HalfLSBOffset | Yes |  |
| RestoreSampleRate | No |  |
| Interpolation | No |  |
| RandomNumberSeed | 0 |  |
| SamplingJitter | Uniform |  |
| RMS | 0.0 |  |
| GainError | 0.0 |  |
| NonlinearityError | 0.0 |  |
| OffsetError | 0.0 |  |
| Active | On |  |


---
## 8. Cadeia de Processamento Digital de Sinais (DSP) do Receptor — DSP_Lib

A recepção coerente do Canal 2 passa por **6 blocos DSP em sequência** (todos instâncias do módulo genérico `DSP_Lib`, cada um configurado com um algoritmo diferente). Ordem funcional na cadeia de recepção:

1. **Compensação de Dispersão Cromática (CD)** — `CDcompensation`
2. **Equalização MIMO adaptativa (recuperação de polarização/ISI)** — `TDE_MIMO`
3. **Demultiplexação de polarização (Stokes Space)** — `StokesSpacePolDemux`
4. **Recuperação de relógio/temporização (Clock & Slice)** — `Clock_and_Slice`
5. **Compensação de desvio de frequência da portadora** — `FreqOffsetComp`
6. **Recuperação de fase da portadora (Blind Phase Search)** — `BlindPhaseSearch`

### 8.1 CD Compensation (Compensação de Dispersão Cromática)  *(instância: `DSP_Lib_vtmg10`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | EqualizationAndPolDemux |  |
| Algorithm | CDcompensation |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | Reset |  |
| Active | On |  |
| AlgTDE | CMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | No |  |
| Format | 16QAM |  |
| ConstMeanPower | 20.0 |  |
| RadiiSq | 2 3 4 6 8 |  |
| Dec | 1 |  |
| Decimate | Yes |  |
| Delay | 1 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0.08e-3 | s/m^3 |
| EffFFTSize | 1 #0.25 |  |
| FFTsize | 0 |  |
| FiberLength | !" expr {Length}" | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 299792458/TxCarrierFrequency | m |
| LambdaRef | 299792458/TxCarrierFrequency | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename |  |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 16 |  |
| NPreSymbols | 3 |  |
| Filtering | MovingAverage |  |
| Ns | 1024 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | Off |  |
| NTaps | 6 |  |
| hxx | (0, 0) |  |
| hxy | (0, 0) |  |
| hyx | (0, 0) |  |
| hyy | (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1 |  |
| PowerdBm | 0 | (dBm) |
| s | 1e-4 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 1 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |

### 8.2 TDE_MIMO (Equalização Adaptativa MIMO)  *(instância: `DSP_Lib_vtmg9`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | EqualizationAndPolDemux |  |
| Algorithm | TDE_MIMO |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | IndependentBlocks |  |
| Active | On |  |
| AlgTDE | MMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | No |  |
| Format | 16QAM |  |
| ConstMeanPower | 170 |  |
| RadiiSq | 2 10 18 26 34 50 58 74 82 90 98 106 122 130 146 162 170 178 194 202 218 226 234 242 250 274 290 306 338 346 394 450 |  |
| Dec | 1 |  |
| Decimate | No |  |
| Delay | 0 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0 | s/m^3 |
| EffFFTSize | 1 |  |
| FFTsize | 0 |  |
| FiberLength | 0 | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 155e-8 | m |
| LambdaRef | 155e-8 | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename |  |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 16 |  |
| NPreSymbols | 3 |  |
| Filtering | MovingAverage |  |
| Ns | 1024 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | RRC |  |
| NTaps | 7 |  |
| hxx | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hxy | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyx | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyy | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1000 |  |
| PowerdBm | 0 | (dBm) |
| s | 5e-4 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 20 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |

### 8.3 Stokes Space Pol. Demux (Demultiplexação de Polarização)  *(instância: `DSP_Lib_vtmg1`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | EqualizationAndPolDemux |  |
| Algorithm | StokesSpacePolDemux |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | IndependentBlocks |  |
| Active | On |  |
| AlgTDE | MMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | No |  |
| Format | 16QAM |  |
| ConstMeanPower | 170 |  |
| RadiiSq | 2 10 18 26 34 50 58 74 82 90 98 106 122 130 146 162 170 178 194 202 218 226 234 242 250 274 290 306 338 346 394 450 |  |
| Dec | 1 |  |
| Decimate | No |  |
| Delay | 0 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0 | s/m^3 |
| EffFFTSize | 1 |  |
| FFTsize | 0 |  |
| FiberLength | 0 | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 155e-8 | m |
| LambdaRef | 155e-8 | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename |  |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 16 |  |
| NPreSymbols | 3 |  |
| Filtering | MovingAverage |  |
| Ns | 1024 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | RRC |  |
| NTaps | 17 |  |
| hxx | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hxy | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyx | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyy | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1000 |  |
| PowerdBm | 0 | (dBm) |
| s | 5e-4 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 20 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |

### 8.4 Clock and Slice (Recuperação de Temporização)  *(instância: `DSP_Lib_vtmg8`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | TimingAndSampling |  |
| Algorithm | Clock_and_Slice |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | IndependentBlocks |  |
| Active | On |  |
| AlgTDE | MMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | No |  |
| Format | 64QAM |  |
| ConstMeanPower | 170 |  |
| RadiiSq | 2 10 18 26 34 50 58 74 82 90 98 106 122 130 146 162 170 178 194 202 218 226 234 242 250 274 290 306 338 346 394 450 |  |
| Dec | 1 |  |
| Decimate | Yes |  |
| Delay | 1 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0 | s/m^3 |
| EffFFTSize | 1 |  |
| FFTsize | 0 |  |
| FiberLength | 0 | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 155e-8 | m |
| LambdaRef | 155e-8 | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename |  |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 16 |  |
| NPreSymbols | 3 |  |
| Filtering | MovingAverage |  |
| Ns | 1024 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | Off |  |
| NTaps | 10 |  |
| hxx | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hxy | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyx | (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| hyy | (0, 0) (0, 0) (1, 0) (0, 0) (0, 0) (0, 0) (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1 |  |
| PowerdBm | 0 | (dBm) |
| s | 5e-5 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 10 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |

### 8.5 FreqOffsetComp (Compensação de Offset de Frequência)  *(instância: `DSP_Lib_vtmg2`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | PhaseAndFrequency |  |
| Algorithm | FreqOffsetComp |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | Reset |  |
| Active | On |  |
| AlgTDE | CMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | Yes |  |
| Format | ReadFromFile |  |
| ConstMeanPower | 20.0 |  |
| RadiiSq | 2 3 4 6 8 |  |
| Dec | 1 |  |
| Decimate | Yes |  |
| Delay | 1 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0 | s/m^3 |
| EffFFTSize | 1 |  |
| FFTsize | 0 |  |
| FiberLength | 0 | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 155e-8 | m |
| LambdaRef | 155e-8 | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename | 256QAM_constellation.txt |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 64 |  |
| NPreSymbols | 21 |  |
| Filtering | MovingAverage |  |
| Ns | 1024 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | Off |  |
| NTaps | 6 |  |
| hxx | (0, 0) |  |
| hxy | (0, 0) |  |
| hyx | (0, 0) |  |
| hyy | (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1 |  |
| PowerdBm | 0 | (dBm) |
| s | 1e-4 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 1 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |

### 8.6 Blind Phase Search (Recuperação de Fase)  *(instância: `DSP_Lib_vtmg7`)*

| Parâmetro | Valor | Unidade |
|---|---|---|
| AlgorithmType | PhaseAndFrequency |  |
| Algorithm | BlindPhaseSearch |  |
| BaudRate | SymbolRate | Bd |
| MultipleBlockMode | Reset |  |
| Active | On |  |
| AlgTDE | CMA |  |
| alpha | 0.2e-3 | (dB/m) |
| BitsPerSymbol | 1 |  |
| FreqShift | 1e9 | Hz |
| CombinedPol | Yes |  |
| Format | ReadFromFile |  |
| ConstMeanPower | 20.0 |  |
| RadiiSq | 2 3 4 6 8 |  |
| Dec | 1 |  |
| Decimate | Yes |  |
| Delay | 1 |  |
| ApplyAutoDelay | Yes |  |
| TypeOfFIR | Complex-valued |  |
| Dispersion | 16e-6 | s/m^2 |
| DispersionSlope | 0 | s/m^3 |
| EffFFTSize | 1 |  |
| FFTsize | 0 |  |
| FiberLength | 0 | m |
| FIR | (0, 0) (1, 0) (0, 0) |  |
| FIRreal | 0.0 1.0 0.0 |  |
| gamma | 1.3e-3 | (1/W-m) |
| InitialPhaseNoise | 0 | rad |
| Lambda | 155e-8 | m |
| LambdaRef | 155e-8 | m |
| Mode | Decision-Directed |  |
| Constellation | 16QAM |  |
| NbSpans | 1 |  |
| ConstellationFilename | 256QAM_constellation.txt |  |
| RotationalSymmetry | 4 |  |
| SignalType | PolMux |  |
| ChannelLabel |  |  |
| LogicalChannels |  |  |
| ChannelIndex | -1 |  |
| ChannelLabelX |  |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY |  |  |
| ChannelIndexY | -1 |  |
| NbTestAngles | 64 |  |
| NPreSymbols | 21 |  |
| Filtering | MovingAverage |  |
| Ns | 0 |  |
| FourthPower | Yes |  |
| Nsteps | 10 |  |
| InitTaps | Off |  |
| NTaps | 6 |  |
| hxx | (0, 0) |  |
| hxy | (0, 0) |  |
| hyx | (0, 0) |  |
| hyy | (0, 0) |  |
| ReadTapsFilename |  |  |
| PilotsAveraging | 4 |  |
| FineBlindEstimation | Yes |  |
| FilterLength | 32 |  |
| JointEstimation | No |  |
| PlotEstimatedPhase | No |  |
| PlotEstimatedError | No |  |
| ErrorSmoothingWindow | 1 |  |
| PowerdBm | 0 | (dBm) |
| s | 1e-4 |  |
| s2 | 0.02 |  |
| ResamplingMethod | FrequencyDomain |  |
| KaiserWinBeta | 5.0 |  |
| SamplingRateOut | 10e9 |  |
| res | 1 |  |
| StepType | log |  |
| RollOff | 0.2 |  |
| NyquistType | raised cosine |  |
| ScramblingLength | 64 |  |
| factor | 2 |  |
| phase | 0 |  |
| ClipRatio | 0.01 |  |
| Vpi | 5 |  |
| linearGain | 1 |  |
| ExportTaps | Off |  |
| WriteTapsFilename |  |  |
| NormalizationMode | Blind |  |
| VolterraMode | ComponentModeling |  |
| CorrectForPilotSymbols | No |  |
| VolterraType | SymmetricKernels |  |
| KernelFilenameToRead1 | kernels1.txt |  |
| KernelFilenameToRead2 | kernels2.txt |  |
| KernelFilenameToRead3 | kernels3.txt |  |
| KernelFilenameToRead4 | kernels4.txt |  |
| MEMs | 50 9 9 |  |
| SignalsSamplingRate | 64e9 | Hz |
| TrainingInputFilename1 | input1.txt |  |
| TrainingInputFilename2 | input2.txt |  |
| TrainingInputFilename3 | input3.txt |  |
| TrainingInputFilename4 | input4.txt |  |
| TrainingOutputFilename1 | output1.txt |  |
| TrainingOutputFilename2 | output2.txt |  |
| TrainingOutputFilename3 | output3.txt |  |
| TrainingOutputFilename4 | output4.txt |  |
| PlotKernels | No |  |
| PlotSignals | Yes |  |
| SaveKernels | Yes |  |
| KernelBaseFilenameToWrite | kernels |  |
| OperationMode | Align |  |
| ScalingFactor | 2 |  |
| View | 30 45 |  |
| SignalFraction | 1 |  |
| SignalFractionToPlot | 0.25 |  |
| LimitSignalRange | On |  |
| SignalRangeValues | 0 10 |  |
| PlotSamplesPerSymbol | 1 |  |
| PlotTitle |  |  |
| DecimatedPlot | Yes |  |
| SampleTime | 0.5 |  |
| SampleRange | 0 |  |
| ClockRecovery | Off |  |
| ConstellationSize | 16 |  |
| InputSamplesPerSymbol | 2 |  |
| OutputSamplesPerSymbol | 1 |  |
| N | 1024 |  |
| SamplingOffset | 0 |  |
| Plot | Off |  |
| FilterType | Baseband |  |
| CenterFrequency | 1e9 | Hz |
| Fc | 2e9 | Hz |
| FFTWindowSize | 128 |  |
| NumberOfCarriers | 1 |  |
| CarrierBandwidth | SymbolRate | Hz |
| CarrierSpacing | 0 |  |
| PlotSpectrum | No |  |
| FilterControl | Dynamic |  |
| FF_coefficients | 0 0 0 1 0 0 0 |  |
| FB_coefficients | 0 1 0 |  |
| ReferenceSignalFileName |  |  |
| Iterations | 1 |  |
| StepSize | 1e-3 |  |
| DecimationPhase | 0 |  |
| Thresholds | 0 |  |
| Levels | -0.5 0.5 |  |
| FFCoefficientsFile | FFCoeffs.vsdb |  |
| FBCoefficientsFile | FBCoeffs.vsdb |  |
| Debiasing | On |  |
| Normalization | Internal |  |
| PlotConvergence | No |  |
| PlotCoefficients | No |  |
| InterfaceType | Python |  |
| LogicalInformation | No |  |
| Path | <self>/ |  |
| UserIfcFile | ExampleUserLib.py |  |
| InitCommand |  |  |
| UserAlgorithm | Example_DSP |  |
| UserParams | 3.0 #Used in Example Procedure |  |
| WrapupCommand |  |  |
| ShareInterface | On |  |


---
## 9. Sincronização — SynchDSP

Sincroniza os dados recebidos com os dados lógicos transmitidos (alinhamento de atraso/polarização) para a estimativa de BER.

| Parâmetro | Valor | Unidade |
|---|---|---|
| Mode | PolMux |  |
| DelayCalculation | Combined |  |
| CompensationType | Individual |  |
| Resynchronize | On |  |
| RepeatLogicalBitSeq | No |  |
| SamplesPerSymbol | 1 |  |
| AutoSwapPolarization | Yes |  |
| ChannelLabel | ChX |  |
| ChannelIndex | -1 |  |
| ChannelLabelX | Ch_X2 |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY | Ch_Y2 |  |
| ChannelIndexY | -1 |  |
| Active | On |  |

---
## 10. Estimador de Taxa de Erro de Bit — BER_4D

Calcula a BER (Bit Error Rate) por contagem Monte Carlo e gera os diagramas de constelação X/Y.

| Parâmetro | Valor | Unidade |
|---|---|---|
| Mode | PolMux |  |
| ConstellationAlignment | On |  |
| DetailedReport | No |  |
| RealNumberFormat |  |  |
| OutputFilename |  |  |
| FileOutputMode | Overwrite |  |
| Outputs | BER_MC |  |
| SignalOutput | Bits |  |
| ChannelLabel | ChX |  |
| ChannelIndex | -1 |  |
| ChannelLabelX | Ch_X2 |  |
| ChannelIndexX | -1 |  |
| ChannelLabelY | Ch_Y2 |  |
| ChannelIndexY | -1 |  |
| RepeatLogicalBitSeq | No |  |
| IgnoreSymbols | ApplyOnce |  |
| StartSymbolToIgnore | 0#s_NumberOfSymbols - 100 |  |
| NumberOfSymbolsToIgnore | 0#300 |  |
| MultipleBlockMode | IndependentBlocks |  |
| InternalVisualizersList | Constellation |  |
| Constellation | Constellation |  |
| VisualizationMode | StartAutomatically |  |
| Title | Constellation Diagram |  |
| TitleX | Constellation Xpol. |  |
| TitleY | Constellation Ypol. |  |
| Active | On |  |

---
## 11. Outros Componentes Auxiliares

### 11.1 Resample (Reamostragem do sinal óptico)

| Parâmetro | Valor | Unidade |
|---|---|---|
| SampledSignalMode | Absolute |  |
| CenterFrequency | TxCarrierFrequency | Hz |
| SampleRate | (4*75e9)*2 | Hz |
| MaximumDynamicRange | 50.0 | dB |
| ResampleFactor | 1 |  |
| NoiseBinMode | Skip |  |
| NoiseBandwidth | 10e12 | Hz |
| NoiseCenterFrequency | 193.1e12 | Hz |
| NoiseBinThreshold | 1e-19 | W/Hz |
| FreqResolutionNB | 100e9 | Hz |
| GridReferenceFrequency | 193.1e12 | Hz |
| NoiseBinDynamic | 3.0 | dB |
| PolarizationTolerance | 1.0 |  |
| GreatestPrimeFactorLimit | GreatestPrimeFactorLimit |  |
| IntegerSamplesPerSymbol | Off |  |
| SymbolRate | BitRateDefault | Bd |
| Active | On |  |

### 11.2 Gain (Ganho — utilizado no ramo de referência/monitoramento)

| Parâmetro | Valor | Unidade |
|---|---|---|
| gain | 0.5 |  |

### 11.3 Const (Fonte constante — referência de potência PowerdBm)

| Parâmetro | Valor | Unidade |
|---|---|---|
| level | PowerdBm |  |

### 11.4 SignalAnalyzer (Analisador de espectro/sinal) — instância 1

| Parâmetro | Valor | Unidade |
|---|---|---|
| SettingsSource | AnalyzerFile |  |
| AnalyzerMode | StartAutomatically |  |
| InitialAnalysisType | Spectrum |  |
| Title |  |  |
| Settings |  |  |
| SettingsFile |  |  |
| InternalVisualizersList |  |  |
| Active | On |  |

### 11.5 SignalAnalyzer — instância 2 (desativado — Active = Off)

| Parâmetro | Valor | Unidade |
|---|---|---|
| SettingsSource | AnalyzerFile |  |
| AnalyzerMode | StartAutomatically |  |
| InitialAnalysisType | Spectrum |  |
| Title |  |  |
| Settings |  |  |
| SettingsFile |  |  |
| InternalVisualizersList |  |  |
| Active | Off |  |

### 11.6 NumericalAnalyzer2D (Analisador numérico 2D — plots de BER/OSNR etc.)

| Parâmetro | Valor | Unidade |
|---|---|---|
| SettingsSource | AnalyzerFile |  |
| AnalyzerMode | StartAutomatically |  |
| InitialViewType | Plot |  |
| Title |  |  |
| UpdateMode | EveryPoint |  |
| UpdatePeriod | 50 |  |
| Settings |  |  |
| SettingsFile |  |  |
| InternalVisualizersList |  |  |
| Active | On |  |

### 11.7 Ferramentas de fiação (sem parâmetros configuráveis)

| Componente | Função |  |
|---|---|---|
| Ground_vtms3 | Terra / referência de sinal não utilizado |  |
| Ground_vtms4 | Terra / referência de sinal não utilizado |  |
| Add_vtms1 | Somador de sinais elétricos/ópticos |  |
| BusCreate4_1_vtms1 | Combina 4 sinais em um barramento (bus) — usado para agrupar os 4 canais Tx |  |
