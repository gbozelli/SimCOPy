# Second Class Notes

Sinais diferentes para diferentes representações (múltiplas representações numéricas do mesmo sinal).

**Modulações:** 4QAM, QPSK, 16QAM, 32QAM, 64QAM.

A variável `simulation_window` (janela de simulação) é necessária pois, se definirmos um número total de símbolos simulados, o `bitrate` dos canais pode ser diferente, o que atrapalha essa configuração. `simulation_window` é o tempo total de simulação, e o `time_step` é derivado da `sampling_rate`. Se houver conflito entre `bitrate` e `sampling_rate` (o que causa aliasing), será necessário verificar que a `sampling_rate` seja suficiente para atender ao critério de Nyquist dos canais usando o `bitrate`.

```mermaid
flowchart TD
    System --> WDM-Number_of_channels
    System --> SamplingRate
    System --> SimulationWindow
    WDM-Number_of_channels --> Channel1
    WDM-Number_of_channels --> Channel2
    WDM-Number_of_channels --> ChannelN
    Channel1 --> IMDD
    IMDD --> DirectModulation
    IMDD --> OOK
    IMDD --> PAM4
    IMDD --> ExternalModulation
    Channel1 --> DCS
    Channel1 --> SEED
    DCS --> SP
    DCS --> QPSK
    DCS --> 16QAM
    DCS --> 32QAM
    DCS --> 64QAM
    DCS --> DP
    Channel1 --> BitRate
```

A partir disso, criar o transmissor, que define quantos bits serão gerados (por isso precisamos das configurações globais). O diagrama está a seguir:

```mermaid
flowchart TD
    WDMChannel_parameters --> Transmitter
    WDMChannel_parameters --> laser_power_output
    WDMChannel_parameters --> insertion_loss_mach_zehnder_if_external_mod
    WDMChannel_parameters --> extinction_ratio_mach_zehnder
    WDMChannel_parameters --> largura_linha_laser
    WDMChannel_parameters --> other_parameters
    laser_power_output --> Transmitter
    insertion_loss_mach_zehnder_if_external_mod --> Transmitter
    extinction_ratio_mach_zehnder --> Transmitter
    largura_linha_laser --> Transmitter
    other_parameters --> Transmitter
    WDMChannel_parameters --> Channel
    Transmitter --> bit_sequence_generator_bits
    Transmitter --> symbol_sequence_generation_symbols
    bit_sequence_generator_bits --> symbol_sequence_generation_symbols
    symbol_sequence_generation_symbols --> pulse_shaping_electrical
    Transmitter --> pulse_shaping_electrical
    Transmitter --> electrical_to_optical_conv
    Transmitter --> directmodulation
    Transmitter --> externalmodulation
    Transmitter --> powercontrol
    Transmitter --> multiplexer
    pulse_shaping_electrical --> electrical_to_optical_conv
    Transmitter --> electrical_to_optical_conv_optical
    electrical_to_optical_conv --> directmodulation
    electrical_to_optical_conv --> externalmodulation
    directmodulation --> powercontrol
    externalmodulation --> powercontrol
    powercontrol --> multiplexer
    multiplexer --> Link

    WDMChannel_parameters --> all_fiber_parameters
    WDMChannel_parameters --> line_compensation
    line_compensation --> Link
    all_fiber_parameters --> Link
    Link --> SplitStepFourierMethod
    SplitStepFourierMethod --> Demux
    Demux --> Receptor
    WDMChannel_parameters --> frequencia_central
    WDMChannel_parameters --> frequencia_nominal
    WDMChannel_parameters --> taxa_de_transmissao
    WDMChannel_parameters --> formato_de_mod
    WDMChannel_parameters --> tudo_necessario_para_recuperar_sinal
    WDMChannel_parameters --> other_parameters
    other_parameters --> Receptor
    frequencia_central --> Receptor
    frequencia_nominal --> Receptor
    taxa_de_transmissao --> Receptor
    formato_de_mod --> Receptor
    tudo_necessario_para_recuperar_sinal --> Receptor
```

O `WDMChannel_parameters` serve para TX, Channel e RX através dos parâmetros definidos no primeiro diagrama.

O `Link` precisa definir se é SP (Single Polarization) ou DP (Dual Polarization). Observação: o PMD (Polarization Mode Dispersion) não está implementado, e a rotação de polarização não está correta (abrange apenas um hemisfério da esfera de Poincaré). O diagrama não explicita essas condicionais (tem um "if" escrito em algum local no código).

Considerar tudo como processos (dispositivo é um processo, algo é um processo, talvez verificar como `dataclass` pode organizar isso).

Assumir que o sistema é coerente e trabalhar com IMDD.