# De onde vem cada coisa no simcopy: implementações e decisões (2026-11-19)

Este documento serve para conferir e validar o projeto. Para cada implementação e cada decisão ele
diz em que se baseou e onde conferir. Ele não repete o conteúdo dos guias: aponta para eles.

Como ler:

- Fontes, pela ordem de prioridade do projeto: [A] o Agrawal, depois [P] o código do professor,
  depois [V] a documentação do VPI e [L] os artigos. [G] é uma decisão do Gabriel (com o número da
  dúvida em `docs/duvidas_abertas.md`), e [M] é uma medida nossa (com o script e o resultado).
- Grau: Fato (está na fonte citada), Inferência (dedução nossa, falta confirmar) ou Escolha
  (decisão de projeto sem uma fonte única que a imponha).
- As decisões de arquitetura com identificador (D-1 a D-28) estão em `docs/decisoes.md`, que é a
  fonte única do que está decidido. Aqui fica a origem de cada uma.

## 1. Bibliografia

### 1.1 Livros

| Sigla | Referência | Uso |
|---|---|---|
| [A] | G. P. Agrawal, *Sistemas de Comunicação por Fibra Óptica* (tradução de *Fiber-Optic Communication Systems*). As equações citadas no código seguem a numeração do exemplar usado no projeto. Conferir a edição: as seções 10.1 a 10.3 (formatos avançados e receptor coerente) indicam a 4ª edição. | Fonte principal de toda a física |
| [N] | G. P. Agrawal, *Nonlinear Fiber Optics* (Academic Press) | Pulsos gaussianos com chirp (seção 3.2) e largura rms com β₂ e β₃ (seção 3.3), nos testes de física da fibra |

### 1.2 Código e material da disciplina

| Sigla | Referência | Uso |
|---|---|---|
| [P] | Código do professor: `draft.py` e os notebooks da disciplina (por exemplo `4_EDFAsAndReceivers`). Funções citadas no código: `QAM_mod`, `QAM_dem`, `QAM_receiver`, `DAC_Nyquist`, `HybridNetwork`, `DPQAM_transmitter`, `DPQAM_receiver`, `EDFA` | Segunda fonte, depois do Agrawal |
| [Aula] | Notas de aula (linhas E11 a E13 sobre o coerente, "IMDD < 5w", filtro de Bessel no IMDD) | Dúvidas 42, 50 e 55 |

### 1.3 Referência de comparação (VPI)

| Sigla | Referência | Uso |
|---|---|---|
| [V] | Documento de parâmetros do esquemático `WDM_of_400G_ZR_4_Channels_Ch2` (VPI TransmissionMaker 11.5): `reference/vpi/parametros_vtmu.md` | Parâmetros do cenário 400ZR, da cadeia de DSP e dos dispositivos |
| [V-dados] | Constelações exportadas pelo grupo (canal 2, 2026, 400gbps_CH2, 2025, 2021), fora do repositório | Comparações e hipóteses das bordas. Resumo em `reference/vpi/catalogo/catalogo_vpi.json` |
| [OIF] | OIF, *400ZR Implementation Agreement* (OIF-400ZR-01.0, 2020) | Contexto do cenário (DP-16QAM a 59,84 GBd). Não é usado no código |

### 1.4 Artigos

DSP coerente:

| Sigla | Referência | Uso |
|---|---|---|
| [Go80] | D. N. Godard, "Self-recovering equalization and carrier tracking in two-dimensional data communication systems", *IEEE Trans. Commun.* 28(11), 1867-1875, 1980 | Equalizador de módulo constante (CMA) |
| [Fa09] | I. Fatadin, D. Ives, S. J. Savory, "Blind equalization and carrier phase recovery in a 16-QAM optical coherent system", *J. Lightw. Technol.* 27(15), 3042-3049, 2009 | Equalizador dirigido por raio (RDE/MMA) para 16QAM |
| [Sz10] | B. Szafraniec, B. Nebendahl, T. Marshall, "Polarization demultiplexing in Stokes space", *Opt. Express* 18(17), 17928-17939, 2010 | Separação de polarização no espaço de Stokes |
| [Pf09] | T. Pfau, S. Hoffmann, R. Noé, "Hardware-efficient coherent digital receiver concept with feedforward carrier recovery for M-QAM constellations", *J. Lightw. Technol.* 27(8), 989-999, 2009 | Busca cega de fase (BPS) |
| [VV83] | A. J. Viterbi, A. M. Viterbi, "Nonlinear estimation of PSK-modulated carrier phase with application to burst digital transmission", *IEEE Trans. Inf. Theory* 29(4), 543-551, 1983 | Estimativa de fase pela quarta potência (hipóteses 1 e 4) |
| [Fa10] | I. Fatadin, D. Ives, S. J. Savory, "Laser linewidth tolerance for 16-QAM coherent optical systems using QPSK partitioning", *IEEE Photon. Technol. Lett.* 22(9), 631-633, 2010 | Quarta potência só com os anéis de dentro e de fora (hipóteses 1 e 4) |
| [MM10] | R. Muller, D. A. A. Mello, "Phase-offset estimation for joint-polarization phase-recovery in DP-16-QAM systems", *IEEE Photon. Technol. Lett.* 22(20), 1515-1517, 2010 | Fase estimada em conjunto nas duas polarizações (hipótese 4) |
| [Zh11] | Q. Zhuge, C. Chen, D. V. Plant, "Low computation complexity two-stage feedforward carrier recovery algorithm for M-QAM", OFC 2011 | Busca de fase em dois estágios (hipótese 1) |
| [Sa08] | S. J. Savory, "Digital filters for coherent optical receivers", *Opt. Express* 16(2), 804-817, 2008 | Compensação de dispersão no DSP e equalizador em borboleta 2×2 (próximos passos: MIMO e PMD) |
| [IK07] | E. Ip, J. M. Kahn, "Digital equalization of chromatic dispersion and polarization mode dispersion", *J. Lightw. Technol.* 25(8), 2033-2043, 2007 | Compensação de PMD pelo equalizador (próximos passos) |
| [Fa08] | I. Fatadin, S. J. Savory, D. Ives, "Compensation of quadrature imbalance in an optical QPSK coherent receiver", *IEEE Photon. Technol. Lett.* 20(20), 1733-1735, 2008 | Ortogonalização de Gram-Schmidt (GSOP) (próximos passos) |
| [Sa10] | S. J. Savory, "Digital coherent optical receivers: algorithms and subsystems", *IEEE J. Sel. Top. Quantum Electron.* 16(5), 1164-1179, 2010 | Visão geral da cadeia de DSP e estimativa de frequência pela quarta potência. Tutorial recomendado |
| [Ki16] | K. Kikuchi, "Fundamentals of coherent optical fiber communications", *J. Lightw. Technol.* 34(1), 157-179, 2016 | Tutorial do receptor coerente |
| [FS17] | M. S. Faruk, S. J. Savory, "Digital signal processing for coherent transceivers employing multilevel formats", *J. Lightw. Technol.* 35(5), 1125-1141, 2017 | Tutorial de DSP para 16QAM e 64QAM |
| [Al15] | A. Alvarado, E. Agrell, D. Lavery, R. Maher, P. Bayvel, "Replacing the soft-decision FEC limit paradigm in the design of optical communication systems", *J. Lightw. Technol.* 33(20), 4338-4352, 2015 | MI e GMI |

Fibra e PMD:

| Sigla | Referência | Uso |
|---|---|---|
| [MMW97] | D. Marcuse, C. R. Menyuk, P. K. A. Wai, "Application of the Manakov-PMD equation to studies of signal propagation in optical fibers with randomly varying birefringence", *J. Lightw. Technol.* 15(9), 1735-1746, 1997 | Fator 8/9 da equação de Manakov |
| [GK00] | J. P. Gordon, H. Kogelnik, "PMD fundamentals: polarization mode dispersion in optical fibers", *PNAS* 97(9), 4541-4550, 2000 | Fundamentos de PMD (DGD maxwelliana) |
| [Ba25] | Banner et al., BIFROST, *Phys. Rev. Applied*, arXiv:2510.01212 | Sorteio uniforme de matrizes SU(2) por quatérnios (eq. 17). O BIFROST é GPL: só foi lido, nada importado (D-19) |

### 1.5 Outros simuladores (consultados, ver `docs/comparacao_simuladores.md`)

OptiCommPy (GPL), QAMpy (GPL), pyofss (GPL), Robochameleon, OpticalCommunications, comnum (MIT/Apache)
e maiman-studio (Apache), na pasta `outros_simuladores/`, fora do repositório. Os GPL nunca são
importados nem copiados. O resultado da busca por uma solução para a fase do VPI está no diário
(`docs/replicacao_vpi_diario.md`, "Busca nos outros simuladores").

### 1.6 Figuras

[GF] Guia de figuras do TCC do Gabriel (`TCC/codigos/guia_figuras.md`, fora do repositório): Times,
ticks para dentro, sem título, legenda acima sem moldura, vírgula decimal, PDF com fonte embutida.

## 2. Implementações, módulo a módulo

"Teste" diz se o módulo tem teste unitário hoje. Os que não têm estão nas issues de teste.

### 2.1 Núcleo, configuração e validação

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `core/units.py` | constantes e conversões (dB, dBm, β₂ e β₃ a partir de D e S, γ a partir de n₂, comprimento efetivo) | [A] Eq. 2.3.13 (β₃) e cap. 2 | Fato | não (issue #5) |
| `core/grid.py` | grade de tempo e frequência compartilhada; frequência óptica de cada bin | [A] Eq. 2.4.3 e 2.4.5, convenção vermelha [G] dúvida 1 | Fato + decisão | não (issue #6) |
| `core/signals.py` | tipos `Bits`, `Symbols`, `Waveform`, `DigitalSignal` | arquitetura, seção 4.2.1 (D-13) | Escolha | sim, parcial (issue #7) |
| `core/context.py` | sementes por canal e por ruído, sem `np.random` global | arquitetura; [G] dúvidas 8 e 49 | Escolha | sim, parcial (issue #6) |
| `config/system.py` | configuração do sistema (uniões discriminadas IMDD e coerente) | arquitetura, seção 4.6 | Escolha | sim, parcial (issue #8) |
| `config/resolve.py` | grandezas derivadas (parâmetros do laser, fibra, guardas) | [A] Eq. 3.3.4, 3.3.6, 3.3.10 e texto após 3.3.3 (R_sp) | Fato | sim, parcial (issue #8) |
| `config/io.py` | leitura de YAML e JSON | Escolha | Escolha | sim, parcial (issue #9) |
| `validation/rules.py` | invariantes antes de alocar (Nyquist, amostras por símbolo inteiras, sobreposição de canais, passo do SSFM, MZM no nulo no coerente, aviso de Brillouin) | arquitetura PR-3; [A] seções 2.6.1 e 6.3 (Brillouin, [G] dúvida 42); [A] seção 10.1.2 e [P] `DPQAM_transmitter` (nulo, [G] dúvida 10); [V] `GreatestPrimeFactorLimit` | Fato + Escolha | sim, parcial (issue #9) |

### 2.2 Fibra

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/ssfm.py` | split-step de Fourier simétrico, escalar e vetorial (Manakov), passo adaptativo pela fase não linear, β₃ | [A] cap. 2 e Eq. 2.4.4 (β₃); [MMW97] e [V] `NonlinearAdjustmentFactors` (8/9, D-2, D-3); [V] passo adaptativo pela fase não linear | Fato | testes de física 01 a 05 (atenuação, dispersão, não linearidade, PMD, sóliton) |
| `kernels/birefringence.py` | PMD de primeira ordem por seções: rotação SU(2) uniforme e DGD por passo | [Ba25] (sorteio), [GK00]; [V] `PMDCoefficient` 0,1 ps/√km e `CorrelationLength` 50 m (D-17, D-18) | Fato | teste de física 04; issue #10 |
| `components/optical/fiber.py` | bloco da fibra; passo máximo com PMD; campo ao longo de z; aviso com uma polarização | [G] dúvidas 32, 33 e 34 | Escolha | não |

### 2.3 Laser

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/laser/rate_equations.py` | equações de taxa em densidade, estado estacionário, forças de Langevin | [A] Eq. 3.3.1, 3.3.2, 3.3.15, 3.3.16, 3.3.25 a 3.3.28; forma em densidade do [P]; limite com poucos fótons [G] dúvida 7 | Fato | não (issue #14) |
| `kernels/laser/small_signal.py` | resposta de pequeno sinal e banda de 3 dB | [A] Eq. 3.3.17 a 3.3.23 | Fato | não (issue #14) |
| `kernels/laser/cw.py` | laser CW com ruído de fase (passeio de Wiener) | [A] Eq. 3.3.33 a 3.3.37 | Fato | não (issue #14) |
| `kernels/laser/conversions.py` | conversão entre a forma em número ([A]) e em densidade ([P]) | [A] texto após Eq. 3.3.3 e Eq. 3.3.8 | Fato | não (issue #14) |
| `components/tx/laser.py` | lasers CW e de modulação direta como blocos | [A] Eq. 3.3.14 e 3.4.1; [G] dúvidas 4 e 5 | Fato + decisão | não (issue #14) |

### 2.4 Moduladores

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/modulator/mzm.py` | MZM com braços de amplitudes diferentes, push-pull, modulador IQ, pré-distorção arco-seno | [A] seção 10.1.2, Eq. 10.1.3 a 10.1.5; braços desiguais do [P] ([G] dúvida 12); pré-distorção [V] `PreDistortion = Yes` ([G] dúvida 14) | Fato | issue #15 |
| `components/tx/modulator.py` | laser CW + MZM do IMDD | [A] Eq. 10.1.3; [G] dúvida 5 | Fato | não |
| `components/tx/iq_modulator.py` | modulador IQ (dois MZM no nulo, defasagem de 90°), divisão em X e Y, skew, desbalanço de fase e de ganho, banda S21 | [A] seção 10.1.2 e Fig. 10.3; [V] e [P] para a divisão em X e Y ([G] dúvida 16); skew ([G] dúvida 15); desbalanço de ganho ([G] dúvida 18); S21 Bessel de ordem 3 ([G] dúvida 13) | Fato + decisão | não |

### 2.5 Elétrica, IMDD e decisão

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/electrical.py` | NRZ e níveis do PAM | [Aula] (Bessel, não Nyquist, [G] dúvida 50); PAM4 ([G] dúvida 51) | Fato + decisão | não |
| `kernels/filters.py` | filtros passa-baixas na FFT (Bessel, retangular), banda equivalente de ruído | [A] Eq. 4.4.4 | Fato | não |
| `kernels/decision.py` | limiar ótimo, fator Q, BER do OOK e do PAM | [A] Eq. 4.6.6 a 4.6.11. A generalização para M níveis não está no Agrawal | Fato (OOK), Escolha (PAM) | não |
| `components/tx/bits.py`, `components/tx/driver.py` | fonte de bits pela semente do canal; acionamento NRZ filtrado | [G] dúvida 49 (bits aleatórios, não PRBS); [Aula] | Decisão | não |

### 2.6 Amplificador, fotodetector e WDM

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/amplifier.py` | ganho e ASE gaussiana circular; n_sp a partir da figura de ruído | [A] Eq. 7.2.11, 7.2.15 (invertida, [G] dúvida 19), 7.4.1 e 7.5.10 (ASE nas duas polarizações, [G] dúvida 20); [P] `EDFA`; ruído circular (C-2 da arquitetura) | Fato | issue #19 |
| `components/optical/amplifier.py` | EDFA com saturação, ganho ou potência alvo por canal, banda do sinal | [A] seção 7.2.2 e Exercício 7.2; [G] dúvidas 3, 21, 22 e 24 | Fato + decisão | não |
| `components/optical/power.py` | controle de potência de lançamento | arquitetura, anexo A.4 | Escolha | não |
| `kernels/photodetector.py` | fotocorrente, ruído de disparo e térmico | [A] Eq. 4.1.1, 4.1.3, 4.4.2 a 4.4.9; [P] notebook 4; padrão de 5,76 pA/√Hz ([G] dúvida 26) | Fato | issue #20 |
| `components/rx/photodiode.py` | fotodiodo como bloco | [A] seção 4.4; [G] dúvidas 6 e 28 | Fato | não |
| `kernels/wdm.py` | deslocamento em frequência, soma dos canais, filtro óptico retangular ou supergaussiano | [P] (soma e máscara retangular); [A] seção 6.3.1 (Eq. 6.3.1); [V] `WDM_MUX_N_1` (supergaussiano) | Fato | issue #21 |
| `components/optical/wdm.py` | mux e demux como blocos | [G] dúvidas 2, 29, 30 e 31 | Decisão | não |

### 2.7 Sistema coerente

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `kernels/mapping.py` | mapeamento e decisão de BPSK, QPSK, 16QAM e 64QAM, Gray por eixo | [A] seção 10.1; [P] `QAM_dem`. O 64QAM estende a regra do 16QAM ([G] dúvida 69) | Fato (até 16QAM), Escolha (64QAM) | não |
| `kernels/pulse_shaping.py` | cosseno levantado e sua raiz, aplicados na FFT | [P] `DAC_Nyquist` (cosseno levantado); [V] `NyquistResponse = squareRootRaisedCosine` (raiz) | Fato | não |
| `components/tx/iq_driver.py` | bits, símbolos, I(t) + iQ(t) | [P], [V]; excursão `drive_depth` ([A] seção 10.1.2, [G] dúvida 56) | Fato + decisão | não |
| `kernels/coherent.py` | híbrida de 90° e detecção balanceada | [A] Eq. 4.5.3 e 10.2.4, seção 8.7.3 e Fig. 8.39; [P] `HybridNetwork` | Fato | não |
| `components/rx/coherent.py` | oscilador local (livre ou cópia do laser), diversidade de polarização, ruído térmico em cada fotodiodo | [A] seções 4.5.2 e 8.7.3; [P] `DPQAM_receiver` (`X_fraction = 0.5`); [G] dúvidas 54 e 57 | Fato + decisão | não |
| `kernels/adc.py`, `components/rx/adc.py` | amostragem a 2 amostras por símbolo e quantização de 12 bits; receptor cego como bloco | [V] `Rx_DualPol` (2 amostras por símbolo, 12 bits); fundo de escala 4,0 ([G] dúvida 65, em aberto) | Fato + Escolha | não |
| `kernels/dsp.py`: compensação de dispersão | filtro na FFT com os valores verdadeiros do enlace (D-7) | [A] seção 8.7.3; [Sa08] | Fato | não |
| `kernels/dsp.py`: desvio de frequência | pico do espectro da quarta potência, com refinamento parabólico | [Sa10] | Fato | não |
| `kernels/dsp.py`: com símbolos conhecidos | atraso e frequência por correlação, equalizador por mínimos quadrados, fase e polarização em blocos | roteiro antigo, tarefa 6 ("data-aided" primeiro); [P] `QAM_receiver` (correlação) | Escolha | não |
| `kernels/dsp.py`: separação no espaço de Stokes | primeira estimativa da matriz de polarização | [Sz10]; [V] `StokesSpacePolDemux` | Fato | não |
| `kernels/dsp.py`: equalizador adaptativo cego | CMA e depois ajuste dirigido por decisão em bloco (mínimos quadrados), 7 coeficientes a 2 amostras por símbolo | [Go80]; [Fa09] (dirigido por raio, que vem desligado por padrão); [V] `TDE_MIMO` (7 coeficientes). A ordem e os padrões estão em aberto ([G] dúvida 66) | Fato + Escolha | não |
| `kernels/dsp.py`: busca cega de fase | BPS com 64 ângulos, janela de 43 símbolos, conjunta nas duas polarizações | [Pf09]; [V] `NbTestAngles = 64`, `NPreSymbols = 21`; janela de 43 medida nos dados do VPI ([M], diário, etapa A) | Fato + medida | não |
| `components/rx/dsp.py` | compensador de dispersão, compensador de frequência e filtro casado como blocos | idem | Fato | não |
| `kernels/demodulation.py` | sincronismo, alinhamento com os dados (polarização, atraso, giro), EVM, BER contada, fator Q, saltos de fase | [A] Eq. 10.3.2 e 10.3.7 (fator Q); [P] `QAM_receiver`; [V] `SynchDSP` (alinhamento). A BER aproximada do QAM quadrado por vizinho mais próximo não está no Agrawal | Fato + Escolha | não |
| `kernels/information.py` | MI e GMI | [Al15] (não está no Agrawal) | Fato | não |

### 2.8 Montagem, execução e gráficos

| Módulo | O que implementa | Base | Grau | Teste |
|---|---|---|---|---|
| `pipelines/imdd.py`, `pipelines/coherent.py` | montagem dos blocos a partir da configuração | arquitetura, seção 6.2 (escolhas na montagem, nunca `if` dentro de um bloco) | Escolha | não |
| `pipelines/runner.py`, `cli.py` | carregar, resolver, validar, montar, rodar, medir; varreduras | arquitetura, seção 6.3 | Escolha | não |
| `plotting.py`, `analysis.py` | estilo das figuras, mapas de densidade das constelações, gráficos de etapa, olho, varredura | [GF]; [A] Eq. 3.4.1 (chirp) | Escolha | não |

### 2.9 Scripts fora do pacote

Ficam em `scripts/`, não em `src/`, e não fazem parte da API.

| Script | O que faz | Base |
|---|---|---|
| `scripts/compara_vpi.py` | simula o cenário do VPI e compara SNR, BER e constelações; `--save-signals` guarda o sinal na saída do ADC | [V], [V-dados] |
| `scripts/catalogo_vpi.py` | lê todos os conjuntos do VPI e mede BER, SER, SNR e fração girada | [V-dados]; critério de fase perdida: mais de 2 % do bloco girado ou BER acima de 0,1 (Escolha) |
| `scripts/medidas_constelacao.py` | medidas diretas nas constelações (janela de fase, imagem de quadratura) | [M], diário, etapa A |
| `scripts/ber_vpi_simcopy.py` | figura de BER de 80 a 140 km | [V-dados] |
| `scripts/experimentos/busca_fase_vpi.py` | variantes da busca de fase das hipóteses 1, 3 e 4 (grade de 256QAM, dois estágios, partição dos anéis, frequência por blocos, desdobramento em 45°, fase conjunta, equalizador fixo) | [Pf09], [VV83], [Fa10], [Zh11], [MM10], [V] `Format = ReadFromFile` |
| `scripts/experimentos/equalizador_gradiente.py`, `equalizador_vpi.py` | equalizador por gradiente, só pelo raio ou com LMS (hipótese 2) | [Go80], [Fa09], [V] `AlgTDE = MMA` |
| `scripts/experimentos/fronteiras_vpi.py` | detector de estados de fase em passos de 45°, teste das fronteiras de bloco | [M] |
| `scripts/experimentos/ruido_radial_tangencial.py` | ruído tangencial contra radial, VPI e simcopy | [M] |
| `scripts/experimentos/compara_hipotese.py`, `figuras_ber_hipotese4.py` | placar das hipóteses contra o VPI; figuras de BER por distância | [M] |

## 3. Decisões, por tema

As decisões numeradas de D-1 a D-28 estão em `docs/decisoes.md`. As respostas completas, com o
raciocínio, estão em `docs/duvidas_abertas.md` (número da dúvida entre parênteses).

### 3.1 Processo

| Decisão | Base |
|---|---|
| Fontes por prioridade: Agrawal, código do professor, documentação. Separar fato de inferência | [G] |
| O Claude implementa. Testes e mutações são do grupo, por issues. Nunca editar um teste para ele passar | [G] |
| GPL nunca é importado (BIFROST, OptiCommPy, QAMpy, pyofss), só lido (D-19) | licenças |
| Dados do VPI, o PDF do artigo e o `CLAUDE.md` ficam fora do repositório (dúvida 38) | [G] |
| Decisão com várias partes vira uma proposta curta em `docs/proposta_*.md` | [G] |
| Todo resultado mostrado vem com figura | [G] |
| Fechamento do projeto (2026-11-19): a replicação do VPI está encerrada ("não precisa simular o VPI"). O DSP ganha MIMO, ortogonalização de Gram-Schmidt e compensação de PMD. O Gabriel escreve os testes, a documentação e o resto das implementações | professor e [G] |

### 3.2 Física e dispositivos

| Decisão | Base |
|---|---|
| Convenção de frequência vermelha do Agrawal, portadora exp(−iω₀t) (1) | [A] Eq. 2.4.3 e 2.4.5, [G] |
| O mux desloca cada canal e soma; o transmissor fica em banda base e grava `f0` (2) | [P], [G] |
| Potência do amplificador por canal (3); `PowerdBm` do VPI é a soma dos 4 canais (61, D-22) | [G]; [M] seção 14.16 da arquitetura |
| Laser direto: desconta o desvio de frequência do estacionário na polarização (4, a confirmar) | [G] |
| Laser de equações de taxa com `bias_ratio` e `modulation_ratio`; IMDD padrão com modulação externa (5) | [G] |
| Ruído de Langevin mantido com poucos fótons (7); fase inicial pela semente do canal (8) | [G] |
| Razão de extinção pelo desbalanço dos braços (12); S21 Bessel de ordem 3, banda infinita por padrão (13); pré-distorção arco-seno (14); skew como fase linear (15); divisão meio a meio em X e Y (16); renormalização da potência do modulador (17, alvo a confirmar); desbalanço de ganho com padrão 0 (18) | [A], [P], [V], [G] |
| ASE: n_sp pela Eq. 7.2.15 invertida (19); nas duas polarizações (20); campo de uma polarização vira dois com Y = 0 (21, a confirmar); modos do amplificador pela configuração (22); hν na frequência central (23); ASE na banda inteira, medida da entrada só na banda dos canais (24) | [A], [G] |
| Densidades bilaterais (25); ruído térmico de 5,76 pA/√Hz pelo Agrawal (26); Bessel de ordem 4 nos dois receptores (27); fotocorrente num `Waveform` elétrico (28) | [A] Eq. 4.4.9, [G] |
| Campo ao longo de z (32); aviso quando a PMD é ignorada (33); passo máximo de 1/10 da fibra com PMD (34); β₃ entra (35, padrão em aberto: 68) | [G] |
| FWM aceito como diferença do VPI (36, D-21); inclinação de dispersão do bloco de compensação do VPI era erro de digitação (37) | [G] |
| Sem limite G ≥ 0 no ganho até decidir (9: "sim" do professor contra "sempre o Agrawal" da 22) | em aberto |

### 3.3 Sistema IMDD

| Decisão | Base |
|---|---|
| Bits aleatórios pela semente do canal, não PRBS (49) | [G] |
| Filtro de Bessel no acionamento, não Nyquist (50) | [Aula], [G] |
| OOK é o padrão, PAM4 é bônus (51) | [G] |
| Span flexível, amplificador opcional (52) | [G] |
| Mux com soma ideal por padrão, filtro e perda opcionais (29); demux supergaussiano (30); saída na posição do canal (31) | [P], [V], [G] |
| Aviso de Brillouin acima de 5 mW por canal (42) | [A] seções 2.6.1 e 6.3 |

### 3.4 Sistema coerente

| Decisão | Base |
|---|---|
| `vbias` é parâmetro; no coerente a regra exige o nulo (10) | [A] seção 10.1.2, [G] |
| Excursão de acionamento 0,25 por padrão no IMDD (11): a de maior extinção. O Gabriel respondeu "não sei", e o valor foi mantido | [A] Eq. 10.1.5 |
| BPSK e uma polarização no coerente (53) | [Aula], [G] |
| Receptor ideal ou real conforme os parâmetros (54, `docs/proposta_ideal_real.md`) | [G] |
| I e Q binários nas linhas E11 a E13 da aula (55) | [G] |
| `drive_depth` pela excursão do Agrawal (56) | [A] seção 10.1.2 |
| Ruído térmico em cada um dos quatro fotodiodos (57) | [G] |
| BPSK pelo modulador IQ; laser e modulador como blocos separados (58) | [G] |
| Gray por eixo no 64QAM; 32QAM fora (69) | [G] |
| ADC em amostras por símbolo (40); potências em dBm no YAML (41) | [G] |
| DSP com símbolos conhecidos primeiro, cego depois, como marco próprio | roteiro antigo, tarefa 6 |

### 3.5 Replicação do VPI (encerrada em 2026-11-19)

| Decisão ou achado | Base |
|---|---|
| BER é a métrica principal; SER junto | [G] |
| As tabelas de BER do grupo contam I e Q separados (60); a constelação "ref" só serve para saber os símbolos (62); sem filtro casado no cenário (63) | [G], [M] |
| Desbalanço de fase de −5° e `drive_depth` 0,51 no YAML do 400ZR (64) | [M] `docs/comparacao_vpi.md`, seção 10 |
| Skew de uma amostra (1,044 ps) no lugar de 0,75 ps | [M] diário, etapa A. Inferência: o VPI arredonda o atraso para amostras inteiras (a confirmar com o grupo) |
| Receptor cego com janela de 43 no YAML do 400ZR | [M] diário, etapa A |
| No meio da faixa a BER do simcopy fica de 0,85 a 1,07 vez a do VPI | [M] `docs/comparacao_vpi.md`, seção 13 |
| Bordas: o DSP do VPI perde a fase caminhando em passos de 45°, com a fase conjunta nas duas polarizações. Perto do limiar a perda é rara e a taxa depende só da SNR. Hipóteses 2 e 3 descartadas; a 4 reproduz três de quatro assinaturas | [M] diário, hipóteses 1 a 4; notebook 09 |
| A "fração girada" das hipóteses 1 a 3 somava as polarizações: o modo de falha da hipótese 1 é o salto de 90° (corrigido no notebook 09) | [M] |

### 3.6 Figuras

| Decisão | Base |
|---|---|
| Estilo do guia do TCC, em português, tamanho dos notebooks por enquanto (`docs/proposta_figuras.md` para os outros destinos) | [GF], [G] |
| Figura gravada em PDF; no notebook, SVG | [G] |
| Constelação recebida como mapa de densidade quando ajuda | [G] |
| Comparação com o VPI: um painel por comprimento, VPI como curva, simcopy contínuo, sem marcar a fase perdida | [G] |
| O guia é ponto de partida: julgar se a figura fica bonita e legível | [G] |
| Estilo novo só do notebook 07 em diante; os notebooks 01 a 06 voltaram à versão anterior até o conserto (legendas acima das figuras, SVG largo) | [G] |

## 4. O que é inferência e falta confirmar

| Item | Onde |
|---|---|
| Edição do Agrawal usada na numeração das equações | seção 1.1 |
| Dúvidas respondidas com "considere" ou "sim" que interpretei: 4, 9, 17 e 21 | `docs/duvidas_abertas.md` |
| O VPI arredonda o skew para uma amostra | diário, etapa A |
| O VPI estima a fase em conjunto e sem distinguir 0° de 45° (o modelo da hipótese 4 reproduz, mas não prova) | diário, hipótese 4 |
| A generalização do fator Q para PAM e a BER aproximada do QAM quadrado não estão no Agrawal | `kernels/decision.py`, `kernels/demodulation.py` |
| O mapeamento do 64QAM estende a regra do 16QAM | `kernels/mapping.py`, dúvida 69 |

## 5. Onde está cada coisa

| Documento | Conteúdo |
|---|---|
| `docs/decisoes.md` | decisões fechadas e abertas (fonte única) |
| `docs/duvidas_abertas.md` | as 71 dúvidas, com resposta e raciocínio |
| `docs/proximos_passos.md` | o que falta, com as issues |
| `docs/arquitetura.md`, `docs/mapa_arquitetura.md` | arquitetura do pacote e lacunas |
| `docs/laser.md`, `docs/moduladores.md`, `docs/amplificadores.md`, `docs/fotodetectores.md`, `docs/mux_demux.md` | guias de implementação de cada dispositivo, com as equações |
| `docs/sistema_imdd.md`, `docs/sistema_coerente.md` | guias dos dois sistemas |
| `docs/comparacao_vpi.md`, `docs/replicacao_vpi_diario.md`, `docs/plano_replicar_vpi.md` | comparação com o VPI e o diário das hipóteses |
| `docs/comparacao_simuladores.md`, `docs/plano_validacao_cruzada.md` | outros simuladores |
| `docs/proposta_*.md` | propostas em discussão (figuras, ideal e real, potência de lançamento, recuperação de portadora e de relógio, API, bloco da fibra) |
| `docs/figuras.md` | figuras de validação (analítico contra implementação) |
| `docs/historico.md`, `docs/revisao_2026-11-14.md` | bugs encontrados e revisão do código |
| `docs/roadmap.md`, `docs/planejamento.md` | planos antigos, substituídos por `docs/proximos_passos.md` |
| notebooks 01 a 10 | um por dispositivo ou sistema; o 09 tem o 16QAM e o VPI, o 10 os outros formatos |
