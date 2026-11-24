# Sistema coerente digital: guia (2026-11-14)

Primeira rodada do coerente: os passos da aula, sem fibra, em BPSK e QPSK de uma polarização. O
código só acrescenta ao que existia: os blocos do IMDD ficaram intactos e os três exemplos IMDD dão
os mesmos números de antes (seção 8). O notebook é o `notebooks/08_sistema_coerente.ipynb`.

Fontes, na ordem de prioridade: Agrawal (*Sistemas de Comunicação por Fibra Óptica*), código do
professor (`draft.py`), documentação do projeto. Cada afirmação diz de onde vem.

---

## 1. O que a aula pede

**[Fato: `docs/04-class/context.md`]**

| Passo | Texto da aula | No notebook 08 |
|---|---|---|
| C1 | "Coerente digital. Inicialmente BPSK. Laser, MZM (modulação em amplitude, phase = 0)" | seção "BPSK: laser e MZM" |
| C2 | "Laser, MZM (modulação QAM, I(t) e Q(t)), se I e Q são OOK, é constante, não há recuperação de fase), photodetector" | seção "Modulador IQ e fotodetector direto" |
| C3 | "... mas no final do MZM tem um laser, photodetector" | seção "Um laser no fim" |
| C4 | "... mas no final do MZM tem um laser, balanceador, photodetector (aí vai conseguir recuperar info)" | seção "Híbrida de 90° e detecção balanceada" |

**[Dedução]** "Se I e Q são OOK, é constante" só vale com dois níveis simétricos no campo (±1):
com 0 e 1, $I^2 + Q^2$ assume 0, 1 e 2. Com o MZM no nulo, um sinal elétrico de dois níveis dá
campo ±1, e é isso que o notebook mostra (potência constante no centro dos símbolos, com desvio
relativo de 8e-16). A dúvida 55 continua aberta para confirmação.

---

## 2. Transmissor

### 2.1 Bits e símbolos (`kernels/mapping.py`)

**[Fato: `draft.py`, `QAM_dem` com M = 4]** Níveis normalizados pelo maior nível (±1 por
quadratura). No QPSK o primeiro bit do símbolo escolhe a quadratura Q (0 → +1, 1 → −1) e o segundo
a fase I, com a mesma regra. É um código de Gray. O BPSK usa a mesma regra num eixo: 0 → +1,
1 → −1.

Só BPSK e QPSK por enquanto. Nos dois, cada bit é uma decisão binária numa quadratura, o que
permite o fator Q do Agrawal (seção 4).

### 2.2 Pulso de Nyquist (`kernels/pulse_shaping.py`, bloco `NyquistDriver`)

**[Decidido: Gabriel, 2026-11-14]** Cosseno levantado, com o `rolloff` que já estava na
configuração (0,18).

**[Fato: `draft.py`, `DAC_Nyquist`]** O professor monta um trem de impulsos e multiplica o espectro
por uma janela plana até $(1-\rho)R_s/2$ que cai em cosseno até $(1+\rho)R_s/2$. É o espectro de
cosseno levantado inteiro no transmissor (não a raiz), aplicado na FFT. O kernel faz a mesma conta.

Consequências, medidas no notebook:
- no centro de cada símbolo o sinal vale o símbolo (diferença de 2e-15);
- entre os símbolos I(t) chega a 2,1 vezes o nível do símbolo;
- a potência média do sinal é $(1 - \rho/4)$ da potência no centro dos símbolos (0,956 medido,
  0,955 pela fórmula).

O filtro é periódico: a janela da simulação tem que conter um número inteiro de símbolos, e o
bloco recusa se não contiver.

### 2.3 Laser e modulador IQ (`components/tx/iq_modulator.py`)

**[Fato: Agrawal, seção 10.1.2 e Fig. 10.3]** Um MZM em push-pull em cada braço de um
interferômetro externo, com defasagem de π/2 entre os braços. O kernel `iq_field` já existia.

O bloco `IQModulator` recebe duas entradas, a portadora e o sinal elétrico. O laser CW é um estágio
próprio (`ch0.tx.laser`), porque o oscilador local homódino usa a mesma saída.

- **Polarização no nulo.** **[Fato: `draft.py`, `DPQAM_transmitter`: `Vpi = Vbias = 2.5`]** Cada
  MZM fica em $V_b = V_\pi$, onde o campo é $\sin(\pi v/V_\pi)$. **[Provisório: dúvida 10, opção
  (b)]** O nulo sai do formato e `mzm.vbias` não é usado no coerente. **Mudou em 2026-11-16:** o
  Gabriel escolheu a opção (a). O bloco continua supondo o nulo, e a validação exige
  `mzm.vbias = mzm.vpi` nos canais coerentes.
- **Pré-distorção.** **[Decidido: M6, 2026-10-22]** $v = (V_\pi/\pi)\arcsin(m\,s)$, e o campo fica
  $m\,s$. Onde $|m\,s| > 1$ (picos entre os símbolos) o sinal é ceifado.
- **`drive_depth` ($m$).** Excursão do campo nos níveis externos da constelação. **[Sem fonte:
  dúvida 56]** O padrão é 1 (o MZM inteiro). Os exemplos usam 0,4, para o pico de 2,1 ficar em 0,85
  e nada ser ceifado.
- **S21, skew, desbalanços de fase e de ganho, extinção finita.** As decisões M4, M5, M7 e a do
  desbalanço de ganho, de 2026-10-22, aplicadas como no IMDD.
- **Potência de saída.** Fica com o `PowerControl`, quando o canal tem `launch_power`, como no
  IMDD (decisão M9, com o alvo ainda a confirmar).

O BPSK passa pelo mesmo modulador IQ, com o braço Q parado no nulo. Com extinção infinita isso
equivale a um MZM só (o passo C1 da aula, que o notebook monta à mão com `mzm_field`), com 6 dB a
menos de potência, que o `PowerControl` repõe.

---

## 3. Receptor

### 3.1 Oscilador local

**[Fato: Agrawal, Eq. 4.5.1 a 4.5.5]** Num fotodetector só,
$P(t) = P_s + P_{LO} + 2\sqrt{P_s P_{LO}}\cos(\omega_{IF}t + \phi_s - \phi_{LO})$.

Duas fontes, escolhidas por `rx.lo_source` (dúvida 54):

| `lo_source` | Bloco | O que é |
|---|---|---|
| `tx_laser` | `LockedLocalOscillator` | Cópia do laser do transmissor, com a potência `rx.lo.power`. O ruído de fase é o mesmo do sinal e se cancela: homódino ideal (Agrawal, seção 4.5.2). Com `lo_frequency_offset` diferente de zero, sobra um desvio de frequência puro. |
| `free` (padrão) | `CWLaser` com o papel de semente `lo_laser_phase` | Outro laser, com largura de linha e desvio de frequência próprios. |

**[Limite desta rodada]** Não há recuperação de portadora. Com `free`, ou com desvio de frequência,
a constelação gira e a BER medida não representa o enlace. A regra de validação
`CoherentNeedsCarrierRecovery` avisa (aviso, não erro: o sinal continua válido para ver o
batimento).

### 3.2 Híbrida de 90° e detecção balanceada (`kernels/coherent.py`, bloco `CoherentFrontEnd`)

**[Fato: `draft.py`, `HybridNetwork`]**
$E_{1,2} = \tfrac{1}{2}(E_s \pm E_{LO})$ e $E_{3,4} = \tfrac{1}{2}(E_s \pm i\,E_{LO})$.

**[Fato: Agrawal, Eq. 10.2.4]** No receptor balanceado a diferença das correntes elimina o termo
constante. Com as quatro saídas, $I + iQ = R\,E_s E_{LO}^*$ (Agrawal, seção 8.7.3 e Fig. 8.39: os
quatro fotodiodos recuperam as partes real e imaginária de $E_x E_{lo}^*$).

**[Fato: `draft.py`, `DPQAM_receiver`]** Cada um dos quatro fotodiodos tem o seu ruído de disparo e
térmico. A diferença passa pelo passa-baixas do receptor (`photodiode.bandwidth_factor`, Bessel de
ordem 4).

A relação vale para envelopes, então não depende da convenção da portadora. Um oscilador local em
$+\Delta f$ (azul) faz a constelação girar como $\exp(+i2\pi\Delta f\,t)$.

Receptor de uma polarização: mistura a componente X do sinal com a componente X do oscilador
local. A potência do sinal em Y não é detectada.

---

## 4. A medida (`kernels/demodulation.py`)

A medida amostra a saída analógica do receptor. Não há ADC nem DSP.

1. **Instante dos símbolos.** **[Decidido: Gabriel, 2026-11-14]** Correlação circular entre a
   corrente e o trem de símbolos transmitido, procurada na janela inteira. O simulador conhece os
   símbolos, então não há sequência de sincronismo.
2. **Ganho complexo.** $g = \sum r\,a^*/\sum|a|^2$ tira a escala da corrente e uma rotação
   constante. Não corrige fase que varia no tempo.
3. **Decisão e contagem.** Ponto mais próximo, e comparação com os bits transmitidos.
4. **EVM.** $\sqrt{\sum|r/g - a|^2/\sum|a|^2}$. Com ruído aditivo, SNR = 1/EVM².
5. **Fator Q.** **[Fato: Agrawal, Eq. 10.3.2 e 10.3.7]** No BPSK, $Q = (I_1 - I_0)/(\sigma_1 +
   \sigma_0)$ com $I_0 = -I_1$, e $\mathrm{BER} = \tfrac{1}{2}\mathrm{erfc}(Q/\sqrt{2})$. No QPSK
   cada quadratura é uma decisão do mesmo tipo, e a medida junta as decisões de todos os bits.

Os símbolos das bordas da janela ficam fora da medida, como no IMDD (`edge_guard_symbols`).

`CoherentMeasurement` tem `q`, `ber_q`, `ber_counted`, `errors` e `n_bits`, com os mesmos nomes de
`QMeasurement`, mais `evm`, `snr_db`, `gain`, `phase` e `delay`. Por isso `sweep`, `plot_sweep` e a
CLI servem aos dois tipos de canal.

### 4.1 Alternativa documentada: sequência de sincronismo

**[Pedido do Gabriel, 2026-11-14: fica documentada, não implementada]**

**[Fato: `draft.py`, `QAM_transmitter` e `QAM_receiver`]** O professor monta o quadro como zeros,
`N_sync` símbolos BPSK (±1, de uma semente própria), a informação e zeros. No receptor, gera a
mesma sequência, passa pelo `DAC_Nyquist`, correlaciona com o sinal recebido (`np.correlate`) e
toma o máximo como início do quadro. Daí saem os índices da sequência e da informação, amostrados
uma vez por símbolo.

O que mudaria no simcopy para adotá-la:
- a fonte de bits passaria a montar o quadro (zeros, sincronismo, informação), e a geometria do
  sinal teria que distinguir símbolos de sincronismo de símbolos de informação;
- a medida correlacionaria só com a sequência de sincronismo, e contaria erros só na informação;
- a sequência de sincronismo também serviria para estimar o ganho complexo sem conhecer os dados.

É mais realista (um receptor de verdade não conhece os dados) e é o que o DSP da próxima rodada
precisa se a recuperação de portadora usar símbolos conhecidos. A correlação com os símbolos
transmitidos dá o mesmo instante com menos peças, e resolve também o ponto 1 de
`docs/revisao_2026-11-14.md` (instante de decisão do IMDD em WDM).

---

## 5. Como rodar

```
simcopy run examples/coerente_qpsk.yaml -o resultados/qpsk
simcopy sweep examples/coerente_qpsk.yaml --param "channels[0].launch_power" --range 50e-9 250e-9 50e-9
```

```python
from simcopy.pipelines.runner import run_file
from simcopy.plotting import plot_constellation
r = run_file("examples/coerente_qpsk.yaml")
print(r.summary())
plot_constellation(r)
```

Estágios de um canal (`r.ids`):

```
ch0.tx.bits        BitSource
ch0.tx.driver      NyquistDriver             <- bits
ch0.tx.laser       CWLaser                   <- driver (só para a grade e o canal)
ch0.tx.modulator   IQModulator               <- laser, driver
ch0.tx.power       PowerControl              <- modulator   (se houver launch_power)
link.*             fibra, amplificadores     (lista vazia nos exemplos)
ch0.rx.lo          LockedLocalOscillator     <- ch0.tx.laser   (ou CWLaser <- sinal, se livre)
ch0.rx.frontend    CoherentFrontEnd          <- sinal, lo
```

Exemplos: `coerente_bpsk.yaml` (−44 dBm no receptor) e `coerente_qpsk.yaml` (−41 dBm), os dois com
BER perto de 1e-3, onde dá para contar erros com 2^15 símbolos.

Campos novos na configuração, todos com padrão: `tx.drive_depth` (1.0), `tx.predistortion` (true),
`rx.lo_source` ("free"). Formato novo: `BPSK`. Sem fibra: `link: {spans: []}`.

---

## 6. O que a montagem ainda recusa ou não faz

Recusa, com `ConfigError`: formatos além de BPSK e QPSK, dupla polarização, mais de um canal
coerente, canais IMDD e coerentes no mesmo sistema.

Não faz:
- **Recuperação de portadora** (frequência e fase). O professor usa compensação grosseira de
  frequência e busca cega de fase em blocos de 16 símbolos, com 21 ângulos (`draft.py`,
  `QAM_receiver`).
- **ADC.** `rx.adc` (2 amostras por símbolo, 12 bits) não é usado: falta o fundo de escala do
  quantizador e a recuperação de relógio (dúvida 40).
- **Compensação de dispersão no DSP.** Uma fibra em `link.spans` roda, mas o sinal chega sem
  compensação.
- **Diversidade de polarização.** A PMD vem ligada por padrão nas fibras e gira a polarização. O
  receptor só vê X: com fibra, desligue a PMD (`birefringence: {enabled: false}`).

---

## 7. Decisões

**De 2026-11-14 (Gabriel):**

| Pergunta | Resposta |
|---|---|
| Escopo da rodada | Passos da aula, BPSK e QPSK, sem fibra, com BER contada e EVM |
| Pulso | Nyquist (cosseno levantado, rolloff 0,18) |
| Código do IMDD | Só acrescentar. No executor, só o despacho por tipo de canal |
| Sincronismo | Correlação com os símbolos transmitidos. A sequência de sincronismo fica documentada (seção 4.1) |

**Tomadas na implementação, a confirmar** (registradas em `docs/duvidas_abertas.md`, seção 11):

1. **Dúvida 10:** MZM no nulo definido pelo formato. `mzm.vbias` não é usado no coerente.
   Respondida em 2026-11-16 (opção (a)): a validação passou a exigir `mzm.vbias = mzm.vpi`.
2. **Dúvida 53:** BPSK entrou nos formatos. Só uma polarização por enquanto.
3. **Dúvida 54:** as duas fontes de oscilador local existem. Híbrida e detecção balanceada ideais.
4. **Dúvida 55:** leitura dos passos C2 a C4 da seção 1.
5. **Dúvida 56 (nova):** `drive_depth`, padrão 1,0 e 0,4 nos exemplos.
6. **Dúvida 57 (nova):** ruído térmico em cada um dos quatro fotodiodos, como o professor, e não um
   por par balanceado.
7. **Dúvida 58 (nova):** o BPSK usa o modulador IQ com o braço Q no nulo, e não um bloco de um MZM.
8. **Dúvida 59 (nova):** `insertion_loss_db` é aplicada em cada MZM interno, como no IMDD.

---

## 8. Verificações feitas (2026-11-14)

- **IMDD intacto.** `imdd_ook_mzm`: Q = 26,608. `imdd_ook_dml`: Q = 1,511, 42 erros.
  `imdd_wdm_4ch`: Q = 34,919 / 32,770 / 33,912 / 34,663. Iguais aos de antes da mudança.
- **Testes existentes:** 73 passam (sem o rascunho `test_birefringence.py`).
- **Camadas:** contrato do import-linter mantido.
- **Pulso de Nyquist:** sinal igual ao símbolo no centro dos símbolos, com diferença de 2e-15.
- **BPSK:** campo dividido pela portadora real (parte imaginária de 1e-16), com ±0,4 no centro dos
  símbolos.
- **QPSK com fotodetector direto:** potência no centro dos símbolos constante, com desvio relativo
  de 8e-16.
- **Um fotodetector com oscilador local:** corrente igual a $R(P_s + P_{LO})/4$ mais o batimento
  (diferença de 4e-19 A). Correlação com I de 1,000 e com Q de 0,005.
- **Detecção balanceada:** $(I_1 - I_2) + i(I_3 - I_4) = R\,E_s E_{LO}^*$ (diferença de 4e-19 A).
  Sem ruído, 0 erros em 8148 bits e EVM de 3,6 %, que vem do filtro do receptor.
- **BER contra o fator Q e contra a conta do ruído** (Agrawal, Eq. 4.5.9 e 4.5.10, somando os
  quatro fotodiodos; banda de ruído de 12,03 GHz; pico do pulso depois do filtro de 0,979):

| Formato | P [dBm] | Contada | Por Q | Conta do ruído |
|---|---|---|---|---|
| BPSK | −46 | 8,09e-3 (265 erros) | 7,85e-3 | 7,94e-3 |
| BPSK | −44 | 1,22e-3 (40) | 1,21e-3 | 1,20e-3 |
| QPSK | −43 | 7,83e-3 (513) | 8,04e-3 | 8,00e-3 |
| QPSK | −41 | 1,45e-3 (95) | 1,25e-3 | 1,21e-3 |
| QPSK | −39 | 1,37e-4 (9) | 7,51e-5 | 6,73e-5 |

  Abaixo de 1e-4 há poucos erros para contar. Nos pontos de BER baixa do QPSK a contagem fica acima
  da estimativa (95 contra 82 esperados, 9 contra 5). Pode ser estatística, ou a interferência
  entre símbolos do filtro, que o fator Q trata como gaussiana. Não foi investigado.
- **Oscilador local livre a 500 MHz:** a constelação vira um anel e a BER vai a 0,49, com o aviso
  da validação.

---

## 9. Insumo para as issues de teste (com o grupo)

- `kernels/mapping.py`: ida e volta `demap(map_bits(b)) == b`, Gray (vizinhos diferem em um bit),
  níveis em ±1.
- `kernels/pulse_shaping.py`: critério de Nyquist (sinal igual ao símbolo no centro), $H(0) = 1$,
  banda nula acima de $(1+\rho)R_s/2$, potência média $(1-\rho/4)$, `rolloff = 0` e `1`.
- `kernels/coherent.py`: soma das quatro potências igual a $P_s + P_{LO}$, e as duas identidades
  das diferenças.
- `kernels/demodulation.py`: atraso conhecido recuperado (inclusive negativo e maior que 2
  símbolos), ganho e fase conhecidos recuperados, EVM contra ruído de variância conhecida, BER
  contra a fórmula.
- `IQModulator`: campo $m\,s$ com pré-distorção, ceifamento, skew, desbalanços, extinção finita.
- `CoherentFrontEnd`: variância do ruído contra $4 i_T^2\Delta f + 2qR(P_{LO} + P_s)\Delta f$.
- Montagem: recusa de formatos, dupla polarização e mais de um canal. Regressão dos números da
  seção 8.

---

## 10. Próximos passos

1. Recuperação de portadora (frequência e fase), para o oscilador local livre.
2. 16QAM (mapeamento, decisão e medida por nível).
3. Fibra com compensação de dispersão no receptor (Agrawal, seção 8.7.3).
4. ADC e reamostragem para a taxa do DSP.
5. Dupla polarização.

---

## 11. Segunda rodada (2026-11-14): 16QAM, dupla polarização, fibra, vários canais e DSP

Pedido do Gabriel: o resto do coerente, com prioridade para o DP-16QAM, que tem parâmetros e dados
do VPI (`reference/vpi/parametros_vtmu.md`). A comparação com o VPI está em
`docs/comparacao_vpi.md`. O notebook é o `notebooks/09_coerente_dp16qam.ipynb`.

A primeira rodada continua como estava: sem `rx.dsp`, os exemplos `coerente_bpsk.yaml` e
`coerente_qpsk.yaml` dão os mesmos números.

### 11.1 Transmissor

- **16QAM** (`kernels/mapping.py`). **[Fato: `draft.py`, `QAM_dem` com M = 16]** Os dois primeiros
  bits escolhem o nível em fase e os dois últimos o nível em quadratura, com o mesmo código de Gray:
  00 → −3, 01 → −1, 11 → +1, 10 → +3. Níveis divididos por 3. O VPI tem `GrayMapping = No`: a taxa
  de erro de símbolo e a SNR não dependem disso, a BER depende.
- **Raiz de cosseno levantado** (`tx.pulse: root_raised_cosine`). **[Fato: VPI, `NyquistResponse =
  squareRootRaisedCosine`, `RollOff = 0.18`]** O padrão continua o cosseno levantado inteiro do
  professor.
- **Dupla polarização.** **[Decidido: M8, 2026-10-22]** O laser é dividido meio a meio e cada metade
  passa pelo seu modulador IQ. O sinal elétrico tem uma linha por polarização.
- **Vários canais.** O mux do IMDD serve como está. O do cenário do VPI é o filtro gaussiano de
  ordem 3, 70 GHz, com 3 dB de perda.

### 11.2 Receptor

- **Filtro do canal.** Com mais de um canal, um filtro supergaussiano de banda
  `rx.optical_bandwidth` e ordem `rx.optical_filter_order` (VPI, `FilterOpt`: 70 GHz, ordem 3).
- **Oscilador local** na frequência do canal mais `lo_frequency_offset`.
- **Diversidade de polarização.** **[Fato: Agrawal, seção 8.7.3 e Fig. 8.39]** O sinal é separado
  em X e Y, o oscilador local é dividido meio a meio (professor: `X_fraction = 0.5`) e cada
  polarização tem a sua híbrida e os seus quatro fotodiodos.

### 11.3 DSP (`rx.dsp`, `kernels/dsp.py`)

**[Decisão registrada no pacote, `docs/decisoes.md`]** O receptor conhece o comprimento e a
dispersão do enlace (D-7) e, nesta fase, usa os símbolos transmitidos, deixando o DSP cego do VPI
como marco próprio.

Blocos, na grade da simulação:

| Estágio | Bloco | O que faz |
|---|---|---|
| `rx.cd` | `DispersionCompensator` | Multiplica o espectro por $\exp[-\tfrac{i}{2}\beta_2 L (2\pi f_g)^2]$, com $f_g$ a frequência na grade óptica (Agrawal, seção 8.7.3). Tira também o atraso do canal em relação ao centro da banda |
| `rx.frequency` | `FrequencyOffsetCompensator` | Opcional (`frequency_offset: true`), desligado por padrão: estima o desvio de frequência pela raia da quarta potência e o desfaz |
| `rx.filter` | `MatchedFilter` | Raiz de cosseno levantado (só com `tx.pulse: root_raised_cosine` e `matched_filter: true`) |

Na medida, com os símbolos conhecidos:

1. instante dos símbolos e desvio de frequência, pelo produto diferencial (abaixo)
2. fase em blocos de `block_length` símbolos e a matriz de polarização 2×2, constante na janela
3. fase fina, pela média móvel de `phase_window` símbolos
4. equalizador de mínimos quadrados com `equalizer_taps` coeficientes a 2 amostras por símbolo
5. ganho complexo por polarização, decisão e contagem

O equalizador é linear nos complexos, como o do VPI: não corrige o desbalanço de quadratura nem o
vazamento de portadora do transmissor, que ficam na constelação.

**Instante e desvio de frequência.** O produto de cada amostra pelo conjugado da amostra um
símbolo antes, somado nas polarizações, não depende da matriz de polarização nem da fase dos
lasers, e transforma um desvio de frequência numa rotação constante. A correlação dele com o mesmo
produto dos símbolos transmitidos tem o pico no instante dos símbolos, e a fase do pico dá o desvio
(vale até metade da taxa de símbolos).

**[Medido em 2026-11-15]** A raia da quarta potência não serve de base para o 16QAM neste cenário:
a altura dela variou de 7 a 84 vezes a mediana do espectro conforme a semente, e com $2^{13}$
símbolos a estimativa errou (−1222 MHz no lugar de 500 MHz) e a medida inteira falhou. Filtrar o
canal antes e suavizar o espectro não resolveram. O produto diferencial acertou em todos os casos
testados ($2^{12}$ a $2^{14}$ símbolos, duas sementes), com a SNR entre 16,10 e 16,18 dB.

Para o 16QAM, `ber_q` é a BER que a SNR medida daria com ruído gaussiano e Gray,
$\tfrac{3}{8}\mathrm{erfc}(\sqrt{\mathrm{SNR}/10})$ (fórmula usual, a mesma do ajuste do pacote, não
está no Agrawal), e `q` é o fator Q equivalente.

### 11.4 O que não é modelado

- **ADC** (2 amostras por símbolo, 12 bits): o `Range = 1.0` do VPI não diz em que escala a corrente
  chega ao conversor.
- **Inclinação da dispersão** (`DispersionSlope` do VPI): o kernel da fibra só tem $\beta_2$.
- **Sequência PRBS de ordem 7** do VPI: os bits são aleatórios pela semente do canal (decisão de
  2026-11-12).
- **`FWM = No` do VPI**: o SSFM do campo agregado inclui a mistura de quatro ondas.
- **DSP cego.**

### 11.5 Exemplos

- `examples/coerente_dp16qam.yaml`: um canal DP-16QAM sem fibra, com os defeitos de transmissor do
  VPI.
- `examples/400zr_4ch_75ghz.yaml`: o cenário do VPI. Ganhou `pulse`, `drive_depth`, `mux`,
  `optical_filter_order` e `dsp`. Com os $2^{18}$ símbolos do arquivo, cada ponto leva muito tempo
  e memória: `scripts/compara_vpi.py` roda com menos símbolos.

---

## 12. Terceira rodada (2026-11-16): DSP cego, ADC, 64QAM, inclinação da dispersão, MI e GMI

Pedido do Gabriel. O PAM4, que é do IMDD, está em `docs/sistema_imdd.md`, seção 11. A comparação do
DSP cego com o VPI está em `docs/comparacao_vpi.md`, seção 9.

Nada do que existia mudou de resultado: sem `rx.dsp.mode: blind`, sem `dispersion_slope` e com os
formatos antigos, os exemplos dão os mesmos números.

### 12.1 DSP cego (`rx.dsp.mode: blind`)

Nenhuma etapa usa os símbolos transmitidos. O receptor ganha dois estágios, na taxa do DSP
(`DigitalSignal`, fora da grade da simulação):

| Estágio | Bloco | O que faz |
|---|---|---|
| `rx.adc` | `ADC` | Amostra a 2 amostras por símbolo e quantiza (12.2) |
| `rx.dsp` | `BlindReceiver` | A cadeia cega de `kernels/dsp.py`, `blind_receiver`, e entrega os símbolos |

A cadeia, com a fonte de cada etapa:

1. **Separação de polarização no espaço de Stokes.** **[Fato: documento do VPI, seção 8.3]** O VPI
   tem a etapa. Os vetores de Stokes das amostras ficam em torno de um plano, e a normal do plano
   dá os dois estados transmitidos (Szafraniec, Nebendahl e Marshall, Opt. Express 18, 17928,
   2010). No VPI ela vem depois do equalizador. Aqui vem antes, para o equalizador partir de
   polarizações já quase separadas. Pode ser desligada (`stokes_demux: false`).
2. **Equalizador adaptativo 2×2**, 7 coeficientes a 2 amostras por símbolo. **[Fato: Agrawal, seção
   8.7.3]** A matriz inversa do enlace pode ser construída "a partir do próprio sinal recebido,
   usando um algoritmo conhecido como algoritmo de módulo constante". **[Fato: VPI]** `AlgTDE =
   MMA`, `NTaps = 7`, coeficientes iniciais com 1 no centro. O algoritmo de módulo constante
   (Godard, IEEE Trans. Commun. 28, 1867, 1980) roda `cma_iterations` passagens pelo bloco. Se as
   duas saídas caírem na mesma polarização, a segunda é refeita ortogonal à primeira. A etapa
   dirigida por raio, que é a minha leitura do `MMA` com `RadiiSq` do VPI (Fatadin, Ives e Savory,
   J. Lightw. Technol. 27, 3042, 2009), existe e vem desligada (`rde_iterations: 0`).
3. **Desvio de frequência**, pela raia da quarta potência dos símbolos já equalizados. Antes do
   equalizador essa raia não servia (seção 11.3). Depois dele acertou em todos os casos testados
   (500,1 a 502,1 MHz para 500 MHz).
4. **Busca cega de fase.** **[Fato: VPI]** `BlindPhaseSearch`, `NbTestAngles = 64`, `CombinedPol =
   Yes`. **[Fato: `draft.py`, `QAM_receiver`]** O professor faz a busca em blocos de 16 símbolos
   com 21 ângulos. Aqui: 64 ângulos, janela deslizante de `phase_window` símbolos (Pfau, Hoffmann e
   Noé, J. Lightw. Technol. 27, 989, 2009), conjunta nas duas polarizações.
5. **Ajuste dirigido por decisão**, `dd_iterations` vezes: os coeficientes são refeitos por
   mínimos quadrados no bloco, tendo como alvo as próprias decisões, e a portadora é recuperada de
   novo. **[Fato: VPI]** O bloco lista `Mode = Decision-Directed`.

O que fica para a medida (`measure_blind`), porque só os dados resolvem: qual saída é qual
polarização, o atraso em símbolos e um giro constante de 90°. É o papel do `SynchDSP` do VPI. A
medida não ajusta ganho nem fase com os dados.

**O que os testes mostraram (16QAM, cenário de 4 canais):**

- **O critério do módulo não chega à solução de mínimo erro.** Só com o equalizador de módulo, a
  SNR fica em 15,4 dB onde o DSP com símbolos conhecidos dá 16,1 dB, mesmo com a fase tirada pelos
  dados. Com o ajuste por decisão, 16,1 dB.
- **A janela da busca cega de fase decide os saltos de 90°.** Com 21 símbolos e busca por
  polarização: 23 saltos numa polarização a 12,9 dB de SNR. Com 81 símbolos: 15 saltos numa
  polarização e nenhum na outra. Com 81 símbolos e busca conjunta: nenhum, de 11 a 16 dB. O ruído
  de fase destes lasers (500 kHz cada, a 59,84 GBd) é lento o bastante para a janela de 81.
- **A etapa dirigida por raio atrapalhou com pouca SNR.** A 13 dB, com passo de 1e-3, uma
  polarização ficou 6 dB abaixo da outra depois dela. Sem ela, e com 4 ajustes por decisão, as duas
  convergem.
- **Os outros formatos.** QPSK e BPSK de uma polarização, com oscilador local livre a 500 MHz e
  janela de 81: o modo cego fica a até 0,25 dB do modo com símbolos conhecidos (−38 e −41 dBm), sem
  saltos. DP-QPSK: 18,60 dB no modo cego e 18,71 dB com símbolos conhecidos. O QPSK só passou a
  funcionar no modo cego depois da correção do ponto 38 de `docs/revisao_2026-11-14.md`.

**Padrões** (`DSPConfig`): `step_size` 1e-3, `cma_iterations` 4, `rde_iterations` 0,
`dd_iterations` 4, `test_angles` 64, `stokes_demux` true. O passo de 5e-4 do VPI não é comparável,
porque a normalização do sinal de lá não está no documento.

O `phase_window` continua 21. **Correção (2026-11-17):** eu dizia que 21 era a janela do VPI, mas
`NPreSymbols = 21` são os símbolos de cada lado: a janela do VPI é 43, medida nas constelações
(`docs/replicacao_vpi_diario.md`, etapa A). Com 43 o receptor cego reproduz a não linearidade do
VPI. No modo cego a escolha é uma troca, medida no cenário de 4 canais (`docs/comparacao_vpi.md`,
seção 9, com 21 e 81):

| `phase_window` | Saltos de fase | Em alta potência (17 dBm, 125 km) |
|---|---|---|
| 21 | abaixo de cerca de 13 dB de SNR (2 de 9 realizações a 12,3 dB, nenhuma de 9 a 13,1 dB) | 15,17 dB |
| 81 | nenhum de 10 a 16 dB | 14,64 dB: não acompanha a fase não linear |

Com 21 o simcopy cego fica de 0,21 a 0,33 dB acima do VPI e reproduz a queda de baixa potência
dele, 1 dB abaixo. Para estudar SNR baixa sem saltos, use 81. Para replicar o VPI, use 43.

A medida conta os saltos de fase (`phase_slips`), comparando o giro de cada bloco de 128 símbolos
com os dados.

### 12.2 ADC (`kernels/adc.py`, bloco `ADC`)

**[Fato: VPI, `Rx_DualPol`]** `SamplingRate = 2*SymbolRate`, `TransferFunctionADC = Bipolar`,
`Resolution = 12 bit`, `Range = 1.0`.

Quantizador uniforme bipolar, na parte real e na imaginária de cada polarização. O `Range = 1.0` do
VPI não diz em que escala a corrente chega, então o fundo de escala é `adc.full_scale_rms` vezes o
valor rms de cada componente (padrão 4,0, sem fonte: dúvida 65). `resolution_bits: 0` só amostra.

Sem ceifamento, $\mathrm{SQNR} = 6{,}02\,b + 4{,}77 - 20\log_{10}(F/\sigma)$ dB. Medido com sinal
gaussiano: 19,0 dB com 4 bits (fórmula 19,3), 28,8 dB com 6 bits (28,9), 40,6 dB com 8 bits (40,9).
Com 12 bits e fundo de escala de 4σ o ceifamento limita em 52,9 dB, muito acima dos 16 dB do sinal.

O ADC só entra no modo cego. No modo com símbolos conhecidos a medida continua amostrando a grade
da simulação, para não mudar os números já comparados com o VPI. Não há recuperação de relógio: o
relógio do ADC é o da simulação, e a fase de amostragem fica com o equalizador fracionário.

### 12.3 64QAM

O `draft.py` não tem 64QAM. O mapeamento segue a regra do 16QAM do professor com três bits por
eixo: a primeira metade dos bits escolhe o nível em fase, a segunda o nível em quadratura, com o
código de Gray refletido, e os níveis são divididos por 7. O 16QAM passou a usar a mesma função e
continua batendo com a tabela do professor.

Para QAM quadrado, `ber_q` é a BER que a SNR medida daria com ruído gaussiano:
$\tfrac{4}{k}(1 - \tfrac{1}{\sqrt M})\,\tfrac12\mathrm{erfc}\sqrt{3\,\mathrm{SNR}/(2(M-1))}$. Conferida
com ruído gaussiano sintético no 64QAM: 2,42e-2 contada contra 2,42e-2 a 18 dB, 8,37e-3 contra
8,49e-3 a 20 dB, 1,77e-3 contra 1,75e-3 a 22 dB.

Exemplo: `examples/coerente_dp64qam.yaml`, com transmissor mais limpo que o do VPI (valores de
exemplo, sem fonte). A −20 dBm: 23,24 dB com os símbolos conhecidos e 23,10 dB cego.

O 32QAM (constelação em cruz) continua fora.

### 12.4 Inclinação da dispersão (`fiber.dispersion_slope`)

**[Fato: Agrawal, Eq. 2.3.13]** $S = (2\pi c/\lambda^2)^2\beta_3 + (4\pi c/\lambda^3)\beta_2$, então
$\beta_3 = (S + 2D/\lambda)\,[\lambda^2/(2\pi c)]^2$. **[Fato: Agrawal, Eq. 2.4.4]** $\beta(\omega)$
tem o termo $\tfrac{\beta_3}{6}(\Delta\omega)^3$. Na convenção do projeto o eixo da FFT é
$-\Delta\omega$, então o termo cúbico entra com o sinal trocado no operador linear da fibra.

- `dispersion_slope: null` (padrão): sem $\beta_3$, como antes.
- `dispersion_slope: 0.0`: $D$ constante no comprimento de onda, o que ainda dá $\beta_3 =
  2D/\lambda$ vezes o fator (0,033 ps³/km com $D$ = 16).
- `dispersion_slope: 0.08e3` (VPI, 0,08 ps/(nm²·km)): $\beta_3$ = 0,163 ps³/km.

A compensação de dispersão do DSP desfaz também $\beta_3$.

Verificações: alargamento rms de um pulso gaussiano só com $\beta_3$ contra a Eq. 2.4.22 do Agrawal
(1,079664, 1,910704 e 4,191352 em 5, 20 e 50 km, iguais nas seis casas). Atraso do centro do pulso
de +0,204 ps em 5 km, contra $\beta_3 L/(4T_0^2)$ = 0,2035 ps: o sinal está certo. A compensação
do DSP desfaz a fibra com erro de 2e-13.

### 12.5 MI e GMI (`kernels/information.py`)

Não estão no Agrawal. Seguem Alvarado, Agrell, Lavery, Maher e Bayvel, J. Lightw. Technol. 33, 4338
(2015): símbolos equiprováveis e um decodificador que supõe ruído gaussiano com a variância medida
nos próprios símbolos. A MI é a taxa alcançável com decodificação por símbolo, a GMI com
decodificação bit a bit, que é o caso dos códigos usados em sistemas ópticos. Saem em bits por
símbolo e por polarização, em `CoherentMeasurement.mi` e `.gmi`.

Verificações com ruído gaussiano sintético:

- MI contra uma quadratura de Gauss-Hermite em duas dimensões, feita à parte: diferença de até
  0,006 bit em QPSK, 16QAM e 64QAM, de 5 a 20 dB.
- GMI do 16QAM contra a soma da informação de cada bit, por integração numérica: diferença de até
  0,003 bit, de 2 a 16 dB.
- No QPSK a GMI é igual à MI, como tem de ser.

No 16QAM de Gray a GMI fica a menos de 0,001 bit da MI acima de 10 dB. Abaixo disso se afasta
(1,93 contra 1,97 a 5 dB).

Só para o coerente. O IMDD não tem MI nem GMI.

### 12.6 O que ainda não existe

- Recuperação de relógio (o relógio do receptor é o da simulação).
- Retropropagação digital, código corretor de erros, 32QAM.
- Testes automáticos de tudo o que entrou nas três rodadas.
