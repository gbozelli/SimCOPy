# Dúvidas abertas do simcopy (2026-10-22)

Todas as perguntas que estão abertas nos guias (`docs/laser.md`, `moduladores.md`,
`amplificadores.md`, `fotodetectores.md`, `mux_demux.md`) e nas notas de trabalho, com o
contexto para responder. Cada uma tem: o que é, como está hoje no código, as opções e o que
cada uma muda, e a minha recomendação. **Responda na linha "Resposta:"** de cada item (pode
ser só "A", "sim", "recomendação", ou um comentário).

As seções 1 e 2 travam o resto: definem como os blocos se ligam. Vale responder primeiro.

---

## 1. Interfaces entre os blocos

### 1. Convenção de frequência (M1 / X2)

**O que é.** O campo óptico real é $\mathrm{Re}[A(t)\,e^{\mp i\omega_0 t}]$; o simulador só guarda o
envelope $A(t)$. O sinal escolhido no expoente decide o que é "frequência positiva" no envelope.
Com $e^{-i\omega_0 t}$ (convenção do Agrawal, física), uma componente $e^{+i2\pi ft}$ do envelope (o
bin $+f$ da FFT do NumPy) está fisicamente em $f_0 - f$: lado **vermelho**. Com $e^{+i\omega_0 t}$
(convenção de engenharia e DSP), está em $f_0 + f$: lado **azul**.

**Como está hoje.** Os blocos foram escritos seguindo fontes diferentes, e as fontes discordam:

| Bloco | Fonte | Convenção |
|---|---|---|
| Fibra (`kernels/ssfm.py`) | Agrawal, Eq. 2.6.18 (cap. 2) | física: +f = vermelho |
| Laser CW (`kernels/laser/cw.py`) | `single_frequency_laser` do `draft.py` | desvio como $e^{+i2\pi\Delta f t}$: +f = azul |
| Laser de modulação direta (`components/tx/laser.py`) | Agrawal, Eq. 3.4.1 (cap. 3) | $\sqrt{P}e^{+i\phi}$ com δν > 0 = azul |
| Modulador IQ (`kernels/modulator/mzm.py`) | — | θ é parâmetro; professor +90°, VPI −90° |
| Mux (`kernels/wdm.py`) | professor | desloca com $e^{+i2\pi ft}$ |

O próprio Agrawal troca o sinal entre o cap. 2 e o cap. 3, e o professor também (`draft.py` usa
`+Freq_offset`, o notebook 3 usa `−Freq_offset`).

**Por que importa.** Cada bloco está certo sozinho. Quando o laser alimenta a fibra: (i) um canal
configurado em +Δf aparece do outro lado da grade; (ii) o chirp do laser de modulação direta chega
à fibra com o sinal trocado. Chirp × dispersão é o efeito que limita o IMDD com modulação direta:
o pulso que devia alargar comprime, ou o contrário.

**Opções.**
- **A. Convenção do Agrawal em todo o projeto.** Muda: `cw.py` (desvio $e^{-i2\pi\Delta ft}$), o bloco
  de modulação direta ($\sqrt{P}e^{-i\phi}$), o sinal de θ no IQ (combinar com o receptor), o mux.
  Poucas linhas. Todas as formas fechadas da fibra (cap. 2, 6 e 9), os testes de física e o
  notebook 01 continuam valendo como estão. Para comparar com o VPI, conjuga-se o campo exportado.
- **B. Convenção de engenharia.** Muda: a equação da fibra no `ssfm.py`
  ($\partial_zA = -\tfrac{\alpha}{2}A + \tfrac{i\beta_2}{2}\partial_T^2A - i\gamma|A|^2A$) e a definição do chirp C
  da gaussiana nos testes de física e no notebook 01. Os blocos de TX ficam como estão. Bate com o
  VPI e com textos de DSP sem conversão; as fórmulas do Agrawal da fibra passam a valer para o
  campo conjugado.

**Minha recomendação:** A (o Agrawal é a fonte principal e a fibra é a parte com mais fórmulas).

**Resposta:** (Gabriel, 2026-10-22) convenção vermelha do Agrawal (opção A). Aplicado: CW, laser direto e mux passam a usar a convenção do Agrawal; `TimeGrid.optical_frequency` dá a frequência óptica relativa de cada bin.

### 2. Quem põe o canal na sua frequência, e o que é `f0` (X1)

**O que é.** No WDM todos os canais vivem num só campo, em banda base em torno do centro da banda
(`SystemConfig.center_frequency`), cada um deslocado do seu desvio (`channel_offsets()`: −112,5,
−37,5, 37,5 e 112,5 GHz para 4 canais a 75 GHz). Alguém aplica esse deslocamento.

**Como está hoje.**
- Professor (`draft.py`): o laser de cada canal já sai deslocado (`Freq_offset`); o mux só soma.
- VPI: cada transmissor tem `EmissionFrequency`; o mux filtra e soma.
- Arquitetura: "multiplexer: soma N canais na grade"; e `Waveform.f0` é "o offset em relação ao
  centro da banda de simulação".
- Meus blocos de laser: põem em `f0` a frequência **absoluta** do canal (193,4 THz) e usam essa
  frequência para hν. Inconsistente com a arquitetura: erro meu, a corrigir de qualquer jeito.

**Opções.**
- **(a) O laser desloca** (`CWLaser(frequency_offset = desvio do canal)`) e o mux só soma, como o
  professor e o VPI. O mux fica trivial.
- **(b) O transmissor trabalha em banda base** (com o desvio guardado em `f0`) e o mux desloca cada
  canal para a grade. O TX fica independente da grade e o mux concentra a lógica da grade (onde a
  C-20 pede o registro da grade).

Nas duas, `f0` passa a ser o desvio, e a frequência absoluta (para hν, λ) sai da configuração.

**Minha recomendação:** (b), por ser a da arquitetura.

**Resposta:** (Gabriel, 2026-11-12) recomendação: o mux desloca (opção b). Aplicado: os transmissores ficam em banda base do canal e gravam o desvio do canal em `f0`; o mux confere esse `f0`, desloca cada canal e soma; depois do mux `f0 = 0`.

### 3. Potência de lançamento: total ou por canal, e quem a impõe (D-22 / X6)

**O que é.** O VPI tem um `PowerdBm` e não se sabe se é a potência total dos 4 canais ou por
canal. A D-22 registra duas verificações independentes (ruído do receptor e coeficiente não linear)
que apontam para **total**; falta confirmar no esquemático do VPI. A diferença é de 6 dB por canal.

**Como está hoje.** `IMDDChannel`/`CoherentChannel.launch_power` (por canal, opcional);
`AmplifierConfig.output_power` ("W, total do campo agregado"). O mux do professor soma sem
normalizar, e a potência total cresce com o número de canais sem ninguém declarar (C-20). A
arquitetura (Anexo A.4) termina o `build_tx` com um bloco `PowerControl(cfg.launch_power)`.

**Perguntas.** (i) Você consegue conferir no esquemático do VPI se `PowerdBm` é total? (ii) Quem impõe
a potência: um `PowerControl` por canal no fim do TX, ou o booster controlado por potência (como o
`AmpSysOpt` do VPI, `PowerControlled`)?

**Minha recomendação:** potência por canal na configuração do canal; o booster, se houver,
controla a potência total = soma das potências dos canais.

**Resposta:** (Gabriel, 2026-10-22) potência no amplificador por canal. Aplicado no comentário de `AmplifierConfig.output_power` (o total é N vezes); o bloco ainda não existe.

### 4. Frequência do laser de modulação direta no estacionário

**O que é.** Nas equações de taxa a fase gira com $d\phi/dt = \tfrac{\alpha_H}{2}[g_c(N - N_{tr}) - 1/\tau_p]$.
Com emissão espontânea e compressão de ganho, isso não é zero no estacionário: o laser polarizado
emite deslocado da frequência de referência, e o desvio depende da corrente. Com o ε do Agrawal:
−0,64 GHz em $I_{th}$, +0,55 GHz em 1,2·$I_{th}$, +2,8 GHz em 2·$I_{th}$, +5,9 GHz em 3·$I_{th}$, +12,7 GHz
em 5·$I_{th}$ (comparar com 75 GHz entre canais).

**Fato.** O Agrawal (texto da Eq. 3.4.1) define o chirp como a mudança da frequência "em relação ao
valor de estado estacionário ν₀".

**Opções.**
- **(a) Não descontar.** A frequência configurada é a da cavidade "fria"; a portadora sai alguns GHz
  fora, conforme a polarização.
- **(b) Descontar o desvio do estacionário em $I_b$.** A frequência configurada é a de emissão na
  polarização (que é como se especifica um canal da grade); o chirp aparece só como variação em
  torno dela.

**Minha recomendação:** (b).

**Resposta:** (Gabriel, 2026-10-22) "considere". Interpretei como adotar a recomendação (b): o bloco do laser direto desconta o desvio de frequência do estacionário na corrente de polarização. **Confirmar.**

### 5. Estrutura do transmissor IMDD

**O que é.** Em `config/system.py`: `IMDDTransmitter(frontend: DirectModulation | ExternalModulation,
rolloff)`. Hoje `DirectModulation.laser` é um `CWLaserConfig`, que não descreve modulação direta (o
bloco `DirectModulatedLaser` usa `RateEquationLaserConfig`, com 13 parâmetros sem padrão). E as
razões de polarização e de modulação ($I_b/I_{th}$, $I_m/I_{th}$) ainda não têm lugar na configuração.

**5a.** Trocar `DirectModulation.laser` para `RateEquationLaserConfig` e acrescentar `bias_ratio` e
`modulation_ratio`? (Sim, para o bloco funcionar a partir da configuração.)

**5b.** Com isso, `IMDDChannel()` sem `tx` deixa de ser possível (o padrão de `frontend` é
`DirectModulation`, que passa a exigir os 13 parâmetros). Isso quebra 2 testes de
`config/tests/test_system.py` (linhas 14 e 21), que criam `IMDDChannel` sem `tx`.
- **A.** `tx` obrigatório no IMDD; os 2 testes passam a informar o transmissor (mudança de
  especificação: precisa do seu aval, regra 1 do CLAUDE.md).
- **B.** O padrão do IMDD vira modulação externa (CW + MZM), que tem padrões.

**5c.** Seu `io_test.py` propõe `tx = {laser, mzm}` direto, sem `frontend`. A união atual
(`direct | external`) faz uma configuração inválida ser irrepresentável: um canal de modulação
direta não tem campo `mzm` (seção 4.6.1 da arquitetura). Com `tx = {laser, mzm}`, um canal de
modulação direta com `mzm` preenchido ficaria ambíguo. Manter a união ou mudar?

**Minha recomendação:** 5a sim; 5b A (coerente com a tarefa 1: campo sem base fica obrigatório);
5c manter a união.

**Resposta:** (Gabriel, 2026-11-12) manter. Aplicado: 5a sim (laser de equações de taxa com `bias_ratio` e `modulation_ratio`, sem padrão); 5b opção B (padrão do IMDD = modulação externa, sem mudar testes); 5c manter a união.

### 6. Chaves de ruído com dois donos (F4)

**O que é.** `SimulationConfig` tem `shot_noise` e `thermal_noise`; `PhotodiodeConfig` também tem
`shot_noise`. A regra do projeto é um dono por grandeza. Qual fica?

**Opções.** (a) O dispositivo é dono (`PhotodiodeConfig.shot_noise`, e ruído térmico desligado com
densidade 0); sai de `SimulationConfig`. (b) `SimulationConfig` vira uma chave geral de depuração
que desliga todos os ruídos (laser, ASE, fotodetector) de uma vez, sem duplicar a chave do
dispositivo.

**Minha recomendação:** (a), e se quiser a chave geral, ela com outro nome (`noise_enabled`).

**Resposta:** (Gabriel, 2026-11-12) manter as duas chaves separadas (sem a chave geral). Aplicado: o ruído do fotodetector entra se o dispositivo e a simulação o pedem.

---

## 2. Laser (`docs/laser.md`)

### 7. Ruído de Langevin com poucos fótons

**O que é.** As forças de Langevin são gaussianas, o que vale com muitos fótons na cavidade. Abaixo
do limiar há menos de um fóton (S ≈ 1e15 m⁻³), e um passo de ruído pode deixar S negativo; além
disso $D_{\phi\phi} = R_{sp}/(4S)$ diverge em S = 0 (laser partindo de zero). Implementei como
proposta: S limitado a ≥ 0 depois de cada passo, e $D_{\phi\phi}$ calculado com S ≥ $R_{sp}\tau_p$
($n_{sp}$ fótons). Acima do limiar isso não atua.

**Opções.** Manter; ou desligar o ruído abaixo de um número mínimo de fótons; ou outra regra.

**Minha recomendação:** manter.

**Resposta:** (Gabriel, 2026-10-22) manter.

### 8. Fase inicial do laser de modulação direta

**O que é.** A fase absoluta de um laser livre é desconhecida. Para o CW você decidiu fase inicial
uniforme em [0, 2π) (2026-10-08). No laser de modulação direta implementei igual quando há ruído, e
φ(0) = 0 sem ruído (para os casos deterministas do notebook). Importa no WDM: a FWM depende das
fases relativas entre canais.

**Minha recomendação:** manter.

**Resposta:** (Gabriel, 2026-10-22) com vários canais, fase inicial aleatória pela semente do canal. Já é assim (stream `ctx.rng_for(canal, "tx_laser_noise")`); sem ruído, φ(0) = 0.

### 9. Limitar o ganho a G ≥ 0

**O que é.** O notebook do professor usa `max(G, 0)`; o Agrawal (Eq. 3.3.3) não limita. G < 0 acontece
abaixo da transparência ($N < N_{tr}$) e representa absorção, que é física. Só muda o regime abaixo
da transparência, por exemplo nos primeiros ~0,3 ns ao ligar o laser a partir de N = 0. Hoje o
kernel **não** limita.

**Minha recomendação:** não limitar (Agrawal).

**Resposta:** (Gabriel, 2026-10-22) "sim". **Não apliquei**: "sim" responde "limitar a G ≥ 0" (professor), mas a resposta 22 diz "considere sempre o Agrawal", que não limita. Qual dos dois?

---

## 3. Moduladores (`docs/moduladores.md`, notebook 03)

### 10. O que é `vbias` (M2)

**O que é.** No push-pull (Agrawal, Eq. 10.1.5, e professor), $T = \cos^2[\pi(2V - V_b)/(2V_\pi)]$:
$V_b = V_\pi/2$ é a **quadratura** (ponto do OOK, T = 1/2 sem sinal) e $V_b = V_\pi$ é o **nulo** (BPSK e
cada braço do QAM, campo bipolar).

**Como está hoje.** `MZMConfig(vpi=5, vbias=2.5)`: quadratura. Serve à modulação externa do IMDD, mas o
`CoherentTransmitter` usa o mesmo padrão, e QAM na quadratura não funciona. O professor, no QAM, usa
$V_\pi = V_{bias}$ = 2,5 V (nulo). Seu `io_test.py` usa `vpi: 3.14, vbias: 1.57` no canal coerente
(quadratura). O módulo do VPI não tem parâmetro de polarização (o modelo já fica no nulo).

**Opções.** (a) `vbias` continua parâmetro com esse significado, e uma regra de validação exige
`vbias = vpi` no coerente. (b) No coerente a polarização sai do formato (nulo), sem campo; no IMDD
externo fica o `vbias`.

**Minha recomendação:** (b).

**Resposta:** (Gabriel, 2026-11-16) opção (a). Aplicado: `vbias` continua sendo parâmetro, e a regra de validação `coherent_mzm_at_null` (erro) exige `vbias = vpi` nos canais coerentes. Os exemplos coerentes ganharam `vbias: 5.0`. O bloco `IQModulator` continua supondo o nulo e não lê o campo.

### 11. Escala do sinal elétrico em volts (M3)

**O que é.** O sinal elétrico depois da formatação é adimensional (símbolos ±1, ±3...). O MZM
precisa de tensão. A excursão decide o compromisso: mais excursão dá mais potência óptica e mais
não linearidade (no notebook 03, seção 6, o 16QAM com excursão inteira comprime: razão externo/interno
2 em vez de 3).

**Como está hoje.** Nenhum campo na configuração. O professor usa $V = K_{mod}\cdot s$ com
$K_{mod}$ = 0,1 V e s normalizado pelo máximo, com $V_\pi$ = 2,5 V: fase de pico de 0,126 rad (muito
linear, pouca potência). No notebook 3 ele usa ±s/2 por braço.

**Proposta.** Campo `drive_amplitude` em múltiplos de $V_\pi$ (pico por braço / $V_\pi$), como fiz com o
laser de modulação direta em múltiplos de $I_{th}$. Sem padrão (sem base na literatura).

**Resposta:** (Gabriel, 2026-11-12) "não sei". Mantido `ExternalModulation.drive_amplitude` = 0,25 (pico por braço / Vπ) como padrão: com V_b = Vπ/2 é a excursão de nulo a pico, a de maior extinção (Eq. 10.1.5). O professor só usa sinal pequeno no coerente (`ind_mod` = 0,1); o VPI não tem IMDD.

### 12. Razão de extinção do MZM (M4)

**O que é.** O MZM do Agrawal é ideal (extinção infinita). A configuração e o VPI têm ER = 30 dB e
perda de inserção de 6 dB. O professor obtém extinção finita com perdas diferentes nos braços:
ER = $[(a_1 + a_2)/(a_1 - a_2)]^2$ (com os padrões dele, 2 e 3 dB, dá 24,8 dB). **[Dedução]** Mapeando:
$a_1 + a_2 = 10^{-IL/20}$ e $(a_1 - a_2)/(a_1 + a_2) = 10^{-ER/20}$. Efeito colateral físico: chirp
residual perto do nulo e vazamento de portadora no IQ (notebook 03, seção 3).

**Minha recomendação:** desbalanço dos braços, como o professor.

**Resposta:** (Gabriel, 2026-10-22) manter o do professor (desbalanço dos braços). Para o bloco do modulador.

### 13. Filtro S21 do modulador (M5)

**O que é.** A banda elétrica do modulador (eletrodos e driver): 40 GHz, ordem 3 no VPI. O tipo do
filtro não está no arquivo do VPI.

**Opções.** Bessel (fase linear), Butterworth ou gaussiano, aplicado à tensão de cada braço antes do
MZM. Dá para conferir no VPI qual é o filtro do `S21_Bandwidth` do `Tx_mQAM`?

**Minha recomendação:** Bessel de ordem 3, salvo se o VPI disser outro.

**Resposta:** (Gabriel, 2026-10-22) Bessel de ordem 3; o VPI assume banda infinita. Para o bloco: S21 Bessel de ordem 3, e sem `s21_bandwidth` (padrão) a banda é infinita, como no VPI.

### 14. `PreDistortion = Yes` do VPI (M6)

**O que é.** O transmissor do VPI tem pré-distorção ligada. **[Dedução, não verificada]** Deve
aplicar arco-seno ao sinal elétrico para compensar o seno do MZM, deixando o campo linear nos
símbolos (no notebook 03 se vê a compressão do 16QAM sem isso). O roadmap já pedia para confirmar.

**Pergunta.** Confere no VPI o que ela faz? E implementamos (como opção do bloco)?

**Resposta:** (Gabriel, 2026-10-22) sim, tem pré-distorção. Aplicado: `arcsin_predistortion` no kernel do MZM e painel no notebook 03.

### 15. Skew de 0,75 ps entre I e Q (M7)

**O que é.** O VPI atrasa o sinal I em relação ao Q em 0,75 ps. Com a grade de 16 amostras por
símbolo a 59,84 GBd (passo de 1,04 ps) é uma fração de amostra. Proposta: atraso por fase linear
na FFT ($e^{-i2\pi f\tau}$), exato para sinais limitados em banda. A D-4 proíbe reamostragem
fracionária; isto é um filtro de atraso, não reamostragem. Está ok?

**Resposta:** (Gabriel, 2026-10-22) considere (fase linear na FFT). Para o bloco do modulador.

### 16. Dupla polarização no transmissor (M8)

**O que é.** O professor divide a potência do laser meio a meio entre X e Y (`X_fraction = 0.5`),
modula cada polarização com um IQ e recombina. O VPI (`Tx_mQAM_DP`) tem a mesma estrutura e ainda
`SplitterExtinctionRatioX/Y` e `CombinerExtinctionRatioX/Y` = 30 dB, que não sei se são a extinção
de polarização do divisor/combinador (vazamento entre X e Y) ou dos acopladores do IQ.

**Perguntas.** Meio a meio, como o professor? E dá para conferir no VPI o que são esses 30 dB?

**Resposta:** (Gabriel, 2026-10-22) "o do VPI é isso mesmo, use dele". Aplicado: divisão meio a meio entre X e Y, como o VPI. **Falta**: os 30 dB de `Splitter/CombinerExtinctionRatio` são vazamento entre X e Y (divisor e combinador de polarização) ou dos acopladores do IQ?

### 17. Potência de saída do modulador (M9)

**O que é.** O VPI tem `Renormalize = On` (provavelmente renormaliza a saída para a potência do
laser; inferência). As opções são saída física (laser × perda de inserção × perda de modulação, que
depende da excursão) ou renormalizada. O roadmap pede "potência média da saída bate com laser ×
perdas".

**Minha recomendação:** saída física; a potência de lançamento é imposta depois (pergunta 3).

**Resposta:** (Gabriel, 2026-10-22) renormalizar. Interpretei: a potência média de saída do modulador volta à potência do laser (`AveragePower` do VPI). **Confirmar** o alvo.

### 18. Desbalanço de ganho I/Q

**O que é.** O VPI tem no TX só desbalanço de fase (5°) e skew; no RX tem
`InsertionLossImbalanceI/Q` (= 0). Seu `io_test.py` tem `iq.gain_imbalance`; o `IQImpairments`
não tem. Acrescentar (padrão 0, neutro)?

**Resposta:** (Gabriel, 2026-10-22) considerar o desbalanço de ganho, padrão 0. Aplicado: `IQImpairments.gain_imbalance_db = 0.0` e parâmetro `q_gain` no `iq_field`.

---

## 4. Amplificador (`docs/amplificadores.md`, notebook 04)

### 19. Figura de ruído → $n_{sp}$ (A1)

**O que é.** Agrawal, Eq. 7.2.15: $F_n = 2n_{sp}(1 - 1/G) + 1/G \approx 2n_{sp}$. Invertendo exatamente,
$n_{sp} = (GF_n - 1)/[2(G - 1)]$, e a NF medida fica igual à configurada em qualquer ganho (notebook
04, seção 4c). O professor usa $n_{sp} = F_n/2$: a NF real cai em ganho baixo (1,8 dB abaixo em
G = 3 dB; igual acima de ~15 dB). Nos trechos de 80–120 km (G de 16 a 24 dB) a diferença é < 0,1 dB.

**Minha recomendação:** a exata.

**Resposta:** (Gabriel, 2026-10-22) Agrawal (Eq. 7.2.15 invertida). Já é o kernel.

### 20. ASE por polarização (A2)

**O que é.** Agrawal, Eq. 7.4.1 e 7.5.10: a ASE total é $2S_{ASE}\Delta\nu$, "o fator 2 leva em conta a
natureza despolarizada da ASE"; ou seja, $S_{ASE}$ em **cada** polarização. O código do professor
divide $S_{ASE}$ entre X e Y (total $S_{ASE}$): metade do Agrawal, OSNR 3 dB melhor.

**Minha recomendação:** Agrawal (é também como se define a OSNR em 0,1 nm na prática).

**Resposta:** (Gabriel, 2026-10-22) Agrawal (ASE em cada polarização). Já é o kernel.

### 21. Campo de uma polarização só (A3)

**O que é.** Os lasers do pacote devolvem (2, N) com Y = 0, então a ASE ortogonal tem onde entrar.
Se chegar um campo (1, N), somar ASE só nessa polarização subestima a ASE total (falta a metade
ortogonal, que entra no batimento ASE-ASE e na potência; o sinal-ASE não muda).

**Opções.** O bloco exige (2, N) e dá erro; ou aceita (1, N) com ASE só nessa polarização.

**Minha recomendação:** exigir (2, N) no bloco.

**Resposta:** (Gabriel, 2026-10-22) "caso chegue (1, N), divida o ASE e na outra polarização". Interpretei: o campo vira (2, N) com Y = 0 e a ASE entra nas duas polarizações. Aplicado no `amplify`. **Confirmar.**

### 22. Modo `target_output_power` (A4)

**O que é.** O amplificador do VPI é controlado por potência: a saída fica em `OutputPower` (a
potência de lançamento). Um amplificador real faz $G = P_{out}/P_{in}$ com a potência **total** de
entrada, que inclui a ASE acumulada dos anteriores: em cadeias longas o sinal ganha um pouco menos
que a perda do trecho ("a ASE rouba ganho"). Alternativa: G calculado só com a potência do sinal.
O VPI também tem `SaturationPower` = 10 mW e `GainMax` = 100 dB (modelo de saturação do ganho).

**Perguntas.** G pela potência total ou só do sinal? A saturação entra agora?

**Minha recomendação:** potência total; sem saturação por enquanto. Ver também a pergunta 24.

**Resposta:** (Gabriel, 2026-10-22) considerar sempre o Agrawal; se houver conflito com o que definimos para o sistema, vale o sistema. Interpretei para o bloco: modos e potência como na configuração (`target_output_power`, por canal: resposta 3), ganho aplicado ao campo total (sinal + ASE), sem saturação (a configuração não tem). **Confirmar.** (Gabriel, 2026-11-12) saturação: "pode colocar". Aplicado: `AmplifierConfig.saturation_power`, g = g0/(1 + P/Ps) (Agrawal, Exercício 7.2), ganho pela potência média de entrada, no modo `fixed_gain`; no controlado por potência é recusada (a confirmar).

### 23. Frequência para hν (A5)

**O que é.** $S_{ASE} \propto h\nu$. Na banda WDM (±150 GHz em torno de 193,7 THz) hν varia ±0,08 %:
desprezível. A ASE é somada ao campo agregado inteiro.

**Minha recomendação:** a frequência central da simulação (`SystemConfig.center_frequency`).

**Resposta:** (Gabriel, 2026-11-16) "pode ser a sua recomendação". Já estava assim: $h\nu$ na frequência central da simulação.

### 24. Banda do ruído de ASE (A6)

**O que é.** A ASE é somada com densidade $S_{ASE}$ na banda inteira da grade ($f_s$), como o
professor. O VPI guarda o ruído fora da banda amostrada à parte ("noise bins"). Medidas filtradas
(OSNR em 0,1 nm, receptor depois do demux) não dependem de $f_s$; mas a potência total do campo
inclui ASE fora da banda dos canais, e um amplificador controlado por potência total (pergunta 22)
passaria a depender de $f_s$.

**Opções.** ASE na banda inteira (e, no controle por potência, usar a potência total mesmo assim, ou
só a potência dentro da banda dos canais); ou limitar a ASE a uma banda configurável.

**Minha recomendação:** banda inteira; no controle por potência, usar a potência na banda ocupada
pelos canais.

**Resposta:** (Gabriel, 2026-11-16) "pode ser a sua recomendação". Aplicado: a ASE continua na banda inteira da grade, e o amplificador controlado por potência passou a medir a entrada só na banda ocupada pelos canais (`Amplifier.signal_band`; `pipelines/imdd.py`, `signal_band`: do canal mais baixo ao mais alto, com margem de meio espaçamento ou de uma taxa de símbolos, a que for maior). Nenhum exemplo mudou de resultado.

---

## 5. Fotodetector (`docs/fotodetectores.md`, notebook 05)

### 25. Densidade do ruído (F1)

**O que é.** Agrawal, Eq. 4.4.2 e 4.4.7: as densidades espectrais **bilaterais** são $qI$ (disparo) e
$2k_BT/R_L$ (térmico). Numa corrente amostrada a $f_s$, a variância por amostra é a densidade bilateral
× $f_s$. O professor usa a unilateral × $f_s$: dobra os dois ruídos (SNR 3 dB pior). O kernel segue o
Agrawal, e o notebook 05 confere contra as Eq. 4.4.5 e 4.4.9 (erro < 4 %, estatístico).

**Minha recomendação:** Agrawal.

**Resposta:** (Gabriel, 2026-11-16) "pode ser a sua recomendação". Já estava assim: densidades bilaterais do Agrawal.

### 26. Densidade do ruído térmico (F2)

**O que é.** `PhotodiodeConfig.thermal_noise_density` = 10 pA/√Hz (o `ThermalNoise` do VPI). Assumi
que é densidade **unilateral** de corrente ($\sigma_T^2 = i_n^2\Delta f$), a convenção usual de A/√Hz.
O professor calcula $4k_BT/R_L$ com 300 K e 50 Ω: 18,2 pA/√Hz. A proposta da tarefa 1 tinha 5,76
pA/√Hz (Agrawal: $R_L$ = 1 kΩ, $F_n$ = 2, T = 300 K suposto).

**Perguntas.** A convenção unilateral está certa? E o padrão: 10 pA/√Hz (VPI, mas a tarefa 1 diz que
valor do VPI não é padrão), a conta do Agrawal, ou obrigatório?

**Resposta:** (Gabriel, 2026-11-16) "faça uma comparação de valores entre as opções, mas eu iria de Agrawal". Aplicado: o padrão de `thermal_noise_density` passou de 10 para 5,76 pA/√Hz (Agrawal, Eq. 4.4.9, com $R_L$ = 1 kΩ, $F_n$ = 2 e T = 300 K). Os exemplos continuam com os 10 pA/√Hz do VPI, escritos no YAML. Comparação (convenção unilateral nas três):

| Opção | Densidade [pA/√Hz] | IMDD, OOK a 10 Gb/s: potência recebida para Q = 6 | Coerente, QPSK a 10 GBd e −41 dBm: SNR |
|---|---|---|---|
| Agrawal ($R_L$ = 1 kΩ, $F_n$ = 2, 300 K) | 5,76 | −24,08 dBm | 11,59 dB |
| VPI (`ThermalNoise`) | 10,00 | −21,71 dBm (2,36 dB pior) | 9,61 dB |
| Professor ($4k_BT/R_L$, 50 Ω, 300 K) | 18,20 | −19,13 dBm (4,94 dB pior) | 6,05 dB |

Sobre a convenção: com a unilateral, o ruído do receptor do simcopy fecha com o do VPI em 0,1 dB (`docs/comparacao_vpi.md`, seção 4.4).

### 27. Filtro do receptor (F3)

**O que é.** O VPI usa Bessel de ordem 4 com 1,15 × taxa de símbolos (o cenário é coerente). O kernel
aplica a resposta analógica na FFT, com o atraso de grupo (~34 ps com 10 GHz; ~5 ps com 69 GHz), que
desloca o instante ótimo de amostragem (o DSP ou a sincronização corrigem). O professor usa um
elíptico com `lfilter` (P-15). No IMDD a escolha clássica é outra: banda de ruído ≈ B/2 (Agrawal), ou
~0,75 × taxa de bits.

**Perguntas.** Bessel aplicado na FFT está ok? No IMDD, `bandwidth_factor` também, sem padrão?

**Resposta:** (Gabriel, 2026-11-16) "eu iria de Bessel, mas gostaria de manter um padrão, tanto IMDD quanto coerente [...] você acha ser a melhor estratégia essa do VPI?". Como está: Bessel de ordem 4 aplicado na FFT é o padrão dos dois receptores. Minha opinião, com os números abaixo: o tipo e a ordem do VPI são uma boa escolha, a banda de 1,15 vezes a taxa não é. Medido em 2026-11-16 (Bessel de ordem 4):

| Banda / taxa de símbolos | IMDD, OOK a −20 dBm: Q | Coerente, QPSK a −41 dBm, com DSP: SNR [dB] |
|---|---|---|
| 0,50 | 10,04 | 12,70 |
| 0,65 | 10,91 | 12,57 |
| 0,75 | 10,24 | 12,46 |
| 1,00 | 8,90 | 11,99 |
| 1,15 (VPI) | 8,32 | 11,60 |
| 1,50 | 7,40 | 10,60 |

No IMDD limitado por ruído, 1,15 custa 2,3 dB em Q contra 0,65. No coerente com DSP custa 1,1 dB de SNR contra 0,5. Um Butterworth de ordem 4 dá até 0,7 dB a mais que o Bessel no coerente com DSP (o equalizador corrige a fase) e é pior no IMDD com banda estreita. **Falta decidir:** o padrão de `bandwidth_factor` (hoje 1,15, o do VPI; eu poria 0,75 e deixaria 1,15 só no YAML do VPI), e se quer a comparação com os outros simuladores (pede instalar o OptiCommPy, que é GPL, num ambiente separado: aguardo o seu aval).

### 28. Saída do fotodetector (F5)

**O que é.** Proposta: a corrente em A, num `Waveform` elétrico (1, N) com parte imaginária nula
(D-5). Alternativa: tensão depois de um amplificador de transimpedância (que ainda não modelamos).

**Minha recomendação:** corrente em A.

**Resposta:** (Gabriel, 2026-11-16) "Sim, para isso que existe o Waveform". Já estava assim: corrente em A, num `Waveform` elétrico.

---

## 6. Mux e demux (`docs/mux_demux.md`, notebook 06)

### 29. Filtros no mux (X3)

**O que é.** O mux do VPI filtra cada entrada (passa-banda "Gaussian" de ordem 3, 70 GHz) e tem 3 dB
de perda. O do professor é uma soma ideal. Filtrar no mux corta o que o transmissor emite fora da
banda (por exemplo o espectro alargado pelo chirp da modulação direta) e distorce um pouco o canal;
a perda reduz a potência lançada (o booster compensa).

**Opções.** Soma ideal; ou filtro e perda configuráveis por porta (como o VPI).

**Resposta:** (Gabriel, 2026-11-12) as duas. Aplicado: `SystemConfig.mux` (`MuxConfig`), soma ideal por padrão; com `optical_filter` e `insertion_loss_db`, filtro por canal e perda, como o WDM_MUX_N_1 do VPI.

### 30. Filtro do demux (X4)

**O que é.** Máscara retangular (professor): devolve o canal exato quando os espectros não se
sobrepõem (notebook 06, seção 2), mas tem resposta ao impulso tipo sinc (*ringing*, P-8 da
arquitetura). Supergaussiana (C-11 da arquitetura; o VPI usa "Gaussian de ordem 3", 70 GHz).
Implementei a supergaussiana como $\exp[-\tfrac{\ln 2}{2}(2(f - f_c)/B)^{2n}]$, com B a banda de 3 dB:
é a minha leitura da definição do VPI. Com ordem 3 e 70 GHz: distorção de 1,8 % no canal e
penalidade de interferência de 0,0006 dB; com ordem 1 (gaussiana): 0,54 dB.

**Perguntas.** Qual filtro? E dá para conferir no VPI a definição do "Gaussian de ordem N" (B é a
banda de 3 dB?)?

**Resposta:** (Gabriel, 2026-11-12) "não sei". Estudo (4 canais a 50 GHz, `docs/sistema_imdd.md`, seção 8): supergaussianos de ordem 3 entre 20 e 50 GHz dão Q de 27 a 35; abaixo de 15 GHz o filtro corta o sinal (Q ~18); retangular de 50 GHz (professor) dá 30,5; retangular de 25 GHz, 20,6 (ringing). Mantido sem padrão; o exemplo usa supergaussiano de ordem 3 com 30 GHz (3× a taxa de bits).

### 31. Saída do demux (X5)

**O que é.** Devolver o canal na posição dele na grade (o `optical_filter` do professor faz isso) ou
trazê-lo para a banda base. No IMDD não importa (o detector vê |E|²). No coerente importa para o
oscilador local: se o canal fica na posição, o LO é um `CWLaser` no desvio do canal mais o
`lo_frequency_offset` (como o professor: `Freq_offset_RX = Freq_offset_TX + offset`).

**Minha recomendação:** manter na posição (fisicamente o demux não muda a frequência; quem faz isso é
o LO).

**Resposta:** (Gabriel, 2026-11-12) manter na posição. Aplicado.

---

## 7. Fibra (notebook 01)

### 32. Campo ao longo de z

**O que é.** `propagate` devolve só o campo no fim da fibra; para mostrar a evolução ao longo de z
(a aula pede "como ele varia ao longo da distância"), o notebook repropaga a cada posição. Opção:
um parâmetro do kernel que devolve o campo em posições escolhidas, no mesmo laço.

**Resposta:** (Gabriel, 2026-11-16) "pode ser". Aplicado: `propagate(..., positions=[...])` devolve também o campo em cada distância pedida, no mesmo laço, e o bloco ganhou `Fiber.along(x, ctx, positions)`. Conferido: sem `positions` nada muda; o campo em cada posição difere da propagação direta até lá em menos de 1e-6.

### 33. Campo de uma polarização ignora a PMD sem aviso

**O que é.** `propagate` só aplica PMD se o campo for (2, N). Um campo (1, N) com
`pmd_coefficient > 0` passa sem PMD, sem aviso (notebook 01, seção 4). Opções: erro, aviso, ou
aceitar (campo escalar = sem PMD por definição).

**Minha recomendação:** erro (ou uma regra de validação).

**Resposta:** (Gabriel, 2026-11-16) "aviso". Aplicado no bloco `Fiber`: com PMD ligada e campo (1, N), ele avisa (`warnings.warn`) e segue sem PMD. Nas montagens pela configuração isso não acontece: os lasers já emitem (2, N), então a PMD entra também no IMDD.

### 34. Passo máximo do SSFM com PMD

**O que é.** No modelo de PMD do pacote, cada passo aplica um DGD de módulo fixo
($D_p\cdot\text{passo}/\sqrt{\max(\text{passo}, l_c)}$) com orientação sorteada. Se a fibra inteira cabe
em um passo (fibra ≤ `ssfm_max_step` = 5 km por padrão), o DGD sai **sempre o mesmo**, sem a
distribuição maxwelliana (notebook 01, seção 4.1). Em trechos de 80 km com passo de 5 km (16 passos)
está ok.

**Opções.** Com PMD ligada, limitar o passo (por exemplo a 1/10 da fibra), por regra de validação
ou automaticamente; ou aceitar.

**Resposta:** (Gabriel, 2026-11-16) "limite o passo seguindo essa regra". Aplicado: com PMD, o passo do SSFM é no máximo 1/10 do comprimento da fibra (`PMD_MIN_STEPS` em `components/optical/fiber.py`). Efeito nos exemplos: PAM4 em 10 km, Q de 44,46 para 44,56; WDM de 4 canais (a compensadora tem 12,8 km), Q de 34,92 para 34,35; OOK em 20 km não muda (o passo já era menor). O cenário do VPI não muda.

### 35. β₃ (tarefa 2)

**O que é.** O kernel só tem β₂; a fibra do VPI tem `DispersionSlope` = 0,08 ps/(nm²·km). Pesa em
grades WDM largas e perto do zero de dispersão. Quando entra? (O sinal de β₃ no eixo da FFT depende
da pergunta 1.)

**Resposta:** (Gabriel, 2026-11-16) "por mim poderia entrar sempre". O $\beta_3$ entrou em 2026-11-16 (`fiber.dispersion_slope`). Se entra sempre ou só quando pedido é a dúvida 68, que ficou para discutir.

### 36. FWM não desligável (D-21)

**O que é.** O VPI tem `FWM = No`; o SSFM do campo agregado inclui a FWM entre canais e não há como
desligar. Opções: aceitar e documentar a diferença nas comparações; ou propagar os canais separados
com equações acopladas (muito mais trabalho).

**Minha recomendação:** aceitar e documentar.

**Resposta:** (Gabriel, 2026-11-16) aceitar. Registrado em `docs/comparacao_vpi.md`, seção 6.

### 37. `DispersionSlope` do bloco de compensação do VPI

**O que é.** No bloco de compensação de dispersão do DSP do VPI, `DispersionSlope = 0.08e-3`; na
fibra, `0.08e3` s/m³: seis ordens de diferença, provavelmente erro de digitação no arquivo. Só
importa se formos replicar a compensação do VPI com inclinação. Dá para conferir?

**Resposta:** (Gabriel, 2026-11-16) "é isso mesmo, foi erro de digitação". A compensação de dispersão do simcopy usa a mesma inclinação da fibra.

---

## 8. Configuração e tarefas antigas

### 38. Tarefa 0: reconciliar com o pacote `simcopy_1`

**O que é.** O pacote tem documentos que não estão no repo (`CLAUDE.md`, `docs/decisoes.md`,
`historico.md`, `roadmap.md`, `arquitetura.md` v0.9, `reference/vpi/parametros_vtmu.md`). Meus guias
citam vários deles, e o grupo não os vê. Você disse "dps vemos". Proposta: comparar e trazer esses
documentos para o repo (mostrando as diferenças antes de sobrescrever); o PDF do artigo fica fora.

**Resposta:** (Gabriel, 2026-11-16) "eu manteria o CLAUDE.md de fora por enquanto, o resto pode incluir". Aplicado: copiados do pacote `docs/decisoes.md`, `docs/historico.md`, `docs/roadmap.md`, `docs/arquitetura.md` e `reference/vpi/parametros_vtmu.md`. Nenhum existia no repo, nada foi sobrescrito. Ficaram de fora o `CLAUDE.md`, o PDF e o `docs/figuras.md` do pacote (o do repo é diferente e foi mantido).

### 39. Tarefa 1: valores padrão

**O que é.** Decidido: padrões de `src/` são genéricos, de literatura, com fonte; valores do VPI só
no YAML de exemplo e nos testes; campo sem base fica obrigatório. Falta: (i) confirmar a tabela
(`memory/tarefa1_padroes.md`: atenuação 0,2 dB/km, D 16, Aeff 80 µm², n₂ 2,6e-20, responsividade 1,
ruído térmico 5,76 pA/√Hz, `correlation_length` 10 m, frequência central 193,1 THz, espaçamento 50
GHz...); (ii) `correlation_length` 10 m (Agrawal, "~10 m") contra 50 m (VPI, D-18); (iii) T = 300 K
(suposição); (iv) testes que mudam de especificação (o de λ derivado espera 1547,715 nm com a
frequência padrão; testes que criam canais sem `tx`/`rx`; o YAML base sem os campos obrigatórios).

**Resposta:** (Gabriel, 2026-11-16) "acho que já temos dados o suficiente para resolver isso". O que os dados resolvem: (i) atenuação, dispersão, área efetiva, $n_2$ e responsividade ficam como estão (Agrawal e VPI coincidem); o ruído térmico passou a 5,76 pA/√Hz (resposta 26). (ii) `correlation_length`: a comparação com o VPI não distingue 10 m de 50 m. Recomendo 10 m (Agrawal) no padrão e 50 m no YAML do VPI, onde já está. (iii) T = 300 K deixa de ser suposição: o Agrawal usa "temperatura ambiente (T = 300 K)" (capítulo 3). (iv) O formato padrão do IMDD passou a OOK (resposta 51). **Não aplicado, falta o seu aval:** `correlation_length` 10 m, frequência central 193,1 THz, espaçamento de 50 GHz, `lo_frequency_offset` 0 e os campos obrigatórios. Mudam os números dos notebooks e quebram testes do grupo que usam os padrões de hoje, e se ligam à forma B de `docs/proposta_ideal_real.md`.

### 40. ADC: amostras por símbolo ou taxa absoluta

**O que é.** Hoje `ADCConfig(samples_per_symbol=2, resolution_bits=12)` (o VPI usa 2 × taxa de
símbolos), com a regra de validação que exige divisão inteira da taxa da grade. Seu `io_test.py`
usa `sampling_rate` absoluta (20 GHz, 64 GHz). A absoluta é mais parecida com a especificação de um
hardware, mas precisa da mesma verificação de divisão.

**Minha recomendação:** manter `samples_per_symbol`.

**Resposta:** (Gabriel, 2026-11-16) "mantenha". Mantido `samples_per_symbol`.

### 41. Potência no YAML: W ou dBm

**O que é.** Hoje em W. Opção: aceitar `power_dbm` como alternativa exclusiva, como já acontece com
`center_frequency` e `center_wavelength` (um ou outro, conferindo a coerência se vierem os dois).

**Resposta:** (Gabriel, 2026-11-16) "pode ser sua opção". Já existia para `launch_power`, `output_power` e `saturation_power`. Aplicado agora: `power_dbm` nos lasers (`tx.laser` e `rx.lo`).

### 42. "IMDD < 5w" nas notas da aula

**O que é.** Em `docs/04-class/context.md` há a linha "IMDD < 5w", antes da cadeia IMDD. Não sei se é
potência (5 W? 5 mW?), número de canais, de trechos, ou outra coisa. O que significa?

**Resposta:** (Gabriel, 2026-11-16) "é potência, se não me engano há uma limitação de 5 mW ou 52 para sistemas IMDD, pode ser até um aviso na montagem". **[Fato: Agrawal, seções 2.6.1 e 6.3]** É o limiar do espalhamento Brillouin estimulado: "da ordem de 5 mW para fibras longas (> 10 km)", e perto de 50 mW em 1 km. Aplicado: regra `imdd_power_below_brillouin` (aviso) quando um canal IMDD entra com mais de 5 mW numa fibra de mais de 10 km. O simulador não tem Brillouin. **Confirmar:** o "52" é os 50 mW das fibras curtas?

---

## 9. Pendências práticas

### 43. `io.py`, `_imdd`: `NameError`

Depois do seu renome de `c` para `channel`, o `_imdd` ainda usa `c["index"]`, `c["symbol_rate"]`
etc. Carregar um canal IMDD de YAML dá `NameError`. Nenhum teste cobre (o YAML de exemplo é
coerente). Corrijo trocando para `channel[...]`?

**Resposta:** (Gabriel, 2026-11-16) "sim, pode corrigir". Já estava corrigido: os exemplos IMDD carregam do YAML.

### 44. `io_test.py` importa `LaserConfig`

Agora se chama `CWLaserConfig`. O arquivo está fora do pytest (`addopts --ignore`), então nada quebra
hoje; atualizo quando o rascunho virar teste?

**Resposta:** (Gabriel, 2026-11-16) "pode atualizar já". Aplicado: `io_test.py` importa `CWLaserConfig`.

### 45. `core/tests/test_birefringence.py`

Arquivo novo (2026-10-16), fora do commit. Quebra a coleta do pytest inteiro: `from birefringence
import dgd_for_step` (o caminho certo é `simcopy.kernels.birefringence`), e a função de teste se chama
`dgd_for_step` (sobrescreve o import e não começa com `test_`). É seu ou de um colega (issue #10)?
Corrigir ou mover para `rascunhos/`?

**Resposta:** (Gabriel, 2026-11-16) "sim, mova, é de um colega". Movido para `rascunhos/test_birefringence.py`. O pytest volta a coletar sem precisar ignorar o arquivo.

### 46. Texto das issues #14 e #15

Elas dizem que os testes "já existem, com N testes", e eu os retirei. Troco para "a criar"?

**Resposta:** (Gabriel, 2026-11-16) "não, não precisa".

### 47. PR do branch `laser` para o `main`

O grupo trabalha a partir do `main`; os arquivos das issues #12 a #21 só aparecem lá depois do merge.
Abrir o PR agora? (O remoto ainda aponta para o nome antigo, `pycosim`; o GitHub redireciona.)

**Resposta:** (Gabriel, 2026-11-16) "já já vamos arrumar isso".

### 48. Issues para os testes de sistema

Na conversa sobre o que testar antes do sistema, sugeri três grupos: (a) regras de validação em
cascata para os blocos novos (passo do laser contra a grade, banda da grade contra o chirp, excursão
do MZM, parâmetros obrigatórios de cada modo do amplificador); (b) contratos entre blocos (formato e
domínio, `f0` e `channel_id`, contabilidade de potência, sementes independentes por canal); (c)
testes de ponta a ponta com forma fechada (BER de OOK contra o fator Q, como no notebook 05;
desvanecimento por dispersão; orçamento de potência; determinismo). Abro as issues?

**Resposta:** (Gabriel, 2026-11-16) "pode abrir assim que mergearmos isso, mas espere". Aguardando o merge.

---

## 10. Novas, do planejamento (`docs/planejamento.md`, 2026-10-25)

### 49. Sequência pseudoaleatória: ordem

**O que é.** A aula pede "criar sequência pseudoaleatória" (E1). O VPI usa `PRBS_Order = 7`
(período de 127 bits; `reference/vpi/parametros_vtmu.md` do pacote). O `decisoes.md` do pacote
registra o problema: na janela do VPI a sequência se repete ~8.258 vezes, o que cria correlação
artificial entre símbolos e entre canais (afeta XPM e a contagem de erros).

**Opções.** Ordem 7, como o VPI; ou ordem maior; ou ordem configurável. E, no WDM, sementes ou
deslocamentos diferentes por canal para não repetir a mesma sequência.

**Minha recomendação:** ordem configurável, sem padrão (campo sem base fica obrigatório, tarefa 1);
nos notebooks curtos, 7.

**Resposta:** (Gabriel, 2026-11-12) manter: bits aleatórios pela semente do canal, não PRBS. Aplicado; a semente do canal define os bits sozinha (resposta 8 do `docs/sistema_imdd.md`).

### 50. Sinal elétrico do IMDD: filtro de Bessel no lugar de Nyquist

**O que é.** A aula: "converter para sinal elétrico (não se usa Nyquist, e sim um filtro de
Bessel)". Hoje `IMDDTransmitter.rolloff = 0.2` (cosseno levantado, Nyquist), e a regra de
validação da banda da grade (`validation/rules.py`) usa esse `rolloff` para calcular a banda do
sinal.

**Perguntas.** (i) No IMDD, trocar `rolloff` por pulso NRZ retangular + Bessel, com ordem e banda
(em múltiplos da taxa de símbolos) na configuração? (ii) Ordem e banda: há valor da aula ou do
professor? Não achei fonte; na resposta 13 a ordem 3 é a do S21 do modulador. (iii) A regra da
banda da grade passa a usar a banda do Bessel?

**Resposta:** (Gabriel, 2026-11-12) "pode revisar e ver". Revisão: a aula pede Bessel sem número; o Agrawal não especifica o filtro do transmissor; o professor não tem IMDD (só o DAC de Nyquist do coerente); o VPI tem no transmissor `LPF_CutoffFrequency = BaudRate` e `LPF_FilterOrder = 4` (o tipo do filtro não está no arquivo). Adotado: Bessel de ordem 4 com −3 dB na taxa de símbolos, nos exemplos e no notebook 07. Os campos continuam obrigatórios no YAML.

### 51. Formato do IMDD

**O que é.** `IMDDChannel.modulation` tem padrão `PAM4`; o notebook 05 e a sequência da aula (E1 a
E9) trabalham com OOK. Fazer as semanas 1 a 3 com OOK e deixar o PAM4 para depois? E o padrão da
configuração muda para OOK ou fica obrigatório?

**Resposta:** (Gabriel, 2026-11-16) "OOK é o padrão, PAM4 fica como bônus". Aplicado: o padrão de `IMDDChannel.modulation` passou a OOK.

### 52. SMF + DCF + amplificador no trecho

**O que é.** A aula pede "duas fibras ao invés de uma (SMF, DCF)" e "um amplificador depois da
fibra DCF" (E6, E7). Hoje `SpanConfig(fiber, amplifier, line_compensation: bool)`: uma fibra só, e
a arquitetura trata a compensação de linha como um bloco que o `build_link` insere ou não, não
como uma fibra.

**Opções.** (a) `SpanConfig.fibers` vira uma tupla (SMF, DCF), e o amplificador vem depois da
última; (b) um campo `dcf: FiberConfig | None` no trecho; (c) `line_compensation = True` cria a DCF
com o comprimento que zera a dispersão acumulada ($D_1L_1 + D_2L_2 = 0$, Agrawal, seção 8.2.1).

**Perguntas.** Qual? E os parâmetros da DCF (D, atenuação, área efetiva): do Agrawal (seção 8.2) ou
do professor/VPI? Até a semana 5 os notebooks montam as duas fibras à mão, então isso só trava a
montagem pela configuração.

**Resposta:** (Gabriel, 2026-11-16) "a opção (a) se alinha mais com a ideia de blocos, mas eu não colocaria o EDFA obrigatoriamente. A função deveria ser mais flexível para que o span também seja flexível". Como está: o amplificador já é opcional (`SpanConfig.amplifier = None`), e fibra padrão + compensadora + amplificador se escreve como dois elementos de `link.spans` (a fibra padrão sem amplificador, a compensadora com). **Não mudei** `SpanConfig.fiber` para uma tupla `fibers`: mexe em `resolve`, nas duas montagens, nas regras e nos testes de configuração do grupo. **Falta decidir:** basta a forma de hoje, quer a tupla `fibers`, ou uma lista ordenada de elementos (fibras e amplificadores em qualquer ordem)?

### 53. BPSK e polarização única no coerente

**O que é.** A aula começa o coerente por BPSK (E10), mas `DCS_FORMATS` tem QPSK, 16QAM, 32QAM e
64QAM, e o padrão do `CoherentChannel` é dupla polarização. Acrescentar BPSK aos formatos? E fazer as
semanas 4 e 5 em uma polarização (`Pol.SP`), deixando a dupla para depois?

**Resposta:** (Gabriel, 2026-11-16) "pode ser, não entendi o que você quis dizer". Em outras palavras: a aula começa o coerente pelo BPSK, que não estava na lista de formatos, e com uma polarização. A pergunta era se eu podia acrescentar o BPSK e fazer os primeiros notebooks do coerente com uma polarização só. Foi o que fiz (notebook 08).

### 54. Modelo do receptor coerente

**O que é.** E12 e E13 precisam de oscilador local, híbrida de 90° e detecção balanceada.
Escolhas para a primeira versão:
- híbrida ideal (sem desbalanço de fase nem de perda) ou já com os desbalanços do VPI
  (`InsertionLossImbalanceI/Q` = 0 no RX, dúvida 18);
- LO derivado do próprio laser do TX (homódino ideal, sem ruído de fase relativo) ou um laser livre
  com `lo_frequency_offset` (500 MHz na configuração, D-20) e largura de linha, que exige
  recuperação de fase e frequência no DSP;
- detecção balanceada ideal (dois fotodetectores iguais, subtração perfeita) ou com responsividades
  diferentes.

**Minha recomendação:** tudo ideal e LO derivado do laser na semana 4; LO livre na semana 5.

**Resposta:** (Gabriel, 2026-11-16) "híbrida, LO e detecção devem ser flexíveis para ser tanto ideal quanto real, dependendo dos parâmetros do usuário, e isso vale para mais componentes [...] você pode criar um .md para discutirmos isso". Criado: `docs/proposta_ideal_real.md`.

### 55. Leitura das linhas E11 a E13 da aula

**O que é.** Minha leitura (`docs/planejamento.md`, seção 1): E11, IQ com fotodetector direto (a
fase se perde: com I e Q binários a potência é constante); E12, o "laser no final do MZM" é o
oscilador local somado antes de um fotodetector só (aparece o batimento, junto com $|E_s|^2$ e
$|E_{LO}|^2$); E13, detecção balanceada (sobra só o batimento, e a informação volta). Está certa?
E "se I e Q são OOK" quer dizer I e Q binários ±1 (QPSK), ou OOK mesmo (0 e 1, que não dá
potência constante)?

**Resposta:** (Gabriel, 2026-11-16) "I e Q binários". Isto é, ±1 em cada quadratura (QPSK). Já está assim.

---

## 11. Do coerente (`docs/sistema_coerente.md`, 2026-11-14)

Respostas de 2026-11-14 (Gabriel), já aplicadas: escopo da primeira rodada (passos da aula, BPSK e
QPSK, sem fibra), pulso de Nyquist (cosseno levantado, rolloff 0,18), só acrescentar código (os
blocos do IMDD ficam intactos) e sincronismo por correlação com os símbolos transmitidos, com a
sequência de sincronismo do professor documentada (`docs/sistema_coerente.md`, seção 4.1).

As dúvidas 10, 53, 54 e 55 foram respondidas em 2026-11-16. A 10 mudou o código (opção (a): a
validação exige `vbias = vpi`). A 54 virou a proposta de `docs/proposta_ideal_real.md`. A 53 e a 55
confirmaram o que estava feito.

### 56. Excursão do acionamento do modulador IQ (`tx.drive_depth`)

**O que é.** Com a pré-distorção, o campo de cada MZM é $m\,s$, com $s$ o sinal elétrico (±1 no
centro dos símbolos) e $m$ = `drive_depth`. O pulso de Nyquist chega a 2,1 entre os símbolos, então
com $m$ acima de 0,47 os picos são ceifados. O VPI tem `PreDistortion = Yes` e `Renormalize = On`,
mas o arquivo não traz a excursão. O professor usa `ind_mod = 0.1` V com `Vpi = 2.5` V, sem
pré-distorção.

**Como está.** Padrão 1,0 (o MZM inteiro) e 0,4 nos exemplos (sem ceifamento).

**Perguntas.** Qual valor de referência? E o padrão deve ser 1,0, um valor sem ceifamento, ou o
campo fica obrigatório?

**Resposta:** (Gabriel, 2026-11-16) "boa pergunta, veja no Agrawal". **[Fato: Agrawal, seção 10.1.2, texto depois da Eq. 10.1.5]** No QPSK "a tensão é alterada de modo que a função de transferência $t_m$ assuma os valores ±1": a excursão inteira do MZM, com sinal binário. O Agrawal não trata o pulso de Nyquist nem o 16QAM, em que o sinal passa do nível do símbolo entre os símbolos. Fica: padrão 1,0 (Agrawal), 0,4 nos exemplos com pulso de Nyquist (sem ceifamento) e 0,51 no cenário do VPI (calibrado, dúvida 64).

### 57. Ruído térmico do receptor balanceado

**O que é.** O professor soma ruído térmico em cada um dos quatro fotodiodos (`DPQAM_receiver`), e
o bloco `CoherentFrontEnd` faz o mesmo: a variância térmica total é $4\,i_T^2\,\Delta f$. Num
receptor com um amplificador por par balanceado seriam duas fontes, não quatro.

**Pergunta.** Manter como o professor?

**Resposta:** (Gabriel, 2026-11-16) "sim". Mantido: ruído térmico em cada um dos quatro fotodiodos.

### 58. BPSK pelo modulador IQ

**O que é.** A aula faz o BPSK com um laser e um MZM. Na montagem pela configuração o BPSK passa
pelo `IQModulator`, com o braço Q parado no nulo: com extinção infinita dá o mesmo campo, com 6 dB
a menos, que o `PowerControl` repõe. Com extinção finita o braço Q vaza portadora em quadratura. O
notebook 08 monta o MZM único à mão.

**Pergunta.** Basta assim, ou o BPSK deve ter um bloco de um MZM só?

**Resposta:** (Gabriel, 2026-11-16) "basta assim, mas acho que essa mescla de MZM e laser vai dar problema no futuro, deveriam ser separados". Mantido. Na montagem o laser e o modulador já são blocos separados; o que está junto é a configuração. Está na seção 4 de `docs/proposta_ideal_real.md`.

### 59. Perda de inserção do modulador IQ

**O que é.** `mzm.insertion_loss_db` é aplicada em cada MZM interno, como no IMDD. O interferômetro
externo já tira 3 dB por si. Não sei se o valor de 6 dB da configuração é do modulador IQ inteiro
ou de cada MZM.

**Resposta:** (Gabriel, 2026-11-16) "eu não sei também". Continua aberta. Hoje a perda é aplicada em cada MZM interno.

---

## 12. Da comparação com o VPI (`docs/comparacao_vpi.md`, 2026-11-14)

### 60. Total das tabelas de BER e SER do grupo

**O que é.** `resultado_SER_por_polarizacao_4canais_125km_canal2.csv` tem `Total_X = 524288`, mas
cada arquivo de constelação tem 262144 pontos. Contando os erros de símbolo a 13 dBm em X, achei
1487 em 262144 (5,7e-3); a tabela tem 1497 em 524288 (2,9e-3). A BER da tabela é a SER da tabela
dividida por 4. Os alvos A e B do pacote (`reference/vpi/ajuste_abc.py`) vieram dessas tabelas.

**Pergunta.** O total está dobrado, ou as tabelas vieram de outro conjunto de pontos?

**Resposta:** (Gabriel, 2026-11-16) "um erro de bit pode ser um erro de símbolo. Se 10 símbolos estão errados, pelo menos 10 bits tiveram que ser alterados para que isso ocorresse". Conferi nos dados (125 km, X): a tabela conta as componentes I e Q em separado. A 13 dBm ela tem 1497 erros em 524288; contei 1491 componentes erradas em 2 × 262144 e 1487 símbolos errados em 262144 (a 8 dBm: 8807, 8821 e 8717; a 17 dBm: 6049, 6031 e 5994). Então o total não está dobrado: a "SER" da tabela é a taxa de erro por componente, metade da taxa por símbolo. A BER é que fica baixa demais: pelo seu raciocínio, 1487 símbolos errados dão pelo menos 1487 bits errados em 4 × 262144, isto é, BER ≥ 1,42e-3, e a tabela tem 7,1e-4, a metade. A comparação do simcopy não usa essas tabelas.

### 61. `PowerdBm` do VPI é a soma dos 4 canais (D-22)

**O que é.** Com a potência lida como total, o termo não linear do simcopy fica a 5 % do VPI em 125
e em 140 km. Lida como potência por canal, ficaria 16 vezes maior. O YAML do 400ZR tem
`output_power_dbm: 13.0` no booster, que no simcopy é por canal: para reproduzir o VPI a 13 dBm é
preciso 13 − 6,02 dBm por canal (o script `compara_vpi.py` faz essa conta).

**Perguntas.** Confirma que é a soma? E o YAML deve passar a ter a potência por canal equivalente
(6,98 dBm), ou o booster ganha um campo de potência total?

**Resposta:** (Gabriel, 2026-11-16) "sim, é a soma". Confirmado. O YAML continua com 13 dBm no booster e o script `compara_vpi.py` faz a conta por canal.

### 62. Constelação de referência do VPI

**O que é.** `dados_ideais/dados_ref_canal2` tem SNR de 21,97 dB contra a grade ideal, sem imagem de
quadratura e sem deslocamento. Não é a saída do transmissor do cenário (que tem os 5° e a extinção
de 30 dB). Usei-a só para saber os símbolos transmitidos.

**Pergunta.** Como ela foi gerada?

**Resposta:** (Gabriel, 2026-11-16) "elas representam os dados ideais, não lembro como foi gerada, mas só use-a para ver os símbolos que seriam". É o uso que ela tem.

### 63. Filtro casado no cenário do VPI

**O que é.** A cadeia de DSP do documento do VPI não tem filtro casado separado (compensação de
dispersão, depois o MIMO de 7 coeficientes). O YAML do 400ZR ficou com `matched_filter: false`, e
assim o ruído do receptor fecha com o VPI em 0,1 dB. Com o filtro casado ligado o simcopy fica
1,85 dB melhor que o VPI nesse termo.

**Pergunta.** Manter o cenário fiel ao VPI (sem filtro casado), como está?

**Resposta:** (Gabriel, 2026-11-16) "pode ser". Mantido sem filtro casado.

### 56 (atualização). Excursão do acionamento

O YAML do 400ZR usava `drive_depth: 0.47`, calibrado para reproduzir o deslocamento constante das
constelações do VPI (0,077 do valor rms). Em 2026-11-16 passou a 0,51, junto com o sinal dos 5°
(dúvida 64). O padrão da configuração continua 1,0.

### 26 (evidência). Densidade do ruído térmico

Com a convenção unilateral e sem filtro casado, o termo de ruído do receptor fecha com o VPI em
0,1 dB (dúvida 63). É evidência a favor da convenção unilateral, não prova: depende de a leitura da
cadeia de DSP do VPI estar certa.

---

## 13. Da terceira rodada (`docs/sistema_coerente.md`, seção 12, 2026-11-16)

Implementado nesta rodada, com efeito sobre dúvidas antigas: a inclinação da dispersão (dúvida 35),
o ADC por amostras por símbolo (dúvida 40) e o PAM4 (dúvida 51). As três foram respondidas em
2026-11-16.

**Estado em 2026-11-16, depois das respostas do Gabriel.** Todas as dúvidas têm resposta. Ficam
com decisão pendente, marcadas com "Falta decidir", "Confirmar" ou "Em aberto" na resposta: 27
(padrão da banda do filtro do receptor e comparação com outros simuladores), 39 (os padrões que
quebram testes do grupo), 42 (o valor "52"), 52 (forma do trecho), 59 (perda de inserção do
modulador IQ), 66, 67 e 68 (para discutir), e a proposta de `docs/proposta_ideal_real.md` (54 e
58).

### 64. Sinal do desbalanço de fase de 5° do VPI

**O que é.** O VPI tem `PhaseShift = -90+5`, isto é, 85° entre os braços. O YAML do 400ZR tem
`phase_imbalance_deg: 5.0`, que no simcopy dá 95°. Os dados do VPI sem os 5° (120 km) mostram que o
deslocamento constante da constelação cai de 0,077 para 0,073 quando o desbalanço sai. Com +5° o
simcopy prevê o contrário (0,077 para 0,080). Com −5° e `drive_depth: 0.51`, bate nos dois
conjuntos (0,076 e 0,073), e a imagem de quadratura vai de 0,042 para 0,044 (VPI 0,044).

**Como está.** Não apliquei: o YAML continua com +5° e `drive_depth: 0.47`, que são os valores da
comparação já publicada em `docs/comparacao_vpi.md`. Trocar muda o piso dessa comparação, que
precisa rodar de novo.

**Pergunta.** Aplico −5° e 0,51, e refaço a comparação de 125 e 140 km?

**Resposta:** (Gabriel, 2026-11-16) "sim, inclusive você tem dados de outros km para testar". Aplicado: o YAML passou a `phase_imbalance_deg: -5.0` e `drive_depth: 0.51`, e a comparação foi refeita em 120, 125, 130, 135 e 140 km, nos dois modos de DSP, mais 120 km sem os 5°: `docs/comparacao_vpi.md`, seção 10.

### 65. Fundo de escala do ADC

**O que é.** O VPI tem `Range = 1.0` no ADC, sem dizer a escala da corrente. O bloco `ADC` usa
`adc.full_scale_rms` = 4,0: o fundo de escala é 4 vezes o valor rms de cada componente (ganho
ideal). Com 12 bits isso deixa o ruído de quantização 53 dB abaixo do sinal, então o valor não pesa
no cenário do VPI. Pesa em estudos com poucos bits.

**Pergunta.** Há um valor de referência? Ou fica 4,0?

**Resposta:** (Gabriel, 2026-11-16) "não sei". Fica 4,0.

### 66. Padrões e ordem do receptor cego

**O que é.** Três escolhas minhas, sem fonte direta:
- a separação no espaço de Stokes vem antes do equalizador, e no VPI vem depois;
- a etapa dirigida por raio (minha leitura do `MMA` do VPI) vem desligada, porque atrapalhou com
  pouca SNR, e quem refina é o ajuste por decisão;
- o `phase_window` padrão é 21, que eu supunha ser o `NPreSymbols` do VPI, e no modo cego deixa a
  fase saltar 90° com pouca SNR. (Correção de 2026-11-17: a janela do VPI é 43, isto é, 21 símbolos
  de cada lado. Medida nas constelações: `docs/replicacao_vpi_diario.md`.)

**Medido depois (`docs/comparacao_vpi.md`, seção 9).** No cenário de 4 canais, com 21 a fase
salta abaixo de cerca de 13 dB de SNR, que é onde o VPI também se perde (1 dB acima em potência).
Com 81 não salta de 10 a 16 dB, mas a SNR em alta potência cai até 0,5 dB, porque a janela longa
não acompanha a fase não linear. Com 21 o simcopy cego fica de 0,21 a 0,33 dB acima do VPI.

**Perguntas.** Mantenho a ordem e a etapa por raio desligada? O padrão de `phase_window` fica 21,
como no VPI, e quem quiser SNR baixa sem saltos escolhe 81?

**Resposta:** (Gabriel, 2026-11-16) "parece complexo, precisamos discutir". Em aberto. **Em palavras simples:** das três escolhas, só uma pede decisão. A ordem das duas primeiras etapas não mudou nenhum resultado, e a etapa que desliguei só atrapalhava. O que falta escolher é quantos símbolos o receptor usa para tirar a média da fase. Com 21 ele acompanha variações rápidas e ganha até 0,5 dB em potência alta, mas se perde quando há muito ruído. Com 81 não se perde, e fica até 0,5 dB pior em potência alta. **Atualização (2026-11-17):** a janela do VPI é 43, não 21, e com 43 o receptor cego do simcopy reproduz a não linearidade do VPI (`docs/replicacao_vpi_diario.md`). Agora eu deixaria 43 como padrão.

### 67. PAM4: níveis de potência e instante de decisão

**O que é.** (i) Os níveis elétricos do PAM4 são igualmente espaçados, e o MZM os comprime nas
pontas. O exemplo usa metade da excursão para os olhos ficarem quase iguais, com extinção de 7,7
dB. Não há pré-distorção dos níveis. (ii) A medida do PAM4 procura o instante de decisão na janela
inteira, e a do OOK continua procurando só de 0 a 2 símbolos (ponto 1 da revisão).

**Perguntas.** Entra a pré-distorção dos níveis? E posso levar a procura na janela inteira para o
OOK, o que corrige a medida em WDM sem DCF e pode mudar os números do notebook 07?

**Resposta:** (Gabriel, 2026-11-16) "parece complexo, precisamos discutir". Em aberto. **Em palavras simples:** são duas perguntas separadas. (i) No PAM4 os quatro níveis de luz não saem igualmente espaçados, porque o MZM comprime as pontas. Dá para corrigir os níveis elétricos antes do modulador. Eu poria como opção, desligada por padrão. (ii) A medida do OOK só procura o melhor instante de amostragem nos dois primeiros símbolos. Num canal WDM sem fibra compensadora o atraso é maior e a medida erra. A do PAM4 procura na janela inteira. Eu levaria isso para o OOK, e alguns números do notebook 07 podem mudar.

### 68. `dispersion_slope`: `null` contra `0`

**O que é.** Com `null` (padrão) a fibra não tem $\beta_3$, como antes. Com `0`, a dispersão $D$ é
constante no comprimento de onda, o que pela Eq. 2.3.13 do Agrawal ainda dá $\beta_3 \neq 0$. O
YAML do 400ZR passou a ter `0.08e3` (o valor do VPI). O efeito medido em 125 km e 13 dBm foi de
−0,02 dB na SNR.

**Pergunta.** O padrão deve continuar `null` (sem $\beta_3$), ou passar a `0`?

**Resposta:** (Gabriel, 2026-11-16) "parece complexo, precisamos discutir". Em aberto. **Em palavras simples:** a pergunta é se a fibra tem dispersão de terceira ordem mesmo quando o usuário não fala nada. Hoje não tem (`null`). Pela sua resposta 35 ("poderia entrar sempre"), o padrão seria um valor típico, por exemplo os 0,05 ps/(nm²·km) do Agrawal. O efeito é pequeno (0,02 dB no cenário do VPI), mas muda a terceira casa dos números dos notebooks com fibra.

### 69. Mapeamento do 64QAM

**O que é.** O `draft.py` só vai até 32 pontos. O 64QAM usa a regra do 16QAM do professor com três
bits por eixo e código de Gray refletido.

**Pergunta.** Serve, ou há uma tabela de referência?

**Resposta:** (Gabriel, 2026-11-16) "código Gray não pode ser generalizado para 32QAM, mas para 64QAM sim". Mantido o mapeamento do 64QAM. O 32QAM continua fora.

### 70. Gráfico de símbolos enviados e recebidos

**O que é.** Anotado em 2026-11-18, a partir do notebook 08. Hoje o `plot_constellation`
(`plotting.py`) mostra a nuvem dos símbolos recebidos com os pontos ideais, e o histograma de cada
quadratura. Ele não mostra, símbolo a símbolo, o que foi enviado contra o que chegou, então um erro
aparece como um ponto na nuvem, mas não dá para ver qual símbolo errou nem em que instante. O
Gabriel quer esse gráfico no notebook 08 e, como o `plot_constellation` e o `plot_chain` são do
pacote, a mudança fica aqui para ajustarmos depois.

**Proposta.** Uma função `plot_symbols(result, channel=0, n_symbols=64)` em `plotting.py`, com dois
painéis (em fase e quadratura): no eixo x, o índice do símbolo. Os símbolos enviados como cruzes e
os recebidos (já divididos pelo ganho complexo medido) como pontos, com um círculo nos que a decisão
errou. Os dados já existem: os enviados vêm de `map_bits` (`kernels/demodulation.py`) e os recebidos
de `sample_symbols` no instante de decisão, como no `plot_constellation`.

**Pergunta.** Serve uma função nova, ou é melhor um argumento no `plot_constellation`? Enquanto
isso, o gráfico fica montado à mão no notebook 08.

**Resposta:**

### 71. Notebook 08: símbolos sem ADC, SNR teórica e controle de potência

**O que é.** Perguntas do Gabriel sobre o notebook 08, anotadas em 2026-11-19 para responder
depois (ele pediu que não fossem respondidas na hora):
1. Como os símbolos são tirados no receptor se o sistema do notebook não tem ADC?
2. Por que a comparação com a SNR teórica, se é um sistema óptico com outras fontes de ruído e não
   um sistema ideal?
3. Essa análise não fica limitada pelo controle de potência?

**Resposta:** (pendente, do lado do Claude)
