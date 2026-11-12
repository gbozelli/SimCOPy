# Moduladores: guia de implementação passo a passo

Documento de trabalho para implementar os moduladores externos no simcopy, seguindo a
ordem da aula (`docs/04-class/context.md`). Marcação igual à de `docs/laser.md`:

- **[Fato: fonte]**: está escrito na fonte indicada (equação, seção ou linha de código).
- **[Dedução]**: obtida a partir de equações citadas; a conta está no texto.
- **[Proposta]**: escolha de implementação; precisa de aprovação.
- **[Dúvida]**: não está nas fontes; precisa de resposta do Gabriel ou do orientador.
- **[Decidido: Gabriel, data]**: resposta do Gabriel a uma dúvida.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014): seção 3.4.2
   (modulação externa), seção 10.1.2 (moduladores de amplitude e de fase, Eq. 10.1.2–10.1.5,
   Fig. 10.3), Ex. 10.4, Eq. 4.7.1 (razão de extinção).
2. Código do professor, `legacy/draft.py`: `MZM` (linha 943), `DPQAM_transmitter`
   (linha 1235), `QAM_transmitter` (linha 327).
3. Notas da aula: `docs/04-class/context.md`.
4. Projeto: `reference/vpi/parametros_vtmu.md` (módulo `Tx_mQAM_DP` do VPI) e `docs/roadmap.md`
   (no pacote `simcopy_1`, ainda fora do repo), `config/system.py` (`MZMConfig`, `IQImpairments`).

---

## 1. O que a aula pede

**[Fato: context.md]** "Moduladores: podem ser descritos por monocromático, pois se usam
ondas contínuas (geralmente não há modulação no laser, exceto para sistemas específicos)."

**[Fato: context.md]** Na cadeia IMDD: "Laser de onda contínua – MZM – fibra – fotodetector".
No coerente: "inicialmente BPSK: Laser – MZM (modulação em amplitude, fase = 0)", depois
"Laser – MZM (modulação QAM, I(t) e Q(t))".

**[Dedução]** A entrada óptica do modulador é o laser do modelo B (`CWLaser`, já
implementado no passo 8 do laser): potência fixa e ruído de fase. O modulador não precisa
das equações de taxa.

---

## 2. Física (Agrawal)

**[Fato: Agrawal, seção 3.4.2]** Acima de ~5 Gb/s o chirp da modulação direta é grande
demais; o laser fica em corrente constante (CW) e um modulador externo codifica os dados.
Dois tipos: LiNbO₃ em Mach-Zehnder (MZ) e eletroabsorção (MEA, razão de extinção de 15 dB
ou mais, algum chirp). A aula usa o MZ.

**[Fato: Agrawal, Eq. 10.1.2 e texto seguinte]** No braço j de um MZ de LiNbO₃ a fase é
$\phi_j(t) = \pi V_j(t)/V_\pi$, com $V_\pi$ tipicamente de 3 a 5 V.

**[Fato: Agrawal, Eq. 10.1.3–10.1.4]**

$$A_t = \tfrac{1}{2}A_i\left(e^{i\phi_1} + e^{i\phi_2}\right),\qquad
t_m = \frac{A_t}{A_i} = \cos\!\left[\tfrac{1}{2}(\phi_1 - \phi_2)\right]
\exp\!\left[\tfrac{i}{2}(\phi_1 + \phi_2)\right]$$

O MZ modula amplitude e fase. **[Fato: Agrawal, texto após a Eq. 10.1.4]** Com
$V_2 = -V_1 + V_b$ (push-pull) a soma das fases é constante e o MZ é um modulador de
amplitude puro (sem chirp):

$$T_m(t) = |t_m|^2 = \cos^2\!\left[\frac{\pi}{2V_\pi}\,(2V_1(t) - V_b)\right]\qquad\text{(Eq. 10.1.5)}$$

Com a mesma tensão nos dois braços ($\phi_1 = \phi_2$) ele é um modulador de fase puro.

**[Fato: Agrawal, Fig. 10.3 e texto]** Modulador em quadratura: dois MZ internos em
push-pull, um em cada braço de um MZ externo, com a tensão $V_3$ dando $\pi/2$ entre os
braços. Para QPSK a saída tem as fases π/4, 3π/4, 5π/4 e 7π/4, $(\pm1 \pm i)/\sqrt{2}$.

**[Fato: Agrawal, Eq. 4.7.1]** Razão de extinção: $r_{ex} = P_0/P_1$ (potência do bit 0 sobre
a do bit 1). O livro não modela o MZ com extinção finita, perda de inserção ou banda
elétrica.

**[Fato: Agrawal, Ex. 10.4]** Exercício: deduzir a função de transferência em função de
$V_1$ e $V_2$ e as condições para modulação de amplitude pura.

---

## 3. Código do professor (`legacy/draft.py`)

**[Fato: linhas 943–954]** `MZM(A, Vrf_upper, Vrf_lower, Splitter_IL, Upper_IL, Lower_IL,
Combiner_IL, Vpi, Vbias)`:
- cada braço recebe $A/2$ vezes a perda do divisor ($10^{-IL/20}$, IL em dB);
- fase do braço de cima $\pi(-V_{rf,up} + V_{bias}/2)/V_\pi$; do de baixo
  $\pi(-V_{rf,low} - V_{bias}/2)/V_\pi$;
- perdas próprias de cada braço (padrão 2 dB e 3 dB) e do combinador; a saída é a soma.

**[Fato: linhas 1235–1326]** `DPQAM_transmitter`:
- laser com `X_fraction = 0,5` (metade da potência em cada polarização);
- em cada polarização, portadora/√2 em dois MZM; o braço Q recebe a portadora vezes
  $e^{i\pi/2}$; a saída é $(E_I + E_Q)/\sqrt{2}$;
- push-pull: $V_{rf,up} = K_{mod}\,s$, $V_{rf,low} = -K_{mod}\,s$, com `ind_mod` = $K_{mod}$ = 0,1;
- $V_\pi = V_{bias} = 2{,}5$ V; perdas 1 / 2 / 2 / 1,1 dB.

**[Fato: linha 349]** Os símbolos são normalizados pelo máximo da parte real (=1) antes da
conversão D/A.

**[Dedução]** Com os braços iguais, o MZM do professor é a Eq. 10.1.3 com polarização
simétrica ($\pm V_{bias}/2$): $\phi_1 - \phi_2 = \pi(V_{bias} - 2V_{rf})/V_\pi$ e
$\phi_1 + \phi_2 = 0$. Mesma transmissão da Eq. 10.1.5; difere da forma do Agrawal só pela
fase constante $\pi V_b/(2V_\pi)$.

**[Dedução]** Braços com amplitudes $a_1 \ne a_2$ (perdas diferentes) não zeram a saída:
$T_{max} = (a_1 + a_2)^2$, $T_{min} = (a_1 - a_2)^2$, razão de extinção
$[(a_1 + a_2)/(a_1 - a_2)]^2$. Com as perdas padrão do professor (2 e 3 dB), **24,8 dB**.
É assim que o código do professor produz extinção finita.

**[Fato: notebook do professor `3_LaserAndModulator`, cópia em `pacote/`]** Mesma função
`MZM` (perda do braço de baixo de 2 dB no padrão). A varredura de polarização imprime
$T = 0{,}5$ para $V_{bias} = 2{,}5$ com $V_\pi = 5$ (quadratura). Exemplos com um braço só,
push-push e push-pull (este com $\pm V/2$ em cada braço, diferença $V$), tom harmônico com
aumento da excursão (harmônicos) e exercícios 3 a 5. O laser desse notebook usa
`phase = 2*pi*(-Freq_offset)*t` (linha 224), com o sinal oposto ao do `draft.py` e do
notebook 6 (`+Freq_offset`): relevante para M1. O `DPQAM_transmitter` usa $\pm K_{mod}s$ por
braço (diferença $2K_{mod}s$): relevante para M3.

**[Dedução]** No `DPQAM_transmitter`, $V_{bias} = V_\pi$: polarização no nulo, campo
$\propto \sin(\pi V_{rf}/V_\pi)$, bipolar. A fase de pico por braço é
$\pi \cdot 0{,}1/2{,}5 = 0{,}126$ rad: regime quase linear (o seno desvia 0,3 % da reta), com
potência de saída baixa.

---

## 4. VPI (`reference/vpi/parametros_vtmu.md`, módulo `Tx_mQAM_DP`)

**[Fato]** `Vpi = 5` V, `InsertionLoss = 6` dB, `ExtinctionRatio = 30` dB,
`S21_Bandwidth = 40e9` Hz, `S21_FilterOrder = 3`, `PhaseShiftX/Y = -90+5` graus ("π/2,
add IQ phase imbalance here"), `SkewX/Y = 0.75e-12` s ("adds delay of I w.r.t. Q driver
signal"), `SplitterExtinctionRatio` e `CombinerExtinctionRatio` = 30 dB, `PreDistortion = Yes`,
`Renormalize = On`, `Gain = 0` dB, laser 1 mW e 500 kHz. **Não há parâmetro de tensão de
polarização** no módulo.

**[Fato: config/system.py]** `MZMConfig(vpi=5.0, vbias=2.5, insertion_loss_db=6.0,
extinction_ratio_db=30.0, s21_bandwidth=40e9, s21_order=3)`, usado por `ExternalModulation`
(IMDD) e por `CoherentTransmitter`; `IQImpairments(phase_imbalance_deg=0, skew=0)`.

**[Dedução]** Na convenção do Agrawal e do professor, `vbias = 2,5` com `vpi = 5` é
polarização em **quadratura** ($T = 1/2$ em $V = 0$): serve ao OOK, mas não ao QAM, que
precisa do **nulo** ($V_b = V_\pi$). O padrão atual do `CoherentTransmitter` ficaria na
quadratura (ver dúvida M2).

---

## 5. Implementado (2026-10-22): kernel `kernels/modulator/mzm.py`

Só física do Agrawal, sem configuração:
- `mzm_field(campo, v1, v2, vpi, a1=1/2, a2=1/2)`: Eq. 10.1.3, com amplitudes de braço
  opcionais (o desbalanço do professor).
- `push_pull(v, v_bias)`: tensões $(v - V_b/2,\,-v + V_b/2)$, a forma simétrica do professor.
- `iq_field(campo, v_i, v_q, vpi, v_bias, quadrature_phase, a1, a2)`: Fig. 10.3,
  $\tfrac{1}{2}[t(v_i) + e^{i\theta}t(v_q)]$. **[Proposta]** θ é parâmetro obrigatório (π/2 no
  professor, −π/2 no VPI: ver dúvida M1).
- Testes: com o grupo (issue #15; **[Decidido: Gabriel, 2026-10-22]** aqui só se implementa).
  Os 12 testes e as 7 mutações que eu tinha escrito estão em `rascunhos/testes_claude/`.
- Notebook `notebooks/03_moduladores.py`: seções 0–6 (ver o próprio notebook).

---

## 6. Proposta de blocos (aguarda as dúvidas)

**[Proposta]** Dois componentes em `components/tx/modulator.py`, com a interface de todo
bloco, `__call__(x, ctx)`, entrada elétrica e saída óptica (2, N):
- `ExternalModulator(laser: CWLaserConfig, mzm: MZMConfig, frequency, ...)` (IMDD): `CWLaser`
  + MZM em push-pull, campo em X.
- `IQModulator(laser, mzm, iq: IQImpairments, frequency, pol)` (coerente): `CWLaser` + IQ;
  em DP, potência do laser dividida entre X e Y, um IQ por polarização (como o professor).

A arquitetura (`docs/03-class/arquitetura-simulador-optico (1).md`, Anexo A.4) já prevê
`ExternalModulator(cfg.tx.laser, cfg.tx.mzm)` e `IQModulator(cfg.tx.laser, cfg.tx.mzm,
n_pol=...)` no `build_tx`.

---

## 7. Dúvidas

**M1. Convenção de frequência (vale para o projeto inteiro).** **[Fato: Agrawal, Eq. 2.4.3,
2.4.5 e 2.4.12]** O campo é $A\,e^{-i\omega_0 t}$ e $\delta\omega = -\partial\phi/\partial t$; a
fibra do pacote segue essa convenção (Eq. 2.6.18). **[Dedução]** Então, no eixo da FFT do
NumPy, +f é fisicamente o lado vermelho. Hoje o laser CW aplica o desvio como
$e^{+i2\pi\Delta f t}$ e o `DirectModulatedLaser` devolve $e^{+i\phi}$ com
δν = +(1/2π)dφ/dt = azul (decisão de 2026-10-15); na convenção da fibra, as duas coisas saem
com o sinal trocado (o chirp da modulação direta chegaria à fibra invertido). No IQ, o VPI usa
−90° e o professor +90°. Opções: (A) manter a convenção do Agrawal na fibra e ajustar laser,
IQ e mux a ela; (B) adotar a convenção de engenharia ($e^{+i\omega_0 t}$, +f = azul), o que
troca os sinais da NLS no kernel ($\partial A/\partial z = -\tfrac{\alpha}{2}A
+ \tfrac{i\beta_2}{2}\partial_T^2 A - i\gamma|A|^2A$) e a definição de C da gaussiana. O
`docs/arquitetura.md` (seção 14.3) já lista isso como ponto crítico de compatibilidade com o
VPI.

**M2. O que é `vbias`?** Tensão de polarização no push-pull simétrico (como o professor),
com quadratura em $V_\pi/2$ e nulo em $V_\pi$? Se sim, o coerente precisa de `vbias = vpi`
(o padrão atual, 2,5 com 5, é quadratura). Ou a polarização deve sair do formato
(nulo para QAM, quadratura para OOK) e não ser parâmetro, como no VPI?

**M3. Escala do sinal elétrico.** Como o `Waveform` elétrico vira tensão? O professor usa
$V = K_{mod}\cdot s$ com $K_{mod}$ = 0,1 V e $s$ normalizado pelo máximo; o `MZMConfig` não
tem esse parâmetro. **[Proposta]** Um campo `drive_amplitude` em múltiplos de $V_\pi$ (como o
laser de modulação direta usa múltiplos de $I_{th}$).

**M4. Razão de extinção.** Modelar pelo desbalanço dos braços, como o professor? **[Dedução]**
Com perda de inserção IL e extinção ER: $a_1 + a_2 = 10^{-IL/20}$ e
$(a_1 - a_2)/(a_1 + a_2) = 10^{-ER/20}$. Isso traz um pequeno chirp residual (o campo deixa de
ser real), que é físico.

**M5. Filtro S21 (40 GHz, ordem 3).** Que tipo (Bessel, Butterworth, gaussiano)? Aplicado à
tensão de cada braço antes do MZM?

**M6. `PreDistortion = Yes` do VPI.** **[Dedução, não verificada no VPI]** Provavelmente
pré-distorce o sinal elétrico com arco-seno para compensar o seno do MZM e deixar o campo
linear nos símbolos. Implementar? (O roadmap já pedia para confirmar.)

**M7. Skew de 0,75 ps.** Atraso fracionário do sinal I aplicado por fase linear na FFT. Isso
conflita com a decisão D-4 (sem reamostragem fracionária)? É um filtro, não reamostragem.

**M8. Dupla polarização.** Dividir a potência do laser meio a meio entre X e Y, como o
professor? E o que são os `SplitterExtinctionRatio` / `CombinerExtinctionRatio` de 30 dB do
VPI: extinção de polarização do divisor e do combinador de polarização, ou dos acopladores Y
do IQ?

**M9. Potência de saída.** O VPI tem `Renormalize = On`. A saída do modulador deve ser
laser × perdas × perda de modulação (o que o roadmap pede como aceitação: "potência média da
saída bate com laser × perdas"), ou renormalizada?

---

## 8. Decisões de 2026-10-22 (`docs/duvidas_abertas.md`)

- **M1 [Decidido: Gabriel]** Convenção do Agrawal. O sinal de θ no IQ só precisa combinar com o
  receptor coerente.
- **M4 [Decidido: Gabriel]** Extinção finita pelo desbalanço dos braços (professor), com
  $a_1 + a_2 = 10^{-IL/20}$ e $(a_1 - a_2)/(a_1 + a_2) = 10^{-ER/20}$.
- **M5 [Decidido: Gabriel]** S21 com Bessel de ordem 3; o VPI assume banda infinita (sem banda
  configurada, o bloco não filtra).
- **M6 [Decidido: Gabriel]** Pré-distorção: o VPI tem. Aplicado: `arcsin_predistortion` no kernel
  (no nulo, $V = (V_\pi/\pi)\arcsin(m\,s)$ dá campo $m\,s$); painel (d) da seção 6 do notebook 03.
- **M7 [Decidido: Gabriel]** Skew por fase linear na FFT.
- **M8 [Decidido: Gabriel, parcial]** Estrutura do VPI: laser dividido meio a meio entre X e Y.
  **[Aberto]** Os 30 dB de `Splitter/CombinerExtinctionRatio`: vazamento entre X e Y ou acopladores
  do IQ?
- **M9 [Decidido: Gabriel]** Renormalizar a saída. **[Minha leitura, a confirmar]** para a potência
  do laser (`AveragePower` do VPI).
- **Desbalanço de ganho I/Q [Decidido: Gabriel]** Entra, com padrão 0. Aplicado:
  `IQImpairments.gain_imbalance_db = 0.0` e `q_gain` no `iq_field`.
