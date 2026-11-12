# Planejamento: 5 semanas para o sistema (TX, enlace e RX)

Plano para entregar a parte da aula que monta o sistema de comunicação a partir dos componentes
(`docs/04-class/context.md`, linhas 27–47). Marcação igual à dos guias: **[Fato: fonte]**,
**[Dedução]**, **[Proposta]**, **[Dúvida]**, **[Decidido: Gabriel, data]**.

Escrito em 2026-10-25. As dúvidas citadas por número estão em `docs/duvidas_abertas.md`.

---

## 1. O que a aula pede

**[Fato: context.md, linhas 27–47]** A sequência da aula, numerada aqui para o cronograma citar:

| Etapa | Aula | Cadeia |
|---|---|---|
| E1 | "Criar sequência pseudoaleatória, converter para sinal elétrico (não se usa Nyquist, e sim um filtro de Bessel)" | bits → sinal elétrico |
| E2 | "Colocar o laser e ver a saída óptica, entrada elétrica" | laser |
| E3 | "Colocar laser e fibra diretos, sem o fotodetector" | laser – fibra |
| E4 | "depois laser, fibra e fotodetector, e monitorar a saída óptica em todos" | laser – fibra – fotodetector |
| E5 | "Laser de onda contínua – MZM – fibra – fotodetector" | CW – MZM – fibra – PD |
| E6 | "Colocar duas fibras ao invés de uma (SMF, DCF)" | ... – SMF – DCF – ... |
| E7 | "Colocar um amplificador depois da fibra DCF" | ... – SMF – DCF – EDFA – ... |
| E8 | "Múltiplos spans" | trechos repetidos |
| E9 | "Múltiplos canais (mux e demux, enlace compartilhado)" | WDM |
| E10 | "Coerente digital, inicialmente BPSK"; "Laser – MZM (modulação em amplitude, phase = 0)" | CW – MZM (BPSK) |
| E11 | "Laser – MZM (modulação QAM, I(t) e Q(t)); se I e Q são OOK, é constante, não há recuperação de fase – photodetector" | CW – IQ – PD |
| E12 | "... mas no final do MZM tem um laser – photodetector" | CW – IQ – (+ LO) – PD |
| E13 | "... mas no final do MZM tem um laser, balanceador – photodetector (aí vai conseguir recuperar info)" | CW – IQ – (+ LO) – detecção balanceada |

**[Fato: context.md]** A aula também registra que os passivos "não se pode testar sozinhos, pois
precisa da cadeia de transmissão TX/RX", e a linha "IMDD < 5w", de sentido ainda em aberto
(dúvida 42).

**[Dedução, a confirmar: dúvida 55]** Leitura de E11 a E13: em E11 o fotodetector direto só vê
$|E|^2$, e com I e Q binários a potência é constante, então a informação (que está na fase) se
perde. Em E12 o oscilador local (o "laser no final do MZM") é somado ao sinal antes do
fotodetector, e o batimento traz a fase para a corrente, junto com os termos $|E_s|^2$ e
$|E_{LO}|^2$. Em E13 a detecção balanceada cancela esses termos e sobra só o batimento: a
informação é recuperada.

---

## 2. Premissas

- **[Premissa, a confirmar]** 1 h por dia útil e 2 h em cada dia do fim de semana: 9 h por semana,
  45 h no total. Se forem 2 h no fim de semana inteiro, são 35 h; ver a seção 6.
- **[Proposta]** Início na segunda, 26/10/2026, e fim no domingo, 29/11/2026.
- Testes unitários ficam fora deste plano: são do grupo, pelas issues (**[Decidido: Gabriel,
  2026-10-22]**). Cada semana termina com o notebook e a issue de revisão dele, como as #12 a #21.
- Os notebooks montam as cadeias à mão, bloco a bloco. A montagem a partir da configuração
  (`build_tx`, `build_link`, `build_rx`, seção 4 da arquitetura) fica para a semana 5.
- Cada etapa precisa de um critério de pronto com forma fechada (Agrawal), como já fazem os
  notebooks 01 a 06.

---

## 3. O que existe e o que falta

**[Fato: código em `src/simcopy`, 2026-10-25]**

| Componente | Kernel | Bloco | Notebook |
|---|---|---|---|
| Fibra (SSFM, PMD) | `kernels/ssfm.py`, `birefringence.py` | falta | 01 |
| Laser (CW, equações de taxa) | `kernels/laser/` | `CWLaser`, `DirectModulatedLaser` | 02 |
| MZM e IQ | `kernels/modulator/mzm.py` | falta | 03 |
| EDFA | `kernels/amplifier.py` | falta | 04 |
| Fotodetector | `kernels/photodetector.py` | falta | 05 |
| Mux e demux | `kernels/wdm.py` | falta | 06 |

Falta por inteiro:
- **Elétrico (TX):** sequência pseudoaleatória (E1), sinal NRZ com filtro de Bessel (E1), mapeador
  de símbolos e formatador de pulso para o coerente (E10).
- **Recepção:** decisão e contagem de BER; oscilador local, híbrida de 90° e detecção balanceada
  (E12, E13); DSP mínimo (recuperação de fase, compensação de dispersão).
- **Monitor:** o `Probe` da arquitetura (seção 5.5), que grava a saída de qualquer bloco sem que o
  bloco saiba. É o que a aula pede em E4 ("monitorar a saída óptica em todos").
- **Montagem:** `build_tx`, `build_link`, `build_rx` e o executor (`pipelines/`, vazio).

**[Fato: `config/system.py`]** Pontos da configuração que não comportam a aula como está (viram
dúvidas):
- `IMDDTransmitter.rolloff = 0.2`: formatação de Nyquist; a aula pede Bessel (dúvida 50).
- O formato padrão do IMDD é `PAM4` (dúvida 51).
- `SpanConfig` tem uma fibra só e `line_compensation: bool`; a aula pede SMF + DCF + EDFA (dúvida 52).
- `DCS_FORMATS` não tem BPSK, e o padrão do coerente é dupla polarização (dúvida 53).

---

## 4. Cronograma

Em cada semana: as dúvidas que precisam de resposta antes de começar, as tarefas com horas, a
entrega e o critério de pronto.

### Semana 1 (26/10 a 01/11): primeira cadeia IMDD, laser direto (E1 a E4)

Antes: dúvidas 2, 4, 6, 9, 25 a 28, 33, 34, 42, 47, 49 e 50.

| Tarefa | Horas |
|---|---|
| Responder as dúvidas da semana | 1,5 |
| Merge do PR `laser` → `main` (dúvida 47); `f0` dos lasers vira o desvio (dúvida 2) | 0,5 |
| Blocos de sequência pseudoaleatória e de sinal NRZ com filtro de Bessel (E1) | 2 |
| Blocos de fibra e de fotodetector | 2 |
| `Probe` para gravar as saídas | 1 |
| Notebook 07: E2 (saída óptica × entrada elétrica), E3 (+ fibra), E4 (+ fotodetector), com as saídas monitoradas | 2 |

Pronto quando: potência média na saída da fibra igual a $P\,e^{-\alpha L}$; corrente média igual a
$R_d\bar P$ mais a corrente de escuro; o alargamento do pulso com chirp na fibra concorda com o
notebook 01. **[Dedução]** É também o primeiro teste da convenção de frequência entre laser e
fibra: com o chirp do laser direto e dispersão anômala, o pulso tem que alargar.

### Semana 2 (02/11 a 08/11): modulação externa, BER e compensação de dispersão (E5, E6)

Antes: dúvidas 3, 10, 11, 17, 51 e 52.

| Tarefa | Horas |
|---|---|
| Responder as dúvidas da semana | 0,5 |
| Bloco do modulador externo (CW + MZM na quadratura; extinção, perda, S21 Bessel de ordem 3) | 2 |
| Bloco de decisão e contagem de BER (limiar ótimo, Agrawal Eq. 4.6.6–4.6.11) | 1,5 |
| Notebook 08: E5; laser direto × externo na mesma fibra; BER × potência recebida e × distância | 3 |
| E6: SMF + DCF montadas à mão no notebook 08 | 1,5 |
| Folga | 0,5 |

Pronto quando: BER back-to-back sobre ½erfc(Q/√2) (como no notebook 05); com
$D_1L_1 + D_2L_2 = 0$ (Agrawal, seção 8.2.1) o sinal na saída da DCF volta ao back-to-back,
a menos da não linearidade.

### Semana 3 (09/11 a 15/11): amplificador, vários trechos e WDM (E7 a E9, fim do IMDD)

Antes: dúvidas 21, 22, 23, 24, 29, 30 e 31.

| Tarefa | Horas |
|---|---|
| Responder as dúvidas da semana | 0,5 |
| Bloco do amplificador (ganho fixo ou potência por canal) | 1,5 |
| Notebook 09: E7 e E8; OSNR e BER × número de trechos | 3 |
| Blocos de mux e demux | 1,5 |
| Notebook 10: E9, canais IMDD num enlace compartilhado | 2 |
| Folga | 0,5 |

Pronto quando: OSNR depois de $N_A$ trechos igual à Eq. 7.4.1–7.4.2 do Agrawal (o notebook 04 já
faz isso só com amplificadores); no WDM, um canal sozinho e o mesmo canal com vizinhos dão a mesma
BER quando a interferência prevista pela Eq. 6.3.1 é desprezível.

### Semana 4 (16/11 a 22/11): coerente, BPSK e QPSK (E10 a E13)

Antes: dúvidas 10, 15, 40, 53, 54 e 55.

| Tarefa | Horas |
|---|---|
| Responder as dúvidas da semana | 0,5 |
| Mapeador de símbolos (BPSK, QPSK, 16QAM) e formatador de pulso | 1,5 |
| Bloco do modulador IQ, uma polarização (pré-distorção, desbalanços, skew) | 1,5 |
| Notebook 11, parte 1: E10 e E11 (fotodetector direto não recupera a fase) | 2 |
| Oscilador local, híbrida de 90° e detecção balanceada | 2 |
| Notebook 11, parte 2: E12 e E13; BER × SNR | 1,5 |

Pronto quando: em E11, a potência detectada é constante para QPSK; em E13, a constelação volta e a
BER segue a forma fechada do receptor coerente (Agrawal, seção 10.3); em E12, os termos
$|E_s|^2$ e $|E_{LO}|^2$ aparecem na corrente e somem em E13. Fontes: Agrawal, seção 4.5
(detecção coerente: 4.5.1, oscilador local; 4.5.2, detecção homódina), 10.2 e 10.3.

### Semana 5 (23/11 a 29/11): 16QAM, coerente com fibra e montagem pela configuração

Antes: dúvidas 5, 35, 37, 39, 41 e 43.

| Tarefa | Horas |
|---|---|
| Notebook 12: 16QAM back-to-back, BER × OSNR | 1,5 |
| Recuperação simples de fase e de frequência (oscilador local livre; Agrawal, seção 10.4.2) | 2 |
| Notebook 13: coerente com fibra e compensação de dispersão no receptor (Agrawal, seção 8.7.3) | 2 |
| `build_tx`, `build_link`, `build_rx` a partir da configuração (YAML) | 2 |
| Revisão final, merge e índice dos notebooks | 1,5 |

Pronto quando: a BER do 16QAM segue a forma fechada; uma configuração YAML reproduz o notebook 08
sem montar nada à mão.

---

## 5. Resumo das entregas

| Semana | Blocos novos | Notebook | Etapas |
|---|---|---|---|
| 1 | sequência, NRZ + Bessel, fibra, fotodetector, `Probe` | 07 | E1–E4 |
| 2 | modulador externo, decisão e BER | 08 | E5, E6 |
| 3 | amplificador, mux, demux | 09, 10 | E7–E9 |
| 4 | mapeador, formatador, IQ, LO, híbrida, balanceado | 11 | E10–E13 |
| 5 | recuperação de fase, compensação no RX, montagem | 12, 13 | coerente completo |

---

## 6. O que cortar se faltar tempo

Nesta ordem:
1. Montagem pela configuração (semana 5): os notebooks seguem montando as cadeias à mão.
2. Recuperação de fase com oscilador local livre: usa-se o LO derivado do próprio laser (homódino
   ideal).
3. Coerente com fibra (notebook 13).

Com 35 h em vez de 45 h, os três saem e a semana 5 vira revisão e folga.

Ficam fora do plano desde já: dupla polarização (dúvida 16) e WDM coerente.

---

## 7. Riscos

- **Dúvidas em aberto.** Cada semana começa com as dúvidas que a travam; se a resposta atrasar, o
  bloco atrasa junto. As da semana 1 são as mais urgentes.
- **Coerente.** É a parte com mais escolhas novas (híbrida, LO, fase) e sem kernel ainda; por isso
  fica com as duas últimas semanas e com os cortes acima.
- **Tempo de execução.** O laser de equações de taxa roda em Python puro, a passos de 0,1 ps; com
  milhares de bits, uma execução pode levar dezenas de segundos. Os notebooks usam sequências
  curtas, e a contagem de BER com muitos bits fica com a modulação externa.
