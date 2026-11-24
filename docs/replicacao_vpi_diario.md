# Replicação do VPI: diário

Registro do laço de `docs/plano_replicar_vpi.md`: o placar, e uma entrada por medida ou hipótese,
com a previsão escrita antes do teste e o veredito depois.

---

## Placar de partida (2026-11-16)

Cenário de `examples/400zr_4ch_75ghz.yaml` (−5°, excursão de 0,51), $2^{14}$ símbolos, semente 0.
Gerado por `python scripts/compara_vpi.py placar`. O critério é o DSP cego.

**DSP cego, janela de fase de 21:**

| Grandeza (simcopy menos VPI) | 120 km | 125 km | 130 km | 135 km | 140 km | Alvo | Estado |
|---|---|---|---|---|---|---|---|
| imagem [%] | −3 a 0 | −2 a −1 | −3 a +1 | −4 a 0 | −3 a −1 | ±2 % | fora |
| deslocamento [%] | −1 a 0 | −1 | −1 a 0 | −1 | −1 a 0 | ±2 % | dentro |
| $A$ [dB] | −0,13 | −0,31 | −0,20 | −0,27 | −0,19 | ±0,1 dB | fora |
| $B$ [dB] | −0,38 | −0,31 | −0,34 | −0,32 | −0,40 | ±0,1 dB | fora |
| resto do piso [dB] | −0,71 | −0,60 | −0,71 | −0,61 | −0,87 | ±0,15 dB | fora |
| $C$ [dB] | −0,56 | −0,58 | −0,62 | −0,59 | −0,54 | ±0,15 dB | fora |
| SNR ponto a ponto [dB] | +0,19 a +0,42 | +0,31 a +0,39 | +0,24 a +0,44 | +0,28 a +0,42 | +0,27 a +0,43 | ±0,1 dB | fora |
| potência ótima [dBm] | 13 e 13 | 14 e 14 | 14 e 14 | 14 e 14 | 15 e 15 | igual | dentro |
| SER sobre a gaussiana [%] | −13 a +6 | −10 a +6 | −3 a +6 | −1 a +5 | −3 a +12 | ±10 % | fora |

**DSP com símbolos conhecidos** (controle): $A$ de −0,05 a +0,08 dB (dentro), $B$ de −0,25 a −0,37
dB, resto do piso de −0,50 a −0,75 dB, $C$ de −0,29 a −0,35 dB, SNR de +0,07 a +0,30 dB.

---

## Etapa A (2026-11-16): medidas nas constelações do VPI

Medidas novas de `scripts/medidas_constelacao.py`, com o mesmo código nos dados do VPI (125, 130,
135 e 140 km, 6 a 18 dBm, $2^{18}$ símbolos) e nos símbolos do simcopy (125 e 140 km, os dois modos
de DSP). O script confere cada medida em dados sintéticos de resposta conhecida
(`python scripts/medidas_constelacao.py sintetico`). Resultados em
`reference/vpi/replicacao/etapa_a_vpi.json` e `etapa_a_simcopy.json`.

### O que a etapa A achou

1. **A janela da média de fase do VPI é de 43 símbolos, não de 21.** Sai dos dados, em todos os
   pontos convergidos dos quatro comprimentos. `NPreSymbols = 21` é meia janela.
2. **O piso que falta no simcopy é imagem de quadratura que depende da frequência, com a forma
   exata de um skew entre I e Q.** No VPI ela vale 0,00233 e no simcopy 0,00113, em todos os pontos.
   A diferença, 0,0012, é o que faltava no resto do piso (0,0009). No modelo do simcopy o valor do
   VPI corresponde a um skew de 1,07 ps, e o arquivo diz 0,75 ps. Uma amostra da grade do VPI são
   1,044 ps.
3. **A diferença na não linearidade está toda no ruído de fase.** Metade do termo $C$ do VPI é
   ruído de fase que a janela de 43 não acompanha. A parte que não é fase bate com o simcopy em 5 a
   13 %.
4. **As quedas do VPI são falhas da recuperação de fase, e a assinatura é uma constelação girada de
   45°, não um salto de 90°.** O equalizador e a separação de polarização continuam funcionando
   nesses arquivos.
5. **Quatro hipóteses caem sem simular nada:** a sequência periódica (D1), a diafonia entre
   polarizações (B2), o processamento em blocos e a assimetria entre X e Y. A da potência de saída
   do transmissor (B1) deixa de ser necessária.

### A0. A sequência de símbolos é periódica?

**Base.** **[Fato: VPI]** `PRBS_Type = PRBS`, `PRBS_Order = 7`. `docs/arquitetura.md` já pedia
para conferir.

**Medida.** Nos símbolos transmitidos do canal 2 (constelação de referência): não há período exato
até 5000 símbolos, a fração de símbolos iguais no atraso 127 é 6,3 % (o acaso dá 6,25 %), a
autocorrelação fora do zero não passa de 0,006, e X e Y não são a mesma sequência deslocada.

**Veredito.** A sequência não se repete. **D1 refutada** para o canal sob teste. Não tenho os
símbolos dos vizinhos.

### A1 e A2. Como o VPI se perde

**Medida.** Giro de cada bloco de 128 símbolos em relação aos dados, sem arredondar, e correlação
do bloco com os dados.

| L [km] | P [dBm] | Fração do bloco girada mais de 22,5° | Trechos | Menor correlação com os dados |
|---|---|---|---|---|
| 125 | 6 | 0,36 | 3 | 0,88 |
| 125 | 18 | 0,52 | 2 | 0,93 |
| 130 | 6 e 7 | 0,33 e 0,11 | 3 e 2 | 0,87 e 0,90 |
| 135 | 8 e 18 | 0,35 e 0,83 | 2 e 2 | 0,73 e 0,91 |
| 140 | 9 e 10 | 0,41 e 0,11 | muitos e 1 | 0,90 e 0,90 |
| 140 | 18 | 0,51 | muitos | 0,90 |
| todos os pontos convergidos | | 0,00 | 0 | 0,97 ou mais |

- Nos arquivos perto do limiar a correlação de cada bloco com os dados fica acima de 0,87, e nem
  nos piores (140 km a 6 dBm) cai abaixo de 0,66. O que se perde é a fase: a constelação fica
  girada por trechos longos, de milhares a centenas de milhares de símbolos.
- O giro fica perto de ±45° (125 km a 6 dBm: entre −45° e +48°; 140 km a 10 dBm: termina em −45°),
  e em alguns arquivos chega a −90° e fica (125 e 135 km a 18 dBm).
- Em 135 km a 6 dBm a fase dá uma volta inteira ao longo do bloco: aí a recuperação de fase está
  solta.

**Veredito.** As duas quedas são da recuperação de fase. Não são ruído, nem divergência do
equalizador, nem perda de sincronismo. O contador de saltos de 90° que eu usava não serve para
estes arquivos: o estado girado de 45° fica na fronteira do arredondamento. A linha 10 do placar
passa a ser a fração do bloco girada.

**No simcopy** o receptor cego não girou em nenhum destes pontos (125 km a 6 e 18 dBm, 140 km a 9,
10 e 18 dBm). Quando ele falhou, na seção 9.2 de `docs/comparacao_vpi.md`, foi com saltos de 90°.
Um estado estável em 45° é um modo de falha que a busca de fase do simcopy não tem.

### A3. Quanto do resto é fase, amplitude ou aditivo

**Medida.** Variância radial e tangencial do erro em cada anel do 16QAM. A parte tangencial que
cresce com o raio é ruído de fase.

Em 125 km, em unidades de $1/\mathrm{SNR}$:

| P [dBm] | Fonte | Resto | Aditivo | Fase | Amplitude |
|---|---|---|---|---|---|
| 8 | VPI | 0,03068 | 0,03044 | 0,00027 | −0,00004 |
| 8 | simcopy cego | 0,02804 | 0,02775 | 0,00012 | 0,00016 |
| 13 | VPI | 0,01693 | 0,01567 | 0,00124 | −0,00001 |
| 13 | simcopy cego | 0,01523 | 0,01443 | 0,00067 | 0,00009 |
| 16 | VPI | 0,02014 | 0,01530 | 0,00454 | 0,00013 |
| 16 | simcopy cego | 0,01783 | 0,01472 | 0,00301 | −0,00003 |
| 17 | VPI | 0,02498 | 0,01721 | 0,00723 | 0,00027 |
| 17 | simcopy cego | 0,02196 | 0,01691 | 0,00489 | −0,00004 |

A parte de fase cresce com o quadrado da potência. Ajustando fase = $f_0 + C_f P^2$:

| | $C_f$ em 125 km | $C_f$ em 140 km | $C$ total em 125 km | $C$ menos $C_f$ |
|---|---|---|---|---|
| VPI | 2,80e-6 | 2,83e-6 | 5,70e-6 | 2,90e-6 |
| simcopy cego, janela de 21 | 1,94e-6 | 1,72e-6 | 4,99e-6 | 3,05e-6 |
| simcopy, símbolos conhecidos | 2,10e-6 | 2,05e-6 | 5,29e-6 | 3,19e-6 |

**Veredito.** Metade da não linearidade do VPI aparece como ruído de fase. A parte que não é fase
é 5 % maior no simcopy em 125 km (13 % em 140 km). A diferença em $C$ (linha 6 do placar) está na
parte de fase, isto é, em quanto a recuperação de fase acompanha. Ruído de amplitude: quase não há
(no máximo 0,0003, no VPI a 17 dBm).

### A4. Quanto do resto é interferência entre símbolos

**Medida.** Regressão do símbolo recebido nos 5 vizinhos de cada lado, lineares e conjugados, da
mesma polarização e da outra.

| | VPI | simcopy cego | simcopy, símbolos conhecidos |
|---|---|---|---|
| imagem dos vizinhos (núcleo conjugado), 125 km | 0,00233 | 0,00113 | 0,00113 |
| imagem dos vizinhos, 140 km | 0,00233 | 0,00115 | 0,00114 |
| coeficiente conjugado do vizinho ±1 | 0,0307 | 0,0211 | 0,0211 |
| interferência linear, 125 km, de 8 a 16 dBm | 0,00123 a 0,00045 | 0,00146 a 0,00065 | 0,00116 a 0,00033 |
| diafonia da outra polarização | 0,00010 | 0,00017 | 0,00017 |

- **Imagem dos vizinhos.** Nos dois simuladores o núcleo conjugado é real e antissimétrico (+ no
  vizinho −1, − no +1, e alternando), sem parte simétrica: é a derivada do pulso, a assinatura de
  um atraso entre I e Q. Não muda com a potência nem com o comprimento. No simcopy o coeficiente
  bate com a teoria para 0,75 ps (0,0218). No VPI é 1,45 vez maior.
- **Calibração**, com o transmissor do simcopy sem fibra e sem ruído:

  | Skew no simcopy | Coeficiente conjugado ±1 | Imagem dos vizinhos |
  |---|---|---|
  | 0,75 ps (o do arquivo) | 0,0215 | 0,00113 |
  | 1,044 ps (uma amostra a 16 por símbolo) | 0,0300 | 0,00219 |
  | 1,061 ps (0,75 vezes a raiz de 2) | 0,0305 | 0,00226 |
  | 1,087 ps | 0,0312 | 0,00237 |
  | VPI, medido | 0,0307 | 0,00233 |

  O VPI corresponde a 1,07 ps. **[Inferência]** O VPI aplica o atraso em amostras inteiras, e 0,75
  ps vira uma amostra, 1,044 ps. Uma amostra dá 0,0300, 2 % abaixo do medido, e a dispersão da
  medida é de 1 %: a explicação é provável, não fechada.
- **Interferência linear.** É do mesmo tamanho nos dois, e cai quando a potência sobe nos dois: é
  o equalizador de 7 coeficientes trocando interferência por ruído. Não é erro em excesso do VPI.
- **Diafonia.** Desprezível nos dois.

**Veredito.** O que falta no resto do piso do simcopy (0,0009) é a imagem dos vizinhos (0,0012 a
mais no VPI). **B2 refutada** (a extinção dos divisores de polarização não aparece). **B1 deixa de
ser necessária**: o piso fecha sem mexer na potência do transmissor. **B3 confirmada na medida**, e
falta o teste no enlace.

### A5. Memória do erro de fase, e a janela

**Medida.** Covariância do erro de fase entre X e Y em função do atraso. Para uma média móvel de
$W$ símbolos ela é negativa até o atraso $(W-1)/2$ e dá um degrau positivo logo depois.

- No VPI o degrau está entre os atrasos 21 e 22 em todos os 38 pontos convergidos, entre X e Y e
  dentro de cada polarização: $W = 43$. Em 125 km a 7 dBm, a covariância vale −21e-5 no atraso 21
  e +15e-5 no 22.
- No simcopy a mesma medida lê 21, que é o valor configurado, em 40 dos 42 pontos. Nos outros dois
  (receptor cego a 6 e 7 dBm) o degrau não passa do ruído com $2^{14}$ símbolos. A medida funciona.
- O ruído de fase não linear tem memória curta: no VPI a 17 dBm a covariância cai de 737e-5 no
  atraso 0 para 203e-5 no atraso 3. Uma janela de 43 símbolos não tem como acompanhá-lo.

**Veredito.** **C1 confirmada na medida:** a janela do VPI é 43. Falta o teste no simcopy.

### A6. A SER é a de ruído gaussiano?

SER medida sobre a SER gaussiana para a mesma SNR, em 125 km: no VPI 0,99 a 8 dBm, 0,84 a 13 dBm e
1,13 a 17 dBm. No simcopy cego 0,96, 0,83 e 1,15. Com símbolos conhecidos 0,90, 0,70 e 1,03.

**Veredito.** O receptor cego do simcopy tem a mesma forma de erro do VPI. A curtose do erro também
bate (+0,59 e +0,50 a 17 dBm). A linha 9 do placar depende de mais símbolos no simcopy para sair
da barra de erro, não de uma hipótese.

### A7 e A8. Ao longo do bloco, e X contra Y

O erro do VPI é uniforme ao longo do bloco (desvio de 1 % entre 32 trechos, o esperado só pela
estatística), e X e Y têm a mesma SNR (diferença de 0,01 dB).

**Veredito.** Não há transiente do equalizador nem processamento em blocos visível, e não há
assimetria de polarização. As duas hipóteses saem.

### Outras duas medidas

- **Rotação que depende do raio.** O anel externo aparece girado em relação ao interno, e o giro
  cresce com a potência: −0,022 rad no VPI e +0,020 rad no simcopy, em 125 km a 13 dBm. O módulo é
  9 % maior no VPI. **[Inferência]** O sinal oposto é a convenção da portadora: o simcopy usa a do
  Agrawal, $\exp(-i\omega_0 t)$, e as constelações do VPI são o conjugado. Vale para o mapa de
  convenções de `docs/plano_validacao_cruzada.md`.
- **Ruído do receptor.** Tirando a interferência, a parte aditiva do VPI é cerca de 5 % maior que
  a do simcopy cego em potência baixa. É a linha $A$ do placar (−0,13 a −0,31 dB) e continua sem
  explicação.

### Estado das hipóteses depois da etapa A

| # | Hipótese | Estado |
|---|---|---|
| B1 | Potência de saída do transmissor 1 dB menor | desnecessária: o piso fecha por B3 |
| B2 | Extinção do divisor e do combinador de polarização | refutada: diafonia de 0,0001 |
| B3 | Skew do VPI diferente do meu | confirmada na medida: 1,07 ps contra 0,75 ps. Falta o teste no enlace |
| B4 | `Renormalize` fixa a potência do transmissor | desnecessária, pelo mesmo motivo de B1 |
| C1 | A janela de fase do VPI é 43 | confirmada na medida. Falta o teste no simcopy |
| C2 | A média móvel do VPI não é a minha | aberta, depende do resultado de C1 |
| C3 | O equalizador do VPI deixa erro em excesso | enfraquecida: a interferência linear é igual nos dois. Sobra como candidata aos 5 % de ruído do receptor |
| C4 | A ordem das etapas | aberta, prioridade baixa |
| C5 | A queda vem de outra etapa que não a busca de fase | reformulada: a queda é da fase, e o modo de falha é um giro de 45° |
| C6 | `ConstellationAlignment` altera a constelação | sem sinal disso: ganho, imagem e deslocamento são constantes |
| D1 | Sequência periódica | refutada |
| D2 a D5 | Enlace | abertas. A parte da não linearidade que não é fase já bate em 5 a 13 % |

### Previsões para o próximo passo (escritas antes de rodar)

**Teste 1: skew de uma amostra (1,044 ps) no transmissor.**
Previsão: a imagem dos vizinhos vai de 0,00113 para 0,0022. O resto do piso sobe 0,0011 nos cinco
comprimentos e fica em 0,0075 (VPI: 0,0073), isto é, de −0,6 dB para +0,1 dB. A imagem de mesma
frequência, o deslocamento, $A$ e $C$ não mudam.

**Teste 2: janela de fase de 43 no receptor cego.**
Previsão: $C_f$ sobe de 1,9e-6 para perto de 2,8e-6, e $C$ de 5,0e-6 para 5,8e-6 (VPI: 5,7e-6),
isto é, de −0,6 dB para +0,1 dB. $A$ muda pouco. A SNR a 17 dBm cai cerca de 0,4 dB.

**Os dois juntos.** Em 125 km a SNR do receptor cego vai de 16,44 para cerca de 16,13 dB a 14 dBm
(VPI: 16,09). A 8 dBm vai de 14,46 para cerca de 14,33 dB (VPI: 14,14): sobram 0,2 dB, que são os 5
% de ruído do receptor ainda sem explicação.

**Depois:** achar que variante da busca de fase fica presa em 45° (a linha 10 do placar), e os 5 %
de ruído do receptor.

### O que preciso de você para seguir

1. Posso rodar os testes 1 e 2? O teste 1 muda o skew do cenário. Eu rodaria primeiro numa cópia,
   e só mudaria `examples/400zr_4ch_75ghz.yaml` se a previsão se confirmar, marcando o valor como
   medido nos dados.
2. O grupo consegue confirmar no VPI como o `SkewX` é aplicado (em amostras inteiras?) e o que é
   `NPreSymbols` no `BlindPhaseSearch` (símbolos de cada lado)?
3. As perguntas 2 e 3 da seção 8 do plano (banco de sinais e código só de experimento) continuam
   abertas. A busca do giro de 45° precisa das duas.

---

## Testes 1 e 2 (2026-11-17): skew de uma amostra e janela de fase de 43

Aprovados pelo Gabriel em 2026-11-16 ("pode seguir").

**Como.** 120, 125, 130, 135 e 140 km com o skew de uma amostra (1,044 ps), e 125 e 140 km com o
skew antigo (0,75 ps) como controle. Em cada ponto, quatro cadeias de DSP sobre o mesmo sinal
propagado: símbolos conhecidos com janela de 21 e de 43, e receptor cego com janela de 21 e de 43.
$2^{14}$ símbolos, semente 0. Resultados em `reference/vpi/replicacao/teste12_*.json`, figura em
`reference/vpi/replicacao/teste12.png`. O sinal de cada ponto na saída do ADC ficou guardado em
`resultados/replica_vpi/` (ignorada pelo git), para testar outros DSP sem propagar de novo.

### Previsões e resultados

| Previsão (escrita antes) | Resultado | Veredito |
|---|---|---|
| Teste 1: imagem dos vizinhos de 0,00113 para 0,0022 | 0,00222 a 0,00224 (VPI: 0,00233) | confirmada |
| Teste 1: resto do piso fica em 0,0075 (VPI: 0,0073) | 0,00725 a 0,00752, receptor cego com janela de 21 | confirmada |
| Teste 1: imagem de mesma frequência, deslocamento, $A$ e $C$ não mudam | não mudaram (diferenças de até 0,05 dB) | confirmada |
| Teste 2: $C_f$ de 1,9e-6 para perto de 2,8e-6 | 2,64e-6 e 2,63e-6 (VPI: 2,80e-6 e 2,83e-6) | confirmada |
| Teste 2: $C$ de 5,0e-6 para 5,8e-6 (VPI: 5,7e-6) | 5,75e-6 a 5,95e-6 | confirmada |
| Teste 2: SNR a 17 dBm cai cerca de 0,4 dB | 0,31 dB em 125 km, 0,27 dB em 140 km | confirmada, um pouco menor |
| Os dois: 125 km a 14 dBm vai para cerca de 16,13 dB (VPI: 16,09) | 16,11 dB | confirmada |
| Os dois: 125 km a 8 dBm vai para cerca de 14,33 dB (VPI: 14,14) | 14,28 dB | confirmada |

**Os dois efeitos são independentes.** O controle com 0,75 ps mostra que a janela de 43 sozinha leva
$C$ de −0,58 para +0,05 dB e quase não mexe no piso. O skew sozinho leva o resto do piso de −0,60
para +0,10 dB e não mexe em $C$.

**SNR ponto a ponto**, nos 37 pontos em que o VPI converge:

| Cadeia | Simcopy menos VPI | Média |
|---|---|---|
| Cego, janela de 21, skew antigo (antes) | +0,19 a +0,44 dB | +0,34 dB |
| Cego, janela de 21, skew de uma amostra | +0,02 a +0,29 dB | +0,16 dB |
| **Cego, janela de 43, skew de uma amostra** | **−0,04 a +0,14 dB** | **+0,04 dB** |
| Símbolos conhecidos, janela de 21, skew de uma amostra | −0,06 a +0,14 dB | +0,03 dB |

Alguns pontos, em dB:

| L [km] | P [dBm] | VPI | Cego 21, 0,75 ps (antes) | Cego 43, 1 amostra |
|---|---|---|---|---|
| 125 | 8 | 14,14 | 14,46 | 14,28 |
| 125 | 14 | 16,09 | 16,44 | 16,11 |
| 125 | 17 | 14,84 | 15,27 | 14,83 |
| 135 | 14 | 15,51 | | 15,50 |
| 140 | 11 | 14,07 | 14,34 | 14,20 |
| 140 | 15 | 15,20 | 15,55 | 15,22 |
| 140 | 17 | 14,45 | 14,88 | 14,46 |

**Com símbolos conhecidos e janela de 43**, $C$ passa 0,2 a 0,5 dB do VPI. Esse modo tira a fase
antes do equalizador, polarização por polarização, e a janela longa não lhe serve. Por isso o YAML
continua com 21 nesse modo.

### Placar depois dos testes 1 e 2

A configuração que replica o VPI: skew de uma amostra, receptor cego, janela de 43.

| Grandeza (simcopy menos VPI) | 120 km | 125 km | 130 km | 135 km | 140 km | Alvo | Estado |
|---|---|---|---|---|---|---|---|
| imagem [%] | −2 a +1 | −2 a 0 | −2 a 0 | −3 a −1 | −3 a −1 | ±2 % | quase |
| deslocamento [%] | −1 a +1 | −1 a 0 | −1 a +1 | −1 a 0 | −1 a 0 | ±2 % | dentro |
| $A$ [dB] | −0,11 | −0,28 | −0,16 | −0,31 | −0,30 | ±0,1 dB | **fora** |
| $B$ [dB] | +0,01 | +0,07 | 0,00 | +0,16 | +0,14 | ±0,1 dB | 3 de 5 |
| resto do piso [dB] | +0,08 | +0,16 | +0,05 | +0,33 | +0,29 | ±0,15 dB | 2 de 5 |
| $C$ [dB] | +0,14 | +0,04 | +0,06 | −0,06 | −0,06 | ±0,15 dB | dentro |
| SNR ponto a ponto [dB] | −0,04 a +0,09 | −0,02 a +0,14 | −0,01 a +0,13 | −0,03 a +0,11 | +0,01 a +0,12 | ±0,1 dB | quase |
| potência ótima [dBm] | 13 e 13 | 14 e 14 | 14 e 14 | 14 e 14 | 15 e 14 | igual | 4 de 5 |
| SER sobre a gaussiana [%] | −5 a +8 | −7 a +3 | −5 a +2 | −7 a +2 | −4 a +7 | ±10 % | dentro |

- **O que sobra é uma diferença só: o ruído do receptor ($A$).** O receptor cego do simcopy tem
  menos ruído que o do VPI em potência baixa: 4,7 % a menos em 125 km a 8 dBm, 1 % a 13 dBm, e
  nenhum de 15 dBm em diante. Como $A$, $B$ e $C$ saem do mesmo ajuste, o $A$ baixo empurra $B$ e o resto do
  piso para cima em 135 e 140 km. Ponto a ponto, em potência alta, o resto do piso bate com o do
  VPI em 0,0002.
- **Potência ótima em 140 km:** o VPI tem 15,16 dB a 14 dBm e 15,20 dB a 15 dBm. A diferença de
  0,04 dB está dentro da barra de erro, e o alvo "igual" é estrito demais para um topo tão plano.

### As duas perguntas ao grupo, investigadas nos dados

**`NPreSymbols`: respondida pelos dados.** A janela é de 43 símbolos, isto é, 21 de cada lado, em
todos os pontos convergidos de todos os conjuntos: 125 a 140 km, 120 km com e sem os 5°, e também a
pasta de 80 km, que é outro cenário. Com 43 o simcopy reproduz $C$ (teste 2). A confirmação do grupo
só fecharia a pergunta formalmente.

**`SkewX`: a explicação mais simples funciona, sem prova.** Nos dados, o skew equivalente no modelo
do simcopy é de 1,05 a 1,08 ps (o coeficiente conjugado dos vizinhos vale 0,0301 a 0,0308 em todos
os conjuntos de 120 a 140 km). Uma amostra da grade (1,044 ps) dá 0,0297: explica 97 % do
coeficiente, e o que sobra vale 0,0001 em $1/\mathrm{SNR}$, menos de 0,1 dB. Outras peças do
transmissor do VPI que eu não modelo podem responder pelos 3 % que faltam, ou até por parte do
resto: `FilterActive = Advanced` com um passa-baixas de cosseno levantado, interpolação cúbica, e uma
reamostragem do sinal óptico a 600 GHz (`Resample`, seção 11.1 do documento). Só o VPI diz como o
`SkewX` é aplicado. Não bloqueia: está abaixo da tolerância.

### O giro de 45° (linha 10): o que testei

**Fatos dos dados do VPI.** As transições são abruptas (dentro de um bloco de 64 símbolos), os
estados ficam em múltiplos exatos de 45° (desvio de 0,7 a 1,9° nos trechos longos) e duram de 2 mil
a 170 mil símbolos.
X e Y giram juntos, e o giro vai e volta nos dois sentidos.

| Hipótese | Teste | Veredito |
|---|---|---|
| A busca de fase do 16QAM tem um mínimo em 45° | Custo da busca em função do ângulo, com ruído gaussiano de 12 a 30 dB | refutada: 45° é o máximo do custo em todas as SNR. Os mínimos secundários ficam em ±27° a ±35°, e somem abaixo de 13,5 dB |
| A decisão é numa grade de 256QAM (o VPI tem `Format = ReadFromFile`, `256QAM_constellation.txt`) | O mesmo custo com a grade de 256QAM, em 11 escalas do sinal | refutada: em nenhuma escala 45° vira mínimo local |
| O desdobramento da fase usa período de 45° em vez de 90° | Aplicado à última busca de fase do simcopy, sem realimentar o equalizador, nos 51 sinais guardados | inconclusiva: produz exatamente estados em múltiplos de 45°, mas a estimativa bruta do simcopy não salta mais de 22,5° no interior de nenhum bloco (os 57 saltos maiores que isso caem nos 60 símbolos de cada borda; no interior o maior é de 19,7°). Só reproduziria as taxas do VPI se a estimativa bruta do VPI fosse muito mais ruidosa |

Como o mínimo do custo nunca está em 45°, o giro tem de vir da contabilidade da fase (desdobramento
ou resolução da ambiguidade de 90°), e não da busca em si. O mecanismo exato fica em aberto.

**Verificação que o grupo pode fazer no VPI:** rodar um ponto em que ele se perde (125 km a 6 dBm)
com o `Format` da busca de fase em `16QAM` em vez de `ReadFromFile`, e, se der, exportar a fase
estimada pela busca nesse ponto.

### A linha $A$: o que testei

Sobre o sinal guardado de 125 km a 10 dBm, outras formas de equalizador cego:

| Variante | SNR [dB] | Aditivo | Interferência |
|---|---|---|---|
| VPI | 15,17 | 0,02217 | 0,00313 |
| simcopy (módulo constante e 4 ajustes por decisão) | 15,25 | 0,02176 | 0,00334 |
| um ajuste por decisão | 15,23 | 0,02201 | 0,00348 |
| módulo constante e dirigido por raio, sem ajuste | 14,98 | 0,02314 | 0,00381 |
| só módulo constante, 20 passagens | 14,83 | 0,02534 | 0,00507 |
| só módulo constante, 4 passagens | 14,64 | 0,02683 | 0,00522 |

O VPI tem menos interferência e mais ruído que o simcopy. Todas as variantes sem o ajuste por
decisão têm mais das duas coisas. Nenhuma reproduz o padrão.

**Candidatas, não testadas:** o equalizador do VPI se adapta por gradiente, símbolo a símbolo, e
pode parar longe da solução de mínimo erro; o filtro elétrico do transmissor (`FilterActive =
Advanced`), que eu não modelo; a reamostragem a 600 GHz. A primeira pede um parâmetro calibrado (o
passo de adaptação); as outras duas pedem saber o que o VPI faz.

### Estado das hipóteses

| # | Hipótese | Estado |
|---|---|---|
| B3 | Skew do VPI maior que 0,75 ps | **confirmada** (teste 1). Aplicada no YAML: uma amostra, 1,044 ps |
| C1 | A janela de fase do VPI é 43 | **confirmada** (teste 2). Fica no modo cego; o YAML continua no modo com símbolos conhecidos |
| C2 | A média móvel do VPI não é a minha | desnecessária: com 43, $C$ fecha |
| C3 | O equalizador do VPI deixa erro em excesso | candidata para $A$, não testada |
| C5 | A queda vem de outra etapa | o giro de 45° vem da contabilidade da fase; mecanismo em aberto |
| D2 a D5 | Enlace | desnecessárias para $C$: fecha em ±0,15 dB |

### O que mudou no repositório

- `examples/400zr_4ch_75ghz.yaml`: skew de 0,75 ps para 1,044 ps, com o comentário de onde saiu.
- Correções de texto onde eu dizia que a janela do VPI era 21: `config/system.py`,
  `docs/sistema_coerente.md`, `docs/comparacao_vpi.md`, `docs/duvidas_abertas.md`, README.
- `scripts/compara_vpi.py simular`: `--second-dsp` pode ser repetido (várias cadeias de DSP no mesmo
  sinal) e `--save-signals` guarda o sinal na saída do ADC.

### O que preciso de você

1. **Passo o YAML do cenário para o receptor cego com janela de 43?** É a configuração que replica o
   VPI (média de +0,04 dB). Muda os números do notebook 09 quando ele for reexecutado.
2. **Linha $A$:** implemento, como experimento fora de `src/`, um equalizador adaptado por gradiente
   como o do VPI, sabendo que o passo vira um parâmetro calibrado? Ou esperamos saber do grupo o que
   o VPI faz?
3. **Giro de 45°:** a verificação no VPI acima é a forma mais barata de avançar. Sem ela, o próximo
   passo seria a etapa E do plano (estatística da queda com 16 sementes), que leva de 2 a 3 horas de
   máquina e só mede, sem explicar.

---

## Respostas de 2026-11-17 e o que veio delas

Respostas do Gabriel às três perguntas acima:

1. **"Sim":** o cenário (`examples/400zr_4ch_75ghz.yaml`) passou ao receptor cego com janela de 43
   (`mode: blind`, `phase_window: 43`). Os 73 testes e as camadas passam. A célula do cenário no
   notebook 09 roda com ele (125 km, 13 dBm, $2^{13}$ símbolos: SNR de 16,02 dB). O notebook não foi
   reexecutado nem editado: os números mudam quando ele for, e o texto dele ainda diz que o DSP usa
   os símbolos conhecidos. Para refazer a seção 10 de `docs/comparacao_vpi.md`, que usava o modo
   antigo, a nota no fim daquela seção diz o que passar ao script.
2. **"Pode fazer":** o equalizador adaptado por gradiente está em
   `scripts/experimentos/equalizador_gradiente.py`. Resultado abaixo.
3. **"Vou pedir pro grupo":** as perguntas ao grupo (o giro de 45°, a configuração do VPI, os dados
   exportados e as dúvidas de modelagem ainda abertas) estão num documento escrito para quem conhece
   óptica e sistemas e não conhece o código:
   <https://claude.ai/artifact/64UkTUsgYcH4wetqK5SaA8>. São 24 perguntas, com uma tabela no fim
   para responder item a item.

### Linha $A$: equalizador adaptado por gradiente (hipótese C3)

**Como.** Sobre os sinais guardados na saída do ADC (skew de uma amostra, $2^{14}$ símbolos), em 125
km de 6 a 17 dBm e em 140 km a 11, 13, 15 e 17 dBm. Cadeia: espaço de Stokes, módulo constante (2
passagens), critério dirigido por raio (o `MMA` do VPI) por 20 passagens com o passo varrido de 1e-4
a 1e-2, e a busca de fase de 43. Uma segunda forma acrescenta 4 passagens dirigidas por decisão,
também por gradiente (LMS), no lugar do ajuste em bloco do simcopy. Cada forma sai de dois jeitos:
com os coeficientes finais aplicados ao bloco inteiro, ou com os coeficientes de cada instante
(como um equalizador que se adapta enquanto processa). Resultados em
`resultados/replica_vpi/gradiente/` (ignorada pelo git).

$A$, $B$ e $C$ ajustados nos mesmos pontos do VPI, em dB do simcopy sobre o VPI, saída com os
coeficientes de cada instante:

| Variante | L [km] | SNR menos VPI [dB] | $A$ | $B$ | $C$ |
|---|---|---|---|---|---|
| receptor cego do simcopy (ajuste em bloco) | 125 | −0,01 a +0,14 | −0,28 | +0,07 | +0,03 |
| só por raio, passo 3e-4 | 125 | −0,40 a −0,14 | −0,17 | +0,55 | +0,38 |
| só por raio, passo 1e-3 | 125 | −0,29 a −0,19 | +0,16 | +0,23 | +0,45 |
| por raio e LMS, passo 1e-3 | 125 | −0,05 a +0,10 | −0,23 | +0,09 | +0,10 |
| por raio e LMS, passo 3e-3 | 125 | −0,09 a +0,02 | −0,13 | +0,13 | +0,11 |
| receptor cego do simcopy (ajuste em bloco) | 140 | +0,01 a +0,13 | −0,36 | +0,21 | −0,12 |
| por raio e LMS, passo 3e-3 | 140 | −0,07 a +0,01 | −0,17 | +0,22 | 0,00 |

- **Só por raio**, com qualquer passo, o equalizador fica de 0,1 a 0,8 dB abaixo do VPI em todas
  as potências, com mais interferência e mais ruído (a 8 dBm em 125 km, interferência de 0,0043 a
  0,0058 contra 0,0036 do VPI). O equalizador do VPI chega a uma solução melhor que essa.
- **Com a etapa LMS**, o passo pequeno reproduz o ajuste em bloco. Um passo maior soma ruído de
  desajuste em todas as potências por igual: com 3e-3 a potência baixa fecha (+0,02 dB a 8 dBm),
  mas a alta abre (−0,09 dB a 16 dBm). O excesso do VPI some de 15 dBm em diante, e o desajuste é
  uma fração do erro em qualquer potência. Não reproduz.
- As duas formas de saída diferem pouco com passo de até 1e-3.

**Veredito.** C3 refutada nas duas formas. O passo de 3e-3 com LMS deixa a SNR ponto a ponto dentro
de ±0,1 dB, mas troca o erro de $A$ por erro em $B$, e não entra no simulador.

### Linha $A$: ruído da busca de fase

**Base.** O excesso do VPI some em potência alta, que é a forma de uma diferença só no termo $A/P$.
O ruído de estimação da fase vem do ruído aditivo e escala com $1/P$. O bloco `BlindPhaseSearch`
do VPI tem `JointEstimation = No` e decide na grade de `256QAM_constellation.txt`, e as duas coisas
aumentariam esse ruído.

**Medida.** Nas medidas do erro, a 8 e a 13 dBm em 125 km a parte de fase é igual nos dois
(simcopy 0,00027 e 0,00126, VPI 0,00027 e 0,00124). O excesso está na parte aditiva: 0,02903
contra 0,03044 a 8 dBm (+4,9 %) e 0,01549 contra 0,01567 a 13 dBm (+1,2 %).

**Teste** (`scripts/experimentos/busca_fase_vpi.py`, receptor cego do simcopy com a busca trocada):

| Busca de fase | Pontos que convergem | Quedas |
|---|---|---|
| por polarização | iguais (±0,03 dB) | nenhuma |
| grade de 256QAM, 16QAM sobre os pontos internos (escala 1) | iguais | nenhuma |
| grade de 256QAM, escalas de 1,41 a 3 | a busca não funciona em nenhum ponto | |
| grade de 256QAM, 16QAM com a potência do 256QAM (escala 4,12), conjunta | ruído de fase 5 a 20 vezes o do VPI, SNR 0,5 a 1 dB abaixo | nenhuma |
| o mesmo, por polarização | ruído de fase ainda maior | giro perto de 45° (média de 42°): 125 km a 6 dBm (49 % do bloco), 140 km a 9 dBm (24 %) e a 11 dBm (87 %) |

**Veredito.** O excesso da linha $A$ não é ruído da busca de fase: está na parte aditiva, que
escala com $1/P$. Fica entre o ruído do receptor antes do DSP (térmico, de disparo, banda do filtro
elétrico) e o ganho de ruído do equalizador. O que separa os dois são as constelações do VPI antes
do DSP (pergunta 11 do documento do grupo). Candidato concreto: a definição de banda do Bessel do
receptor (pergunta 12). O simcopy põe −3 dB em 1,15 vezes a taxa, e pela tabela da dúvida 27 os
0,28 dB de $A$ equivaleriam a uma banda de ruído uns 10 % maior. Não ajustei: seria um parâmetro
calibrado.

**Achado lateral, sobre o giro de 45°.** A grade de 256QAM com a fase estimada por polarização é a
primeira variante que produz estados girados perto de 45°, o modo de falha do VPI. Mas, na escala
em que isso acontece, os pontos que convergem têm muito mais ruído de fase que os do VPI, e 140 km
a 11 dBm, que converge no VPI, falha. Não é a cadeia do VPI do jeito que implementei. Reforça a
pergunta 1 do documento (rodar o VPI com `Format = 16QAM` na busca de fase).

### Estado das hipóteses

| # | Hipótese | Estado |
|---|---|---|
| C3 | O equalizador do VPI deixa erro em excesso | refutada: o gradiente soma ruído em todas as potências, e o excesso do VPI some em potência alta |
| C7 | A busca de fase do VPI é mais ruidosa (por polarização, grade de 256QAM) | refutada para a linha $A$: a parte de fase do erro já é igual à do VPI |
| C8 | A grade de 256QAM produz o giro de 45° | parcial: produz, mas não com a qualidade do VPI nos pontos que convergem. Aguarda a pergunta 1 |
| A1 | O ruído do receptor do VPI antes do DSP é maior (banda do filtro, convenção do ruído térmico) | aberta, aguarda as perguntas 11 e 12 |

---

## Hipótese 1 das bordas (2026-11-19): a busca de fase do VPI

Contexto: o catálogo de todos os conjuntos (`scripts/catalogo_vpi.py`, seção 13 de
`docs/comparacao_vpi.md`) mostrou que nas bordas quem decide a BER do VPI é o DSP dele perdendo a
fase, com a constelação girada perto de 45°. O notebook 09 passou a ter uma seção por hipótese.

**Como.** Variantes da busca de fase (`scripts/experimentos/busca_fase_vpi.py`) sobre os 51 sinais
guardados na saída do ADC (120 a 140 km, $2^{14}$ símbolos), com o resto do receptor cego igual.
Placar ponto a ponto contra o VPI (`scripts/experimentos/compara_hipotese.py`): nos 9 pontos em que
o VPI perde a fase, e nos 42 em que converge. Resumo em `reference/vpi/hipoteses/h1_busca_de_fase.json`.

**Mecanismo testado.** `FourthPower = Yes` e `FineBlindEstimation = Yes` do bloco do VPI lidos como
uma busca em dois estágios: fase grossa pelo ângulo de −(soma de y⁴ na janela), dividido por 4, e
busca fina de ±22,5° em torno dela. Se o ruído virar o sinal da soma, a fase grossa erra 45° e a fina
não corrige.

| Variante | Bordas | Meio | BER sobre a do VPI no meio |
|---|---|---|---|
| busca de hoje (conjunta, decisão no 16QAM) | 0/9 | 42/42 | 0,97 |
| grade de 256QAM, conjunta / por polarização | 0/9 e 4/9 | 42/42 e 22/42 | 2,67 e 3,19 |
| dois estágios, conjunta | 0/9 | 42/42 | 1,06 |
| dois estágios, por polarização | 6/9 | 28/42 | 1,24 |
| dois estágios, por polarização, só anéis de dentro e de fora, janela grossa de 43 / 37 / 32 / 29 / 25 / 21 | 0, 2, 5, 8, 8, 9 de 9 | 42, 39, 38, 35, 24, 20 de 42 | 1,12 a 1,38 |

**Veredito.**
- A grade de 256QAM fica descartada: custa 2,7 a 3,2 vezes a BER no meio.
- Os dois estágios por polarização, com o estágio grosso nos anéis de dentro e de fora, reproduzem o
  modo de falha do VPI (45°, trechos longos, frações giradas e BER perto das dele).
- Janela grossa de 32 (o `FilterLength = 32` do VPI): acerta as 5 bordas de potência baixa e 38 dos
  42 pontos do meio, sem perder a fase em potência alta. Janela de 21: acerta as 4 de potência alta
  quase nos números (125 km a 18 dBm: 55 % girado e BER 0,19, contra 52 % e 0,18), mas falha em 22
  pontos do meio. Os erros do meio ficam colados no limiar, onde a falha depende do bloco.
- A borda de potência baixa está explicada. A de potência alta (no VPI com SNR de 13,8 a 15,8 dB)
  ainda não. Próximo: o equalizador do VPI (hipótese 2).

## Hipótese 2 das bordas (2026-11-19): o equalizador do VPI

**Como.** `scripts/experimentos/equalizador_vpi.py`: cada equalizador combinado com três buscas de
fase (a de hoje e as de dois estágios da hipótese 1, janela grossa de 32 e de 21), nos mesmos 51
sinais. Equalizadores: o de hoje; o critério por raio com 20 passagens, passo 1e-3 ou 3e-4, saída com
os coeficientes de cada instante ou com os do fim; o critério por raio seguido de 4 passagens por
decisão (LMS). Resumo em `reference/vpi/hipoteses/h2_equalizador.json`.

**Veredito.** Refutada como causa das bordas. Só pelo raio, a BER no meio vai a 1,2 a 1,5 vez a do
VPI (se afasta dele), e a borda de potência alta aparece em no máximo 1 de 4 pontos. Com o LMS depois,
fica igual ao de hoje (BER no meio de 0,99 vez a do VPI).

**Ressalva sobre a hipótese 1.** `FourthPower = Yes`, `FineBlindEstimation = Yes`, `FilterLength = 32`
e `ClockRecovery = Off` têm os mesmos valores em todos os blocos `DSP_Lib` do VPI (equalizador,
relógio, frequência, busca de fase): parecem padrões do componente, que a busca de fase pode não usar.
O grupo mudou nesse bloco só `NPreSymbols`, `NbTestAngles` e `Format`. A busca em dois estágios
reproduz o modo de falha, mas a leitura dos parâmetros não está confirmada. O bloco `Clock_and_Slice`
parece só decimar (2 para 1 amostra por símbolo): provavelmente não há laço de relógio a imitar.

## Hipótese 3 das bordas (2026-11-19): a recuperação de frequência do VPI

**Como.** O VPI tem `LO_Offset = 500e6` e o bloco `FreqOffsetComp` com `Ns = 1024` e
`MultipleBlockMode = Reset`. Em `scripts/experimentos/busca_fase_vpi.py` (variantes `freq_bloco_*`),
o desvio de frequência passou a ser estimado por blocos de 1024, 512 e 256 símbolos, com a fase
contínua entre blocos, combinado com a busca de hoje e com os dois estágios (janela grossa de 32 e
de 43). Resumo em `reference/vpi/hipoteses/h3_frequencia.json`. Teste direto nos dados do VPI:
`scripts/experimentos/fronteiras_vpi.py` acha as trocas de estado (alinhado e girado de 45°) em todos
os pontos em que o VPI perde a fase e mede se elas caem nas fronteiras de algum período de bloco
(`reference/vpi/hipoteses/h3_fronteiras_vpi.json`).

**Veredito.** Refutada. O espalhamento da frequência entre blocos de 1024 fica em 4,0 MHz até 12 dBm
e 4,5 MHz de 16 dBm para cima: não cresce com a potência, e 4 MHz giram menos de 1° em 32 símbolos.
O placar não muda (busca de hoje 0/9 e 42/42, dois estágios a um ponto dos da hipótese 1). As 8236
trocas de estado do VPI não se concentram em nenhum período (em 1024: 89, acaso 105 ± 10): o VPI
não reinicia fase nem frequência por bloco. `Ns = 1024` aparece igual em cinco dos seis blocos de
DSP (a busca de fase tem `Ns = 0`): provavelmente é padrão.

## Perto do limiar, a perda de fase é um evento raro (2026-11-19)

**Composição do ruído.** `scripts/experimentos/ruido_radial_tangencial.py` decompõe o erro de cada
símbolo em radial e tangencial nos 34 pontos em que o VPI converge (2026, 125 a 140 km). A razão
tangencial/radial do VPI e a do simcopy coincidem: de 15 dBm para cima, diferença de no máximo 1,9 %
no anel de fora (130 km a 18 dBm: 1,885 contra 1,878). A SNR fica a no máximo 0,15 dB. O ruído de
fase não linear está reproduzido: a diferença das bordas não vem do canal.

**Taxa de entradas no estado girado.** Contando as trocas de estado de cada ponto do VPI ($2^{18}$
símbolos): perto do limiar o VPI entra no estado girado uma ou duas vezes no bloco inteiro e fica
girado por dezenas de milhares de símbolos (125 km a 18 dBm: uma entrada, 52 % do bloco girado). Os
9 pontos de borda do placar têm de 0,03 a 0,19 entrada por bloco de $2^{14}$: um receptor igual ao do
VPI, com blocos de $2^{14}$, perderia a fase em 0,7 deles, em média. Em função da SNR (a do simcopy no
mesmo ponto), a taxa do VPI é uma curva só para potência baixa e alta: não há mecanismo próprio da
potência alta. A borda alta parecia ter limiar de SNR mais alto por causa dos eventos raros.

**Consequências.**
- O placar "perde ou não perde num bloco de $2^{14}$" não decide perto do limiar. A janela grossa de
  21 (9/9 bordas) e a de 32 (5/5 bordas baixas) entram no estado girado muito mais vezes que o VPI:
  com a taxa do VPI, perder nas 5 bordas baixas teria chance de 1,7e-6.
- O receptor de hoje entra muito menos vezes que o do VPI: nos 14 pontos com SNR de 9,8 a 10,8 dB o
  VPI entra de 0,5 a 6,5 vezes por bloco de $2^{14}$, e o simcopy não perde a fase em nenhum.
- Os limiares de SNR da seção 13 de `docs/comparacao_vpi.md` (13,0 a 14,2 dB em potência baixa, 13,8
  a 15,8 dB em potência alta) ficam corrigidos: o limiar é uma taxa que cai depressa com a SNR.

**Próximo passo.** Medir no simcopy a taxa de entradas contra a SNR, para a busca de hoje e as
variantes da hipótese 1: com $2^{14}$ símbolos longe do limiar (0 a 5 dBm e 19 a 20 dBm, sinais a
simular) e com $2^{18}$ perto dele. O alvo é a curva do VPI.

## O modo de falha: o VPI caminha em passos de 45°, as variantes da hipótese 1 saltam 90° (2026-11-19)

**Detector de estados.** `fronteiras_vpi.phase_states` segue, em cada polarização, a fase local
contra os símbolos conhecidos (soma de y a* numa janela de 64), em passos de 45°. Conta as trocas de
estado que duram 256 símbolos e separa o tempo em múltiplos ímpares de 45° (o estado girado) do tempo
em múltiplos pares diferentes do principal (saltos de 90° ou 180°). O mesmo detector roda no VPI
(`fronteiras_vpi.py --estados`, `reference/vpi/hipoteses/estados_fase_vpi.json`) e no simcopy
(`busca_fase_vpi.py`, campos `changes`, `frac_45` e `frac_slipped`).

**VPI.** Nos 57 pontos em que perde a fase: mediana de 50 % do tempo em 45° e 35 % saltado de 90° ou
180°. A fase caminha ao acaso em passos de 45° por todos os estados. No estado girado a constelação é
limpa: em 125 km a 18 dBm a SNR com a fase corrigida pelos símbolos conhecidos é 13,92 dB nos trechos
girados e 13,86 dB nos alinhados. Nos 48 pontos convergidos, nenhuma troca.

**Correção da hipótese 1.** A "fração girada" usada nas hipóteses 1 a 3 (`rotation` em
`busca_fase_vpi.py`) soma as duas polarizações. Com a fase estimada por polarização, uma polarização
alinhada e a outra saltada de 90° somam um vetor a 45°, e o bloco parece girado de 45° sem estar. Com
o detector novo, as variantes de dois estágios por polarização (e a da grade de 256QAM por
polarização) falham por saltos de 90° (20 a 45 % do tempo saltado, 0 a 1 % em 45°). Elas não
reproduzem o modo de falha do VPI, ao contrário do que a hipótese 1 concluiu.

## Hipótese 4 das bordas (2026-11-19): uma busca de fase que não distingue 0° de 45°

**Assinaturas do VPI** (detector de estados acima): a fase caminha em passos de 45° (62 % do tempo
errado em 45°); as duas polarizações trocam de estado juntas (5688 de 6199 trocas com uma troca da
outra polarização a menos de 300 símbolos, tipicamente 3 a 5): a fase é estimada em conjunto, apesar
do `JointEstimation = No`, que é igual nos seis blocos; o estado de 45° é limpo (SNR com fase
conhecida igual à dos trechos alinhados).

**Modelo.** `busca_fase_vpi.py` (variantes `conjunta45_*`, `--fixed-equalizer`): fase comum às duas
polarizações (a diferença entre elas sai de sum y1^4 conj(y0^4)), estágio grosso pela quarta potência
dos anéis de dentro e de fora desdobrado com período de 45°, busca fina na grade do 16QAM em ±22,5°
com janela de 43. Equalizador de hoje calculado uma vez por sinal (`base_equalized`), como no VPI,
onde o equalizador não depende da busca de fase. Sinais: 99 pontos de 120 a 140 km, 0 a 20 dBm,
2^14 símbolos (48 simulados para isso, `resultados/replica_vpi/skew1amostra/longe_*.json`). Resumo em
`reference/vpi/hipoteses/h4_estados.json`. Dois defeitos achados e corrigidos no caminho: o ramo
inicial da grossa (agora escolhido por comparação local com a estimativa de 90°, num trecho limpo) e a
diferença de fase entre as polarizações (que vinha das buscas por polarização).

**Resultado.** Janela grossa de 25 a 27: taxa de trocas por faixa de SNR perto da do VPI, nos dois
lados (total 0,90 vez com 25, 0,76 vez com 27; na potência alta, com 27, 6,3 e 0,64 trocas por 2^14
em 10 e 12 dB contra 5,4 e 0,72), 56 a 57 % do tempo errado em 45° (VPI 62 %), BER no meio de 0,98 a
1,00 vez a do VPI, nenhum ponto perdido a mais no meio. A borda de potência alta sai da mesma regra
que a baixa. Não reproduz: o estado de 45° limpo (a busca fina treme, SNR com fase conhecida 1,4 a
1,6 dB abaixo nos pontos perdidos). Só a quarta potência, sem busca fina, é limpa, mas a BER no meio
vai a 1,08 vez a do VPI e a taxa cai à metade. A parte fina da busca do VPI deve ser outra (candidata:
a grade de `256QAM_constellation.txt`, que não temos).

**Próximos passos.** Pedir ao grupo `256QAM_constellation.txt` e o trecho do manual sobre o
`BlindPhaseSearch`; testar uma busca fina que não distinga 0° de 45°; simular alguns pontos perto do
limiar com 2^18 símbolos.

## Busca nos outros simuladores (2026-11-19)

Pergunta: algum dos simuladores de `outros_simuladores/` tem uma busca de fase que produza o estado de
45° do VPI, ou uma solução para ele? Resposta: não. Todos assumem a ambiguidade de 90° do QAM
quadrado, desdobram a fase com período de 90° e falham por salto de quadrante (90°), como as
variantes da hipótese 1. Lidos (os GPL só lidos, nada copiado):
- QAMpy (GPL): `bps` (ângulos em [−45°, 45°), desdobra 4φ), `bps_twostage` (Zhuge, Chen e Plant,
  OFC 2011: BPS grossa e fina em torno dela), `partition_16qam` e `phase_partition_16qam` (quarta
  potência dos anéis de dentro e de fora, e o anel do meio girado de ±(45° + atan(1/3)); Fatadin,
  Ives e Savory, PTL 2010; Muller e Mello, PTL 2010, fase conjunta em DP-16QAM). Em blocos fixos.
- OptiCommPy (GPL): `bps` (Pfau, Hoffmann e Noé 2009, ângulos em [0°, 90°)), `viterbi`
  (desdobra com período 2π/M), `cpr` desdobra 4φ de novo.
- comnum (MIT/Apache): `bps`, `viterbi_viterbi` e `viterbi_viterbi_for` (o deslocamento da quarta
  potência sai da constelação; com o deslocamento fixo de π/4 uma constelação de outra orientação
  fica girada de 45° constantes, o que não é a caminhada ao acaso do VPI), `QuadrantQam` (quadrante
  codificado diferencialmente, Pfau 2009, seção III-A).
- maiman-studio (Apache): BPS num quarto de volta, documenta o salto de quadrante como modo de falha.
- Robochameleon (só DDPLL), pyofss (só fibra), OpticalCommunications (usa o OptiCommPy).

Leitura: a caminhada em passos de 45° não é comportamento de nenhum algoritmo de livro. Reforça que é
detalhe da implementação do VPI (a grade de `256QAM_constellation.txt` ou o desdobramento dele), e o
pedido ao grupo do arquivo e do manual do `BlindPhaseSearch`.
