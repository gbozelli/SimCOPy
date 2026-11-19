# Plano para replicar o comportamento do VPI (2026-11-16)

Pedido do Gabriel: um plano para testar com os dados do VPI, hipótese por hipótese, até o simcopy
reproduzir o comportamento dele.

**Onde estamos.** No cenário de 4 canais DP-16QAM, em 120, 125, 130, 135 e 140 km, o simcopy acerta
o ruído do receptor (±0,08 dB), a imagem de quadratura, o deslocamento e a potência ótima. Sobram
cinco diferenças: o piso (0,25 a 0,37 dB abaixo), a não linearidade (0,3 a 0,6 dB abaixo), o
receptor cego (0,2 a 0,4 dB acima), e as duas quedas do VPI, em potência baixa e em potência alta,
que o simcopy não reproduz (`docs/comparacao_vpi.md`, seção 10).

**Como.** Um laço: escolher a maior diferença, escrever a previsão de uma hipótese antes de rodar,
testar no dado mais barato que a decide, registrar o veredito, atualizar o placar. As hipóteses
estão ordenadas do teste mais barato para o mais caro: primeiro medidas novas nos dados do VPI que
já temos (sem simular nada), depois o transmissor (segundos), depois o DSP sobre sinais guardados
(minutos), e por último o enlace (horas).

**Antes de começar** preciso de quatro confirmações suas (seção 8).

---

## 1. O que conta como "replicado"

O critério é fechar com o DSP cego, que é o do VPI. O DSP com símbolos conhecidos serve de
controle. As tolerâncias são o dobro da variação entre sementes medida no simcopy (±0,04 dB com
$2^{14}$ símbolos).

| # | Grandeza | Hoje (simcopy menos VPI) | Alvo |
|---|---|---|---|
| 1 | Imagem de quadratura | −5 % a 0 % | ±2 % |
| 2 | Deslocamento, com e sem os 5° | −1 % e +1 % | ±2 % |
| 3 | Ruído do receptor, $A$ | ±0,08 dB (conhecidos), −0,13 a −0,31 dB (cego) | ±0,1 dB |
| 4 | Piso, $B$ | −0,25 a −0,40 dB | ±0,1 dB |
| 5 | Resto do piso (sem imagem e sem deslocamento) | −0,50 a −0,87 dB | ±0,15 dB |
| 6 | Não linearidade, $C$ | −0,29 a −0,35 dB (conhecidos), −0,54 a −0,62 dB (cego) | ±0,15 dB |
| 7 | SNR ponto a ponto, onde o VPI converge | +0,07 a +0,30 dB (conhecidos), +0,19 a +0,44 dB (cego) | ±0,1 dB |
| 8 | Potência ótima | igual nos cinco comprimentos | igual |
| 9 | SER para a mesma SNR | não medido | ±10 % |
| 10 | Queda de baixa potência | VPI se perde abaixo de 13,5 dB; simcopy não se perdeu com 13,1 dB | mesma taxa de saltos por símbolo, dentro de um fator 2 |
| 11 | Queda de alta potência | VPI se perde em 17 ou 18 dBm em 4 de 5 comprimentos; simcopy não | reproduzir em pelo menos 3 dos 5 |
| 12 | Forma do ruído (seção 4, P2) | não medido | definido depois da etapa A |

**Replicado** quer dizer os itens 1 a 9 dentro do alvo nos cinco comprimentos, e os itens 10 e 11
em estatística. Se sobrar um item que só fecha com um dado que o VPI não exportou, ele sai do
critério e vai para a lista de pedidos ao grupo.

---

## 2. Os dados e o papel de cada conjunto

| Conjunto | Conteúdo | Papel |
|---|---|---|
| 125 km, com os 5° | 0 a 20 dBm, $2^{18}$ símbolos, X e Y | calibração |
| 120 km, com e sem os 5° | 0 a 20 dBm e 7 a 17 dBm | calibração do transmissor |
| 130, 135 e 140 km | 0 a 20 dBm | **validação**: nenhum parâmetro é ajustado olhando para eles |
| 80 km (`400gbps_CH2`) | outro cenário, configuração desconhecida | fora, até sabermos a configuração |
| Tabelas de BER do pacote, 80 a 110 km | só BER, e pela metade (dúvida 60) | conferência fraca da tendência com o comprimento |

Só existem constelações depois do DSP. Não há sinal óptico, espectro nem potência em pontos
intermediários. Por isso a etapa A tira o máximo das constelações antes de qualquer simulação.

---

## 3. Regras do laço

1. **Uma hipótese por vez**, com a previsão escrita antes de rodar: o que muda, quanto, em qual
   linha do placar, e o que não pode mudar.
2. **Três vereditos:** confirmada, refutada ou inconclusiva. Refutada também é resultado e fica no
   registro.
3. **Uma hipótese só é aceita** se fecha a diferença para a qual foi proposta, não piora outra linha
   do placar, e acerta os comprimentos de validação sem novo ajuste.
4. **Origem de cada parâmetro**, nesta ordem: está no documento do VPI, sai de um dado independente,
   ou é calibrado. No máximo um parâmetro calibrado por grandeza independente medida. Hoje há um
   calibrado, a excursão do modulador.
5. **O cenário mora no YAML.** Nenhum padrão de `src/` muda para aproximar do VPI.
6. **Registro:** `docs/replicacao_vpi_diario.md`, criado quando o laço começar, com uma entrada por
   hipótese (data, previsão, teste, números, veredito, arquivos). O placar vai no topo desse
   arquivo.
7. **Parada:** seção 7.

**O que eu faço sem perguntar:** rodar simulações e análises, criar scripts em `scripts/`, guardar
resultados em `reference/vpi/` e sinais em `resultados/` (ignorada pelo git), escrever o diário.

**O que eu pergunto antes:** mudar o YAML do cenário, mudar qualquer coisa em `src/` que não seja
correção de defeito, criar bloco novo, mexer em notebook.

---

## 4. Preparação

| # | O que | Para quê |
|---|---|---|
| P1 | Subcomando `placar` em `scripts/compara_vpi.py`: a tabela da seção 1 para qualquer conjunto de resultados, com dentro ou fora do alvo | Ver o efeito de cada hipótese em todas as linhas de uma vez |
| P2 | Medidas novas nas constelações, com o mesmo código para o VPI e para o simcopy (lista abaixo) | Separar o piso em partes sem depender de hipótese |
| P3 | Banco de sinais: a saída do receptor de cada ponto, antes do DSP, guardada em `resultados/banco_sinais/` | Testar uma hipótese de DSP em segundos, sem propagar de novo |
| P4 | Barra de erro: 8 sementes em três potências por comprimento | Saber quanto cada linha do placar varia sozinha |
| P5 | Transmitir no simcopy a sequência de símbolos do VPI (lida da constelação de referência) | Comparar o erro símbolo a símbolo. Opcional, depende da pergunta 3 da seção 8 |

Medidas de P2:

- erro separado em parte radial e parte tangencial, por anel da constelação (o ruído de fase
  aparece como excesso tangencial que cresce com o raio)
- parte do erro que é função linear dos símbolos vizinhos (interferência entre símbolos), e o que
  sobra
- autocorrelação do erro, e correlação do erro entre X e Y
- caudas da distribuição do erro
- detector de saltos de fase: quantos, e em que posição do bloco
- SER contra SNR
- erro ao longo do bloco (transiente no começo, periodicidade)

---

## 5. As hipóteses

Cada linha traz a base, a previsão e o teste. O custo é o tempo de máquina de um teste.

### Etapa A. Medidas nos dados do VPI que já temos (minutos, sem simular)

| # | Pergunta | O que decide |
|---|---|---|
| A1 | Nos arquivos em que o VPI se perde, a constelação gira 90° em trechos, ou vira ruído? | Se as quedas são saltos da busca de fase (itens 10 e 11) ou falha de outra etapa (equalizador, sincronismo) |
| A2 | Quantos saltos por símbolo o VPI tem em cada arquivo marginal? | A curva de taxa de saltos contra SNR, que é o alvo do item 10 |
| A3 | Quanto do resto do piso do VPI é tangencial (fase) e quanto é radial? | Se o que falta no simcopy é ruído de fase residual ou ruído branco |
| A4 | Quanto do resto do piso do VPI é interferência entre símbolos? | Se o que falta é do transmissor ou do equalizador (determinístico) ou é ruído |
| A5 | O erro do VPI tem memória? X e Y erram juntos? | Fase não linear comum às duas polarizações, e a janela efetiva da busca de fase |
| A6 | A SER do VPI é a de ruído gaussiano com a SNR medida? | Item 9 |
| A7 | O erro do VPI muda ao longo do bloco? | Transiente do equalizador, processamento em blocos |
| A8 | X e Y têm a mesma SNR no VPI? | PMD, extinção dos divisores de polarização |

As mesmas medidas rodam nos símbolos do simcopy. A diferença em cada uma vira linha do placar e
corta a lista abaixo: se A3 e A4 mostrarem que o que falta é ruído branco, B1 é a candidata; se for
tangencial, C1 e C2; se for interferência, B3 e C3.

### Etapa B. Transmissor (um canal, sem fibra, segundos por teste)

| # | Hipótese | Base | Previsão | 
|---|---|---|---|
| B1 | A potência de saída do transmissor do VPI é cerca de 1 dB menor que a do simcopy | **[Fato medido]** quase metade do resto do piso é ASE do booster. **[Dúvida 59]** os 6 dB de `InsertionLoss` podem ser de cada MZM ou do modulador inteiro | resto do piso sobe 0,0009 nos cinco comprimentos; $A$, $C$, imagem e deslocamento não mudam |
| B2 | O simcopy não modela a extinção do divisor e do combinador do transmissor | **[Fato: VPI]** `SplitterExtinctionRatioX/Y` = `CombinerExtinctionRatioX/Y` = 30 dB, além dos 30 dB de cada MZM | muda o deslocamento ou acrescenta piso; a excursão calibrada muda junto |
| B3 | O skew, o tipo do filtro S21 ou o pulso do VPI diferem do que implementei | **[Inferência]** li "Gaussian de ordem 3" e Bessel do documento sem conferir a definição | muda só a parte de interferência entre símbolos (A4) |
| B4 | `Renormalize = On` fixa a potência de saída do transmissor do VPI | **[Fato: VPI]** o parâmetro existe; não sei para que valor ele normaliza | a ASE do booster deixa de depender da excursão; equivale a B1 com outro número |

### Etapa C. DSP, sobre o banco de sinais (segundos a minutos por teste)

| # | Hipótese | Base | Previsão |
|---|---|---|---|
| C1 | `NPreSymbols = 21` é meia janela: a média de fase do VPI usa 43 símbolos | **[Fato medido]** $C$ do VPI fica entre o do simcopy com janela de 21 e o com janela de 81 | $C$ sobe para perto do VPI; a taxa de saltos cai (o que iria contra o item 10) |
| C2 | A média móvel do VPI não é a minha (forma da janela, uso das duas polarizações) | **[Fato: VPI]** `Filtering = MovingAverage`, `CombinedPol = Yes` | muda $A$ e $C$ do receptor cego, não muda o de símbolos conhecidos |
| C3 | O equalizador do VPI deixa um erro em excesso que o meu não deixa | **[Fato: VPI]** `AlgTDE = MMA`, `Mode = Decision-Directed`, `Iterations = 20`, `StepSize = 1e-3`, `InitTaps = RRC`. **[Inferência]** ele adapta por gradiente, sem o ajuste por mínimos quadrados em bloco que o simcopy faz. **[Fato medido]** só com o critério do módulo o simcopy perde 0,7 dB | $A$ e o piso sobem no cego; a SNR cai 0,2 a 0,4 dB |
| C4 | A ordem das etapas do VPI muda o resultado | **[Fato: VPI]** equalizador, depois Stokes, depois desvio de frequência e fase. No simcopy o Stokes vem antes | pequena, ou nenhuma |
| C5 | A queda vem do estimador de desvio de frequência ou do equalizador, não da busca de fase | depende de A1 | reproduz a posição da queda sem mexer na janela de fase |
| C6 | `ConstellationAlignment = On` do `BER_4D` altera a constelação exportada | **[Fato: VPI]** o parâmetro existe | muda ganho, rotação ou deslocamento medidos, não a SNR do resto |

O documento do VPI tem `ClockRecovery = Off`: a recuperação de relógio não entra na lista.

### Etapa D. Enlace (15 minutos a 1 hora por teste)

| # | Hipótese | Base | Previsão |
|---|---|---|---|
| D1 | A sequência PRBS de ordem 7 do VPI (período de 127 bits) muda a não linearidade | **[Fato: VPI]** `PRBS_Order = 7`. No simcopy os bits são aleatórios | muda $C$ e a forma do ruído não linear (A5) |
| D2 | O valor de $C$ do simcopy ainda depende do passo, das amostras por símbolo ou do tamanho do bloco | conferência numérica, sem hipótese física | $C$ não muda mais de 0,05 dB ao dobrar cada um |
| D3 | O filtro do mux e do demux do VPI não é a minha supergaussiana | **[Inferência]** "Gaussian de ordem 3, 70 GHz": supus que 70 GHz é a banda de 3 dB | muda a diafonia e a parte dos vizinhos em $C$ |
| D4 | A mistura de quatro ondas, que o VPI desliga, pesa em $C$ | **[Fato: VPI]** `FWM = No`. O sinal da diferença é o contrário do esperado | só dá para medir de forma indireta. Fica por último |
| D5 | A realização de PMD pesa | uma realização só até agora | variação entre sementes de PMD menor que 0,05 dB |

### Etapa E. As quedas (estatística)

| # | O que | Como |
|---|---|---|
| E1 | Taxa de saltos por símbolo contra SNR no simcopy | 16 sementes de $2^{14}$ símbolos em cada ponto marginal, que somam os $2^{18}$ do VPI |
| E2 | A queda de 17 e 18 dBm | o mesmo, em potência alta, com a configuração de DSP que sair da etapa C |

---

## 6. Ordem e custo

| Passo | Conteúdo | Tempo de máquina |
|---|---|---|
| 1 | P1, P2 e a etapa A inteira | menos de 1 hora |
| 2 | Relatório da etapa A: placar com as linhas novas e a lista de hipóteses que sobrou | parada para você ler |
| 3 | P3 e P4 | cerca de 2 horas |
| 4 | Etapas B e C, uma hipótese por vez | segundos a minutos cada |
| 5 | Rodada completa nos cinco comprimentos com o que foi aceito | 1 hora |
| 6 | Etapa D, só para as linhas que ainda estiverem fora | 1 a 4 horas |
| 7 | Etapa E | 2 a 3 horas |
| 8 | Fechamento: seção nova em `docs/comparacao_vpi.md` e proposta de teste de regressão para o grupo | |

Paradas para você: depois do passo 2, a cada hipótese aceita que mude o YAML, e no fim.

---

## 7. Quando parar

- **Replicado:** todas as linhas do placar dentro do alvo. O cenário é congelado e a comparação
  vira proposta de teste de regressão (issue para o grupo).
- **Esgotado:** as hipóteses acabaram e sobra diferença. O resultado é a lista do que só fecha com
  um dado do VPI que não temos, com o tamanho de cada diferença que sobrou.
- **Preso:** três hipóteses seguidas inconclusivas na mesma linha do placar. Paro nessa linha e
  sigo para a próxima.

---

## 8. O que preciso de você antes de começar

1. **As regras da seção 3:** posso seguir o laço sem perguntar a cada hipótese, parando só nos
   pontos marcados?
2. **O banco de sinais** em `resultados/banco_sinais/` (ignorada pelo git, cerca de 100 MB): pode?
3. **Código só de experimento:** as hipóteses C3, D1 e P5 pedem um equalizador como o do VPI, uma
   fonte de bits PRBS e uma fonte de bits lida de arquivo. Posso criá-los em `scripts/experimentos/`,
   fora de `src/`, e só propor a entrada em `src/` se a hipótese for confirmada?
4. **Do lado do VPI,** se o grupo conseguir: a potência óptica na entrada do booster, uma rodada
   sem fibra e uma com um canal só (`docs/comparacao_vpi.md`, seção 10.5). Cada um desses elimina
   uma etapa inteira deste plano.
