# Plano de validação cruzada com outros simuladores (2026-11-16)

Pedido do Gabriel: um plano para verificar o simcopy contra os outros simuladores e validar os
dados. É o segundo plano. O primeiro, `docs/plano_replicar_vpi.md`, trata do VPI.

**A ideia.** Para cada bloco do simcopy, dar a mesma entrada a uma implementação independente e
comparar a saída, número por número, com uma tolerância escrita antes. Onde existir forma fechada
no Agrawal, ela decide quem está certo.

**O que isto não é.** Não é uma disputa de recursos (essa está em `docs/comparacao_simuladores.md`)
e não usa os outros simuladores como verdade: dois programas podem errar igual.

**Antes de começar** preciso de quatro decisões suas (seção 8). A principal é a da licença.

---

## 1. Três referências, da mais forte para a mais fraca

| Referência | O que prova | Limite |
|---|---|---|
| Forma fechada (Agrawal) | o bloco faz a conta certa | só existe para casos simples |
| Outro simulador, com a mesma entrada | duas implementações independentes concordam | os dois podem ter o mesmo erro, ou convenções diferentes |
| Dados do VPI | o sistema inteiro se comporta como o de referência do grupo | só há constelações depois do DSP |

Quando o simcopy e o outro simulador discordam, a ordem de investigação é: convenção, passo e
convergência, parâmetro, e só então defeito. Quem decide é a forma fechada. Se não houver, o caso
fica registrado como divergência aberta, sem vencedor.

---

## 2. Licença e isolamento

**[Fato: `docs/comparacao_simuladores.md`, seção 2]**

| Projeto | Licença | Linguagem | Serve para |
|---|---|---|---|
| OptiCommPy | GPL-3.0 | Python | quase tudo: é o de maior sobreposição |
| QAMpy | GPL-3.0 | Python | DSP cego e métricas |
| comnum | MIT ou Apache-2.0 | Rust, com interface Python | fibra, amplificador, DSP, modelos GN e EGN |
| Maiman Studio | Apache-2.0 | Python | laser por equações de taxa, IMDD, amplificador |
| PyOFSS | GPL-3.0 | Python 2, parado desde 2012 | propagação de um pulso |
| Robochameleon | GPL-3.0 | MATLAB | DSP |
| OpticalCommunications | GPL-3.0 | notebooks | usa o OptiCommPy: não acrescenta referência |

Você já tem a regra de nunca importar o BIFROST por ser GPL. Proponho a mesma regra, escrita para
todos:

1. Nenhum simulador externo é importado por código do simcopy, nem em `src/`, nem em `scripts/`,
   nem em notebook.
2. Os externos ficam instalados num ambiente próprio em `rascunhos/externos/`, que o git ignora.
   Os scripts que os chamam ficam lá também e não vão para o repositório.
3. Entre os dois lados só passam arquivos de dados: a entrada que o simcopy gerou e a saída que o
   externo calculou.
4. No repositório ficam os vetores de entrada, as saídas, o script de comparação (que só lê
   arquivos) e a tabela de resultados.

Começo pelos de licença permissiva (comnum e Maiman) onde eles cobrem o bloco. O OptiCommPy e o
QAMpy entram se você aprovar.

---

## 3. Qual referência para qual bloco

Tolerâncias: saída determinística, erro quadrático relativo do campo; saída com ruído, estatística
dentro da barra de erro de N realizações.

| Bloco do simcopy | Forma fechada | Simulador | O que se compara | Tolerância |
|---|---|---|---|---|
| Mapeamento e decisão (BPSK a 64QAM, PAM4) | não precisa | OptiCommPy, QAMpy, comnum | mesmos bits, mesmos símbolos, a menos da convenção de bits | exato |
| Métricas: BER, SER, EVM, MI, GMI | BER teórica com ruído gaussiano | OptiCommPy, QAMpy | os mesmos símbolos recebidos e transmitidos nos dois | 1e-6 relativo; 0,001 bit em MI e GMI |
| Pulsos de cosseno levantado e raiz | resposta analítica | OptiCommPy | resposta em frequência | 1e-6 |
| MZM e modulador IQ | Agrawal, Eq. 10.1.5 | OptiCommPy | campo de saída para o mesmo sinal elétrico | 1e-9 |
| Laser CW com ruído de fase | espectro lorentziano | OptiCommPy | largura de linha medida | 3 % |
| Laser por equações de taxa | Agrawal, capítulo 3 | Maiman Studio | curva de potência por corrente, oscilação de relaxação, chirp | 2 % |
| Fibra linear ($\beta_2$, $\beta_3$, atenuação) | Agrawal, seção 2.4 e Eq. 2.4.22 | OptiCommPy, comnum, PyOFSS | campo de saída para o mesmo pulso | 1e-9 |
| Fibra não linear, uma polarização | Agrawal, seção 2.6.2 (fase não linear) | OptiCommPy, comnum, PyOFSS | campo de saída, com o passo reduzido até convergir | 1e-4 |
| Fibra com duas polarizações (Manakov) e PMD | estatística do atraso diferencial | OptiCommPy, comnum | campo (PMD desligada); média e distribuição do atraso | 1e-4; 3 % |
| Amplificador | Agrawal, Eq. 7.2.15 | OptiCommPy, comnum, Maiman | densidade da ASE e OSNR | 1 % |
| Mux e demux | resposta do filtro | OptiCommPy | resposta em frequência | 1e-6 |
| Fotodetector | Agrawal, Eq. 4.4.5, 4.4.9 e 4.6.11 | OptiCommPy, Maiman | variâncias dos ruídos e fator Q | 2 % |
| Receptor coerente sem ruído | $R\,E_s E_{LO}^*$ | OptiCommPy | correntes para os mesmos campos | 1e-9 |
| Compensação de dispersão | inversa exata da fibra linear | OptiCommPy | campo | 1e-9 |
| ADC | fórmula do ruído de quantização | OptiCommPy | ruído acrescentado | 1 % |
| DSP cego (equalizador e busca de fase) | não há | QAMpy, OptiCommPy | SNR depois do DSP, sobre as mesmas amostras guardadas | 0,1 dB |
| Sistema sem fibra, 16QAM | Agrawal, Eq. 10.3.2 e 10.3.7 (BPSK e QPSK); BER de QAM quadrado | OptiCommPy | curva de BER contra SNR | barra de erro |
| IMDD com OOK e PAM4 | Agrawal, Eq. 4.6.10 | OptiCommPy, Maiman | BER contra potência recebida | barra de erro |
| Não linearidade em WDM (o termo $C$) | modelo GN, que não está no Agrawal | comnum (modelos GN e EGN) | coeficiente de ruído não linear | 0,5 dB |

Na comparação de recursos várias células destes projetos ficaram com "?": eu só li as páginas, não
rodei nada. O primeiro passo de cada linha é confirmar que o externo tem mesmo o bloco.

---

## 4. Protocolo

### 4.1 Antes de qualquer bloco: o mapa de convenções

A maior parte das divergências entre simuladores é convenção. Para cada externo, um caso trivial
fixa cada item, e nenhum bloco é comparado antes disso:

| Convenção | No simcopy |
|---|---|
| Sinal da portadora | $\exp(-i\omega_0 t)$, a do Agrawal: o campo de outro simulador pode ser o conjugado |
| Eixo de frequência da FFT | a frequência óptica relativa é o negativo do eixo da FFT |
| Sinal de $\beta_2$ e de $D$ | $D$ positivo é dispersão anômala, $\beta_2$ negativo |
| Potência | soma das duas polarizações; campo em $\sqrt{\mathrm{W}}$ |
| Densidade de ruído | ruído de disparo pela densidade bilateral do Agrawal; ruído térmico por densidade unilateral de corrente, em A/√Hz; ASE por polarização |
| Fator não linear | 8/9 de Manakov, sempre |
| Mapeamento de bits | Gray, com a ordem de bits do código do professor |
| Instante de amostragem | achado por correlação com os símbolos transmitidos, sem recuperação de relógio |

### 4.2 Vetores de teste

- O simcopy gera a entrada de cada caso e grava em `reference/cruzada/<bloco>/` com um manifesto:
  parâmetros em unidades do SI, data e versão do simcopy.
- O externo lê a entrada, roda o bloco dele e grava a saída com a versão dele.
- Blocos determinísticos comparam vetores. Blocos com ruído comparam estatísticas de N realizações,
  porque as sementes não são compatíveis.

### 4.3 Comparação

Um script, `scripts/validacao_cruzada.py`, lê os arquivos, calcula a diferença, aplica a tolerância
da seção 3 e escreve a tabela em `docs/validacao_cruzada_resultados.md`: bloco, referência,
diferença medida, tolerância, dentro ou fora.

### 4.4 Quando discordar

1. Convenção (seção 4.1).
2. Passo, amostras por símbolo e tamanho da janela: reduzir até a diferença parar de mudar.
3. Parâmetro: conferir o manifesto dos dois lados.
4. Defeito. Se for do simcopy, corrijo e vira issue de teste para o grupo. Se for do externo, fica
   anotado, com o caso que o mostra.

---

## 5. Validar os dados

São três conjuntos de dados diferentes.

**a. Os números que o simcopy já publicou** (notebooks e guias). Cada número de destaque ganha uma
linha numa tabela: onde aparece, variação entre sementes, se ainda muda ao dobrar o passo, as
amostras por símbolo e o tamanho da janela, e contra o quê foi conferido. Um número sem nenhuma
conferência fica marcado.

**b. Os dados do VPI.** Uma ficha do conjunto, com o que já achamos (as tabelas contam I e Q em
separado, a BER delas é a metade, a pasta de 140 km duplicada, a constelação de referência não é a
saída do transmissor, a pasta de 80 km é outro cenário) e o que falta conferir:

- a sequência de símbolos é a mesma em todos os arquivos, e tem o período que a ordem 7 do PRBS dá
- X e Y têm a mesma potência e a mesma SNR
- a normalização das constelações é a mesma em todas as potências
- quais arquivos são de DSP convergido e quais não são

**c. Triângulos.** Para as poucas grandezas em que há três fontes independentes, as três lado a
lado:

| Grandeza | Forma fechada | simcopy | Terceira fonte |
|---|---|---|---|
| BER contra SNR sem fibra | BER de QAM quadrado | sim | OptiCommPy |
| Ruído do receptor ($A$) | ruído térmico e de disparo dos oito fotodiodos | sim | VPI |
| Não linearidade ($C$) | modelo GN | sim | VPI |
| Sensibilidade do IMDD | Agrawal, seção 4.6 | sim | Maiman ou OptiCommPy |

---

## 6. Ordem

| Fase | Conteúdo | Por que nesta posição |
|---|---|---|
| 0 | Sua decisão sobre a licença, o ambiente em `rascunhos/externos/`, o mapa de convenções | sem isso nada compara |
| 1 | Funções puras: mapeamento, métricas, pulsos, filtros | baratas e exatas: qualquer diferença é convenção ou defeito |
| 2 | Blocos determinísticos: MZM e IQ, fibra linear, compensação de dispersão, receptor sem ruído | a mesma entrada dá a mesma saída |
| 3 | Fibra não linear, uma e duas polarizações | depende de convergência no passo |
| 4 | Blocos com ruído: amplificador, fotodetector, laser, ADC | pede estatística |
| 5 | DSP sobre amostras guardadas | usa o banco de sinais do plano do VPI |
| 6 | Sistemas: curva sem fibra, um ponto de IMDD, o termo $C$ contra o modelo GN | soma tudo |
| 7 | Validação dos dados (seção 5) e lista de issues de teste para o grupo | fecha |

As fases 1 e 2 cabem numa sessão cada. A fase 5 se liga à etapa C de `docs/plano_replicar_vpi.md`:
rodar o DSP do QAMpy sobre o mesmo sinal é uma terceira opinião sobre os 0,2 a 0,4 dB que separam o
receptor cego do simcopy do VPI.

---

## 7. O que sai no fim

- `docs/validacao_cruzada_resultados.md`: uma linha por bloco, com a diferença medida.
- `reference/cruzada/`: os vetores e as saídas, para qualquer um repetir.
- A ficha dos dados do VPI e a tabela dos números publicados (seção 5).
- A lista de issues para o grupo: cada caso conferido vira um teste automático, que é do grupo.
- As divergências abertas, se houver, cada uma com o caso que a reproduz.

O plano termina quando toda linha da seção 3 tem um resultado: dentro, fora com causa achada, ou
sem referência disponível.

---

## 8. O que preciso de você antes de começar

1. **Licença:** posso instalar o OptiCommPy e o QAMpy (GPL-3.0) em `rascunhos/externos/`, com as
   quatro regras da seção 2? Sem eles o plano fica com o comnum e o Maiman, e perde as linhas de
   métricas, DSP e receptor coerente.
2. **Robochameleon** pede MATLAB. Você tem licença? Se não, fica de fora.
3. **PyOFSS** é Python 2 e está parado desde 2012. Proponho deixar de fora: o OptiCommPy e o comnum
   cobrem a mesma propagação.
4. **Modelo GN:** não está no Agrawal. Aceita como referência para o termo $C$, com a fonte citada?
