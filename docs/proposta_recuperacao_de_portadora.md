# Recuperação de portadora: proposta para discussão (2026-11-18)

Nasceu do notebook 08, que diz que a recuperação de portadora "ainda não existe no simcopy". Ao
medir, o quadro é outro: o DSP (`rx.dsp`) já recupera a frequência e a fase da portadora, e funciona
com o oscilador local livre. O que falta é mais estreito do que parecia.

**Proposta em uma frase:** não escrever um algoritmo novo. Mostrar o que já existe no notebook 08,
ampliar a faixa de desvio de frequência do modo cego (hoje Rs/8) e reduzir os saltos de fase com
SNR baixa, e ajustar o aviso da validação.

Nada aqui foi implementado. A seção 5 lista o que precisa ser decidido.

---

## 1. O problema, em palavras simples

No receptor coerente o sinal é misturado com um oscilador local. Se o oscilador é outro laser (não
uma cópia do laser do transmissor), ele tem uma frequência um pouco diferente e uma fase que vai
passeando. O resultado é que a constelação gira: no notebook 08, a 500 MHz e a 10 GBd, ela dá 0,05
volta por símbolo e vira um anel.

Para desfazer isso o receptor precisa estimar, sem saber o que foi transmitido:

1. **o desvio de frequência**, que gira a constelação sempre no mesmo sentido. O método do simcopy
   eleva os símbolos à quarta potência, o que apaga a modulação do QPSK e deixa uma raia no
   espectro, em 4 vezes o desvio;
2. **a fase**, que oscila devagar por causa do ruído de fase dos lasers. A busca cega de fase testa
   vários ângulos e fica com o que deixa os símbolos mais perto dos pontos ideais, médio numa janela.

## 2. Como está hoje

**[Fato: código em `src/simcopy`]**

| Peça | O que faz | Limite |
|---|---|---|
| `carrier_recovery` (`kernels/dsp.py`), usado no modo `blind` | raia da quarta potência para a frequência e busca cega de fase (64 ângulos, janela de 21 símbolos) | frequência até **Rs/8**; saltos de fase com SNR baixa (a docstring do `DSPConfig` já avisa, abaixo de ~13 dB) |
| `timing_and_frequency` (`kernels/dsp.py`), no modo `data_aided` | produto diferencial com os símbolos conhecidos | usa os dados transmitidos; vale para |desvio| < Rs/2 |
| `FrequencyOffsetCompensator` (`components/rx/dsp.py`) | a mesma raia da quarta potência, antes do filtro casado | desligado por padrão; procura até Rs/16; frágil no 16QAM |
| `measure_coherent` (caminho sem DSP, o do notebook 08) | tira só um ganho complexo | uma rotação constante: não segue a frequência nem a fase |

A regra de validação `CoherentNeedsCarrierRecovery` avisa quando há oscilador livre ou desvio de
frequência, mas **só para canais sem `rx.dsp`** (`validation/rules.py`).

## 3. O que medi (2026-11-18)

QPSK e BPSK de um canal, dos exemplos `coerente_qpsk.yaml` e `coerente_bpsk.yaml`, a 10 GBd, com o
oscilador local livre (largura de linha de 500 kHz). BER contada, com os erros entre parênteses.

**3.1 Sem DSP, o anel.** A −30 dBm, com 500 MHz de desvio, a BER é 0,493 no QPSK e 0,488 no BPSK.
No BPSK as medidas do resumo perdem o sentido: EVM de 4845 % e instante de decisão de 115715
amostras, porque a constelação gira e o ganho complexo sai quase zero.

**3.2 Com DSP, 500 MHz não custa nada.** Comparando o oscilador livre com o travado **com o mesmo
DSP** (o DSP já melhora a BER sozinho, por causa do equalizador, então a comparação justa é essa):

| Caso | QPSK −44 dBm | QPSK −42 dBm | BPSK −47 dBm | BPSK −45 dBm |
|---|---|---|---|---|
| travado, `data_aided` | 3,2e-3 (207) | 3,5e-4 (23) | 2,5e-3 (82) | 4,6e-4 (15) |
| livre 500 MHz, `data_aided` | 2,9e-3 (192) | 3,4e-4 (22) | 3,0e-3 (97) | 2,7e-4 (9) |
| travado, `blind` | 3,0e-3 (197) | 1,4e-4 (9) | 2,6e-3 (86) | 1,8e-4 (6) |
| livre 500 MHz, `blind` | 3,0e-3 (199) | 2,6e-4 (17) | 2,3e-3 (76) | 2,1e-4 (7) |

As diferenças ficam dentro da estatística de contagem.

**3.3 A faixa do modo cego é Rs/8.** QPSK a −40 dBm, 10 GBd (Rs/8 = 1,25 GHz):

| Desvio | `data_aided` | `blind` |
|---|---|---|
| 0,5 e 1,0 GHz | 0 erros | 0 erros |
| 1,2 GHz | 3 erros | 0 erros |
| 1,4 GHz | 1 erro | **BER 0,50, 195 saltos de fase** |
| 2, 3 e 4 GHz | 0, 1 e 3 erros | **BER 0,50, ~190 saltos** |

O cego cai logo depois de Rs/8. O `data_aided` segue até 4 GHz (o limite do código é Rs/2). Com 2 GHz
de desvio o cego falha em todas as potências medidas (−48 a −42 dBm).

**3.4 Saltos de fase com SNR baixa, no modo cego.** QPSK a −48 dBm: 20 saltos com o oscilador
travado e 28 com o livre (BER 0,38 e 0,48). BPSK a −51 dBm: o travado vai bem (3,6e-2), mas o livre
perde a portadora (BER 0,50, 136 saltos). A −46 dBm (QPSK) e −49 dBm (BPSK) o livre já recupera.

## 4. Proposta

**[Proposta]** Quatro etapas, cada uma útil sozinha.

1. **Mostrar o que já existe no notebook 08.** Rodar o mesmo oscilador livre de 500 MHz com
   `dsp=DSPConfig(mode="blind", cd_compensation=False, stokes_demux=False)` e mostrar o anel virando
   pontos. Trocar a frase "ainda não existe" por "existe no DSP, e o caminho sem DSP não a tem".
2. **Ajustar a validação.** O aviso diz hoje que "o receptor ainda não recupera a portadora" e só
   sugere o oscilador travado. Passa a sugerir também `rx.dsp`, e um segundo aviso para o modo cego
   quando o desvio passar de Rs/8.
3. **Ampliar a faixa do modo cego.** **[Dedução]** A raia da quarta potência é medida a 1 amostra
   por símbolo, onde 4 vezes o desvio precisa caber em ±Rs/2, daí o Rs/8 (o código procura até
   `1/8` da taxa em `carrier_recovery`). Estimar a 2 amostras por símbolo (o `ADC` já entrega isso)
   dobraria a faixa para Rs/4, e um estimador diferencial pode ir além. Alvo a combinar (seção 5).
4. **Menos saltos com SNR baixa.** **[Dedução]** A 10 GBd e 1 MHz de largura de linha somada, a fase
   anda cerca de 0,025 rad por símbolo, então uma janela bem maior que 21 símbolos ainda acompanha o
   ruído de fase e reduz os saltos. Medir a taxa de saltos contra a SNR para várias janelas
   (21, 43, 81, 161) antes de escolher o padrão. Se não bastar, a codificação diferencial no
   transmissor elimina o efeito de um salto, mas muda o mapeamento e fica para outra discussão.

## 5. O que decidir

1. **Escopo.** Minha recomendação: as etapas 1 e 2 agora, que são pequenas, e as 3 e 4 só se o grupo
   precisar de desvios acima de Rs/8 ou de operar com SNR baixa.
2. **Faixa de desvio.** Qual desvio de frequência o simcopy deve aguentar no modo cego? O VPI usa
   500 MHz (dentro da faixa atual). Um laser livre real pode passar de 1 GHz, o que está acima de
   Rs/8 em taxas baixas como 10 GBd. Vale perguntar ao professor.
3. **Texto do notebook 08.** Ajustar a frase e, se quiser, acrescentar a figura com o DSP (etapa 1).
4. **`frequency_offset` do `DSPConfig`.** Continua desligado por padrão (é frágil no 16QAM e a
   medida já estima o desvio com os símbolos)?
5. **Janela da busca de fase.** O padrão de 21 símbolos vem do VPI (43 de janela medida nos dados
   dele, a 59,84 GBd). Para o 10 GBd dos exemplos de um canal, o valor deve depender da taxa?
6. **Ordem com o relógio.** A recuperação de relógio (`docs/proposta_recuperacao_de_relogio.md`) vem
   antes da portadora na cadeia. Os dois laços mexem na mesma janela, e vale decidir os dois juntos.

## 6. Como testar

Insumo para a issue de testes, que fica com o grupo. A vantagem aqui é que o simulador conhece o
desvio verdadeiro. **[Proposta]**

- **Erro da estimativa de frequência:** com o desvio imposto (`rx.lo_frequency_offset`) conhecido,
  o desvio estimado deve ficar perto dele, sem ruído e com ruído.
- **Livre contra travado, mesmo DSP:** a BER deve ser a mesma dentro da estatística (a tabela 3.2).
- **Varredura de desvio** de 0 a Rs/2: o cego deve funcionar até o limite escolhido e falhar além
  dele de forma clara (BER 0,5 e muitos saltos).
- **Saltos de fase contra SNR**, por janela, no cego.
- **Sem sorteio extra:** a mesma entrada dá a mesma saída.

## 7. Esforço

Estimativa minha, ainda não medida: etapa 1, 0,5 h. Etapa 2, 0,5 h. Etapa 3, 3 a 4 h. Etapa 4, 2 a
3 h (a maior parte é a varredura).
