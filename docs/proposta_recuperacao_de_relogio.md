# Recuperação de relógio: proposta para discussão (2026-11-18)

Nasceu de uma pergunta sobre o notebook 08: por que o instante de decisão do receptor coerente fica
5 amostras depois do centro do símbolo, e como o simulador sabe disso. A resposta é que ele sabe
porque conhece os símbolos transmitidos, e um receptor real não conhece.

**Proposta em uma frase:** manter a medida de hoje como está e acrescentar um laço de recuperação
de relógio (Gardner) entre o ADC e o receptor cego, que descubra o instante de decisão só com o
sinal recebido, com o desvio de relógio do ADC como parâmetro para ter o que recuperar.

Nada aqui foi implementado. A seção 6 lista o que precisa ser decidido.

---

## 1. O problema, em palavras simples

Para decidir qual símbolo chegou, o receptor olha o sinal em um instante por símbolo. Se olhar um
pouco antes ou depois do ponto certo, o ponto recebido sai deslocado e aparecem mais erros.
Recuperar o relógio é achar esse instante, e continuar achando se ele andar, usando só o próprio
sinal.

No simulador o transmissor e o receptor usam a mesma grade de tempo, então o instante certo não
anda: só existe um atraso fixo (o filtro do receptor, a fibra). Num sistema real os dois relógios
são diferentes, e o instante certo desliza aos poucos.

## 2. Como está hoje

**[Fato: código em `src/simcopy`]**

| Peça | O que faz | Limite |
|---|---|---|
| `symbol_timing` (`kernels/demodulation.py:210`) | correlaciona o sinal com os símbolos transmitidos e devolve o atraso | usa os dados transmitidos; devolve um número **inteiro** de amostras, um só para a janela toda |
| `timing_and_frequency` (`kernels/dsp.py:59`) | a mesma ideia, no receptor com símbolos conhecidos | idem |
| `BlindReceiver` (`components/rx/adc.py`) | equalizador adaptativo a 2 amostras por símbolo, que aceita "qualquer fase de amostragem" e "absorve o atraso" (docstring de `adaptive_equalizer`) | resolve um atraso parado, não um que anda |
| `ADC` (`components/rx/adc.py`) | amostra a cada `sps/2` amostras da grade, a partir de um `sampling_offset` inteiro | o relógio do ADC é o da simulação; "não há recuperação de relógio" (docstring) |

**[Fato medido, 2026-11-18, notebook 08, sem ruído]** O filtro do receptor (Bessel de ordem 4,
11,5 GHz) atrasa o sinal 4,7 amostras. O `symbol_timing` arredonda para 5, e a EVM é 3,63 %. Se a
amostragem usar 4,7, a EVM cai para 1,88 %. Então metade da EVM do notebook era arredondamento do
instante, não interferência entre símbolos.

**[Fato: documentos]**
- O documento de parâmetros do VPI tem `ClockRecovery = Off` (`docs/plano_replicar_vpi.md`): não há
  referência do VPI para comparar.
- `docs/sistema_coerente.md` §6 e `docs/comparacao_simuladores.md` já listam a recuperação de
  relógio como o que falta no DSP.
- **D-4** (`docs/decisoes.md`): sem reamostragem fracionária na v1. O precedente é o skew do IQ (M7),
  aceito por fase linear na FFT, que é um filtro de atraso e não uma reamostragem.

## 3. Os dois algoritmos, em palavras simples

Os dois fazem o mesmo trabalho em três passos: (1) calculam, a cada símbolo, um **erro de tempo**
(o quanto estamos adiantados ou atrasados), (2) passam esse erro por um **filtro de laço**, que age
como um termostato e muda o instante de amostragem aos poucos, e (3) usam um **interpolador** para
ler o sinal entre duas amostras. Só o passo (1) muda.

| | Gardner | Mueller-Müller |
|---|---|---|
| Amostras por símbolo | 2 (uma no símbolo, uma entre dois símbolos) | 1 |
| Precisa dos símbolos? | não | precisa de decisões (estimativas dos símbolos) |
| Erro de tempo | $e_k = \mathrm{Re}\{(y_k - y_{k-1})\,y_{k-1/2}^*\}$ | $e_k = \mathrm{Re}\{a_{k-1}^*\,y_k - a_k^*\,y_{k-1}\}$ |
| Ponto forte | não depende da fase da portadora | funciona com uma amostra por símbolo |
| Ponto fraco | o sinal do detector diminui com o roll-off do pulso (o do notebook é 0,18) | pede portadora já recuperada e é sensível à interferência entre símbolos |

Aqui $y_k$ é a amostra no instante do símbolo $k$, $y_{k-1/2}$ a amostra do meio entre dois símbolos
e $a_k$ a decisão. **[Fato: fora do repositório, conferir na fonte]** Gardner, "A BPSK/QPSK
timing-error detector for sampled receivers", IEEE Trans. Commun., 1986. Mueller e Müller, "Timing
recovery in digital synchronous data receivers", IEEE Trans. Commun., 1976.

## 4. Proposta

**[Proposta]** Quatro etapas, cada uma útil sozinha.

1. **Instante fracionário na medida de hoje** (pequena). O `symbol_timing` passa a interpolar o pico
   da correlação e a amostragem usa o atraso fracionário. Não é recuperação de relógio: é a medida
   do notebook ficar certa (3,6 % para 1,9 %). O deslocamento fracionário pode usar a fase linear
   na FFT, como no skew.
2. **Laço de Gardner.** Funções em um kernel novo (`kernels/timing.py`): erro de tempo, filtro de
   laço (proporcional e integral) e interpolador. Um bloco `ClockRecovery` em
   `components/rx/clock.py`, com a interface de todo bloco, que recebe o `DigitalSignal` a 2 amostras
   por símbolo do ADC e devolve outro, a 2 amostras por símbolo:

   ```
   CoherentFrontEnd -> ADC (2 sps) -> ClockRecovery (2 sps) -> BlindReceiver -> medida
   ```

3. **Desvio de relógio e jitter no ADC**, para o laço ter o que recuperar. `adc.clock_offset_ppm`
   e `adc.jitter_rms` (este já está na lista de `docs/proposta_ideal_real.md`), com valor neutro 0.
   Sem esta etapa o laço só veria um atraso parado, e o `BlindReceiver` já resolve isso.
4. **Mueller-Müller** (opcional), como `algorithm="mueller_muller"` do mesmo bloco. Só faz sentido
   se o ADC passar a 1 amostra por símbolo, ou se o laço for depois do equalizador.

## 5. Como saber se funciona

Insumo para a issue de testes, que fica com o grupo. **[Proposta]**

- **Curva S:** o erro de tempo médio contra um deslocamento imposto passa por zero no deslocamento
  certo e tem o sinal certo ao redor dele. É o gráfico para o notebook.
- **Deslocamento fixo conhecido** (por exemplo 0,3 amostra): o laço o encontra e a EVM fica perto da
  do instante ideal.
- **Desvio de relógio conhecido** (por exemplo 100 ppm): o laço acompanha, e a EVM fica igual à sem
  desvio. Sem o laço, ela cresce ao longo da janela.
- **Sem desvio:** o laço não piora a medida.
- **Determinismo:** sem sorteio, a mesma entrada dá a mesma saída.
- **Erros:** entrada a mais de 2 amostras por símbolo, ou elétrica em vez de digital.

## 6. O que decidir

1. **Escopo.** Só a etapa 1, as etapas 1 a 3 (minha recomendação) ou todas?
2. **Onde na cadeia.** Antes do equalizador (minha recomendação: não mexe no `BlindReceiver`) ou
   depois, como na cadeia do VPI (dispersão, equalizador, Stokes, relógio, frequência, fase). O item
   C4 de `docs/plano_replicar_vpi.md` já trata diferenças de ordem como pequenas.
3. **Interpolador.** O cúbico (Farrow) acompanha um atraso que muda no tempo. A fase linear na FFT
   só serve para um atraso parado. Recomendo o cúbico no laço e a FFT na etapa 1.
4. **D-4.** O interpolador é reamostragem do lado do DSP, não da grade da simulação. Vale como
   exceção, como o skew?
5. **Tamanho da saída com desvio de relógio.** O `DigitalSignal` exige `n_símbolos × amostras por
   símbolo`. Com 100 ppm e 32768 símbolos, o relógio anda cerca de 3 símbolos na janela, e o laço
   precisaria descartar ou repetir amostras. Mantemos o tamanho (cortando ou completando as pontas)
   ou aceitamos um tamanho diferente?
6. **Valores de ppm e jitter.** O VPI não tem. Padrão neutro 0 nos blocos e valores reais só nos
   exemplos em YAML, como na proposta de componentes ideais e reais?
7. **Validação.** Sem o VPI como referência, aceitamos validar só contra verdades impostas (os
   deslocamentos e ppm da seção 5)?

## 7. Esforço

Estimativa minha, ainda não medida: etapa 1, 1 h. Etapa 2, 4 a 5 h (kernel 2 h, bloco 1 h, notebook
1,5 h). Etapa 3, 1 h. Etapa 4, 2 h.
