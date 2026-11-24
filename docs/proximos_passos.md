# Próximos passos (2026-11-19)

Contexto. O professor comentou: incluir no DSP o MIMO, a ortogonalização de Gram-Schmidt e a
compensação de PMD, e não é preciso simular o VPI. Com isso o projeto entra em fechamento: o
Gabriel e o grupo escrevem os testes unitários, a documentação e o resto das implementações, pelas
issues. Este documento substitui `docs/roadmap.md` e `docs/planejamento.md`.

As fontes de cada implementação e de cada decisão estão em `docs/referencias.md` (as siglas entre
colchetes abaixo são as de lá). As decisões ficam em `docs/decisoes.md` (D-26: replicação do VPI
encerrada, D-27: DSP do professor, D-28: fechamento).

## 1. DSP pedido pelo professor (D-27)

A cadeia do receptor coerente, depois dos três itens, fica nesta ordem ([Sa10], [FS17]):

    ADC -> GSOP -> compensação de dispersão -> MIMO 2x2 adaptativo (polarização e PMD)
        -> desvio de frequência -> busca de fase -> decisão

O que já existe em `kernels/dsp.py`: separação de polarização no espaço de Stokes ([Sz10]) e um
equalizador adaptativo 2×2 de 7 coeficientes a 2 amostras por símbolo (módulo constante, [Go80]),
seguido de um ajuste dirigido por decisão em bloco, por mínimos quadrados. Ele já é um MIMO, mas
não é um bloco próprio, não atualiza símbolo a símbolo e nunca foi testado com PMD.

### 1.1 MIMO 2×2 adaptativo como bloco (#32)

- O que fazer: um bloco `MIMOEqualizer` em `components/rx/dsp.py`, com a borboleta de quatro
  filtros FIR ([Sa08], [IK07]). Algoritmos: módulo constante ([Go80]), dirigido por raio para 16QAM
  e 64QAM ([Fa09]) e LMS dirigido por decisão. Número de coeficientes, passo e passagens vêm do
  `DSPConfig`.
- Aceitação: (a) sem ruído, desfaz uma matriz de Jones conhecida (rotação e atraso), com erro abaixo
  de uma tolerância declarada no teste. (b) Em ruído branco, a BER fica a menos de 0,2 dB de SNR da
  do receptor com símbolos conhecidos. (c) Com um coeficiente, equivale a uma matriz 2×2 constante.
- Material de partida: `scripts/experimentos/equalizador_gradiente.py` tem protótipos de
  equalizador por gradiente, símbolo a símbolo (hipótese 2). São experimentos: revisar antes de usar.

### 1.2 Ortogonalização de Gram-Schmidt (GSOP) (#33)

- O que fazer: corrigir o desbalanço de fase e de amplitude entre I e Q do receptor (híbrida e
  fotodiodos), logo depois do ADC ([Fa08]). Com I e Q de média nula:
  I' = I/√P_I, Q'' = Q − ρ·I/P_I com ρ = E[I·Q], e Q' = Q''/√P_Q''.
- Pré-requisito: hoje o `CoherentReceiver` não tem desbalanço IQ. É preciso acrescentar um modelo
  de front-end imperfeito: erro de fase da híbrida (graus) e diferença de responsividade entre os
  fotodiodos de I e de Q. Sem isso a GSOP não tem o que corrigir.
- Aceitação: com erro de fase e desbalanço de ganho conhecidos, a saída tem correlação I·Q nula e
  potências iguais (tolerância declarada), e a EVM volta ao valor sem desbalanço.
- Atenção (inferência, conferir em [FS17]): o desbalanço do transmissor (os −5° do 400ZR) passa pela
  dispersão antes de chegar ao receptor, e a GSOP do receptor não o corrige. Para ele seria preciso
  um equalizador que trate I e Q separados (linear amplo).

### 1.3 Compensação de PMD (#34, depende da #32)

- O que fazer: a PMD é compensada pelo próprio MIMO, que inverte a matriz de Jones do enlace
  dependente da frequência, desde que tenha coeficientes suficientes para o DGD ([IK07], [Sa08]). O
  trabalho é verificar, dimensionar o número de coeficientes e documentar.
- Teste: com o padrão de 0,1 ps/√km (D-17), o DGD médio em 140 km é 0,1·√140 ≈ 1,2 ps, ou 0,07
  símbolo a 59,84 GBd: quase nada. Para ver o efeito, acrescentar um caso com DGD de primeira ordem
  fixo (por exemplo 10 e 30 ps) além do caso estatístico ([GK00]).
- Aceitação: uma curva de penalidade de BER contra DGD, com o MIMO ligado e só com o espaço de
  Stokes. Com o MIMO, a penalidade fica abaixo de um limite declarado até o DGD que o número de
  coeficientes cobre.

## 2. Testes unitários (Gabriel e grupo)

Cada issue de teste segue o modelo de `memory/issues/v2/_comum.md`: valor esperado escrito como
número, com a conta e as unidades; `simcopy.testing.approx` com tolerância explícita; uma mutação em
`scripts/mutacoes.py` que o teste detecta; nunca editar um teste para ele passar. O `pytest` passa
hoje com 73 testes. As mutações não foram rodadas de novo nesta etapa.

Já abertas: #5 a #10 (núcleo, configuração, fibra), #14 (laser), #15 (MZM), #19 (amplificador),
#20 (fotodetector), #21 (WDM). Novas:

| Issue | Módulos |
|---|---|
| #25 mapeamento e formatação de pulso | `kernels/mapping.py`, `kernels/pulse_shaping.py`, `components/tx/iq_driver.py` |
| #26 modulador IQ | `components/tx/iq_modulator.py`, `iq_field` e `arcsin_predistortion` de `kernels/modulator/mzm.py` |
| #27 receptor coerente e ADC | `kernels/coherent.py`, `components/rx/coherent.py`, `kernels/adc.py` |
| #28 DSP | `kernels/dsp.py`, `components/rx/dsp.py`, `components/rx/adc.py` |
| #29 medida | `kernels/demodulation.py`, `kernels/information.py`, `kernels/decision.py` |
| #30 transmissor IMDD | `kernels/electrical.py`, `kernels/filters.py`, `components/tx/driver.py`, `components/tx/bits.py`, `components/tx/modulator.py` |
| #31 sistema | `pipelines/imdd.py`, `pipelines/coherent.py`, `pipelines/runner.py` (de ponta a ponta, contra as curvas teóricas) |

## 3. Documentação

- Conferir e validar `docs/referencias.md` e `docs/decisoes.md`, em especial o que está marcado
  como inferência (seção 4 das referências).
- Revisar os notebooks 07 a 10 (#35 a #38), como nas issues #11 a #13 e #16 a #18 dos anteriores.
- Validar as fontes (#41).
- Responder a dúvida 71 (perguntas do notebook 08: símbolos sem ADC, SNR teórica, controle de
  potência).
- README: incluir o sistema coerente e apontar para os documentos acima.

## 4. Decisões pendentes (#40)

| Item | Onde |
|---|---|
| Dúvidas sem resposta: 59, 65, 66, 67, 68, 70 e 71 | `docs/duvidas_abertas.md` |
| Respostas interpretadas, a confirmar: 4, 9, 17 e 21 | `docs/duvidas_abertas.md` |
| Valores do 400ZR que sobram em `src/`: `center_frequency = 193,7 THz` e a docstring da CLI | `docs/decisoes.md`, abertas |
| Propostas em discussão: figuras, ideal e real, potência de lançamento, recuperação de portadora, recuperação de relógio, API | `docs/proposta_*.md` |

## 5. Figuras (#39)

- Notebooks 01 a 06: ficaram na versão anterior ao estilo novo (D-25), porque as figuras saíram com
  a legenda em cima e mal formatadas. Aplicar o estilo neles figura por figura, conferindo cada uma.
- Notebooks 07 a 10: as figuras em SVG ficam mais largas que a tela e pedem rolagem lateral. Limitar
  a largura (`figsize`) ou mostrar PNG.
- Guardar na memória do projeto o que já foi pedido: comparação com uma referência num painel por
  caso, referência como curva, sem marcas que poluem.

## 6. Integração e publicação

- O branch `laser` já entrou no `main` (PRs #22 e #24). O trabalho desta etapa vai num branch novo, com o seu PR.
- CI com `pytest`; mutações em série (`scripts/mutacoes.py`).
- `CITATION.cff`, `CHANGELOG.md`, `__all__` explícito, documentação com mkdocs (D-12).

## 7. Desempenho (depois)

Um ponto do 400ZR com 2^14 símbolos leva de 15 s (potência média) a 15 minutos (20 dBm), porque o
passo da fibra encurta com a não linearidade. Começar por um perfil (passos da fibra, FFT,
equalizador, busca de fase). Opções depois do perfil: FFT na GPU, o laço do equalizador em código
compilado (Numba ou Rust), pontos em paralelo sem disputar a CPU.

## 8. Replicação do VPI (encerrada, D-26)

Estado final: no meio da faixa a BER do simcopy fica de 0,85 a 1,07 vez a do VPI, e o ruído que
chega à decisão tem a mesma composição. Nas bordas, o modelo da hipótese 4 (fase conjunta,
estimativa grossa sem distinguir 0° de 45°) reproduz a taxa de perda de fase contra a SNR, o modo de
falha e a BER do meio, mas não o estado de 45° limpo. Tudo no diário e no notebook 09.

Se for retomada: pedir ao grupo o `256QAM_constellation.txt` e o trecho do manual do
`BlindPhaseSearch`; testar uma busca fina que não distinga 0° de 45°; simular pontos perto do limiar
com 2^18 símbolos; testar os conjuntos de configuração desconhecida (hipótese 6).

## 9. Fora do escopo por enquanto

Retropropagação digital, séries de Volterra, FEC, 32QAM e recuperação de relógio
(`docs/proposta_recuperacao_de_relogio.md`).

## 10. Issues abertas nesta etapa

Abertas em 2026-11-19 em gbozelli/SimCOPy, a partir dos rascunhos de `memory/issues/v3/` (fora do
git). As de teste seguem o modelo da versão 2 (`memory/issues/v2/_comum.md`).

| Issue | Assunto |
|---|---|
| #25 a #31 | testes: mapeamento e pulsos, modulador IQ, receptor coerente e ADC, DSP, medida, transmissor IMDD, sistema |
| #32, #33, #34 | DSP do professor: MIMO, GSOP, compensação de PMD |
| #35 a #38 | revisão dos notebooks 07, 08, 09 e 10 |
| #39 | figuras: estilo nos notebooks 01 a 06 e largura das figuras |
| #40 | decisões pendentes |
| #41 | validar as fontes (`docs/referencias.md`) |
