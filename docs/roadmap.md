# Roadmap

> Substituído por `docs/proximos_passos.md` (2026-11-19). Fica como registro do plano original.

Tarefas em ordem de dependência. Cada uma tem critério de aceitação verificável.
Os parâmetros de referência vêm de `reference/vpi/parametros_vtmu.md` (arquivo
VPI 11.5 `WDM_of_400G_ZR_4_Channels_Ch2-_versao_editavel.vtmu`).

Regra geral para todas as tarefas: `pytest` passa e `python scripts/mutacoes.py`
continua com todas as mutações detectadas. Toda capacidade nova ganha teste
contra forma fechada quando ela existe, e uma mutação em `scripts/mutacoes.py`
que o teste precisa pegar.

---

## 0. Reconciliar o repositório com este pacote

O código deste pacote foi integrado antes ao repositório por outro agente
(Antigravity/Gemini), que reescreveu trechos. Um deles introduziu λ fixo em
1550 nm (`docs/historico.md`, item 4).

- Comparar o repositório com este pacote arquivo a arquivo e listar as
  diferenças para o Gabriel **antes** de sobrescrever qualquer coisa. Pode haver
  mudanças dele que precisam ser mantidas.
- Aceitação: `pytest` com 73 testes passando e `scripts/mutacoes.py` com 10/10.

## 1. Remover o cenário 400ZR de `src/`

Os valores padrão em `src/simcopy/config/system.py` e `config/io.py` são os do
400ZR (largura de linha 500 kHz, roll-off 0,18, ruído térmico 10 pA/√Hz, banda
1,15·R_s, 70 GHz, S21 40 GHz, `center_frequency = 193,7 THz`). Isso viola a
regra de que o cenário mora só em `examples/` e nos testes.

- **Perguntar ao Gabriel** entre: campos obrigatórios (sem padrão) ou padrões
  genéricos documentados. Não decidir sozinho.
- Trocar também a docstring de `cli.py`, que cita o arquivo do 400ZR.
- Aceitação: `grep -rn "400ZR\|400zr\|59.84\|193.7" src/` vazio; o YAML de
  referência continua carregando e validando.

## 2. Inclinação de dispersão (β₃)

A fibra do VPI tem `DispersionSlope = 0.08e3 s/m³` = 0,08 ps/(nm²·km). O kernel
só tem β₂. Com 4 canais ocupando 1,8 nm, D varia ~1 % na banda, e a compensação
de CD do canal sob teste depende disso.

- Acrescentar `dispersion_slope` a `FiberConfig`, β₃ a `ResolvedFiber` e o termo
  `(β₃/6)·ω³` ao operador linear de `kernels/ssfm.py`.
- Teste de física contra a forma fechada da largura rms de gaussiana com β₂ e β₃
  (Agrawal, *Nonlinear Fiber Optics*, seção 3.3):
  `σ²/σ₀² = (1 + Cβ₂z/2σ₀²)² + (β₂z/2σ₀²)² + (1 + C²)²·(β₃z/4√2σ₀³)²`.
- Acrescentar ao teste 02 um caso com β₂ = 0 e só β₃ (pulso assimétrico, com
  oscilação na cauda cujo lado depende do sinal de β₃).
- Mutação: sinal de β₃ invertido.
- Observação: o bloco de compensação de CD do DSP no arquivo VPI tem
  `DispersionSlope = 0.08e-3`, seis ordens de grandeza abaixo do da fibra.
  Provável erro de digitação no esquemático; perguntar ao Gabriel.

## 3. Transmissor (`src/simcopy/components/tx/`)

Interface de todos os blocos (seção 5 da arquitetura):
`__call__(self, x, ctx) -> saída`, sem estado global, sem plot, sem
`np.random` global. Aleatoriedade por `ctx.rng_for(...)`.

| Bloco | Entrada → saída | Parâmetros do VPI |
|---|---|---|
| `BitSource` | — → `Bits` | seed do canal. O receptor regenera os bits pela mesma seed |
| `SymbolMapper` | `Bits` → `Symbols` | 16QAM; **suportar com e sem Gray** (`GrayMapping = No` no VPI); vetorizado |
| `PulseShaper` | `Symbols` → `DigitalSignal` | RRC, roll-off 0,18 |
| `DAC` | `DigitalSignal` → `Waveform` elétrico | interpolação para a grade global; LPF do TX |
| `IQModulator` | `Waveform` elétrico → óptico | laser 1 mW e 500 kHz; `Vpi = 5`, IL 6 dB, ER 30 dB, S21 40 GHz ordem 3; desbalanço de fase 5°; skew 0,75 ps; `PreDistortion = Yes` (confirmar com o Gabriel o que isso faz no VPI) |

Aceitação:
- RRC seguido de filtro casado tem ISI nula nos instantes de símbolo (erro < 1e-6).
- `SymbolMapper` seguido do desmapeamento devolve os bits exatos, com e sem Gray.
- ER medida na saída do MZM = ER configurada (± 0,1 dB).
- Potência média da saída bate com laser × perdas.
- Desbalanço IQ sem ruído produz o EVM previsto pela forma fechada do desbalanço
  de fase (conferir a fórmula e deixá-la escrita no teste).

## 4. Óptica (`src/simcopy/components/optical/`)

| Bloco | Parâmetros do VPI |
|---|---|
| `Multiplexer` | 4 canais, espaçamento 75 GHz, gaussiana de ordem 3, BW 70 GHz, IL 3 dB |
| `Amplifier` | modo potência de saída alvo (`OutputPower = PowerdBm`), NF 5 dB; ASE gaussiana circular `sqrt(P_ASE/2)·(randn + j·randn)` |
| `Fiber` | envelope fino de `kernels/ssfm.py`; RNG da birrefringência por `ctx.rng_for("link", span, "birefringence")`; `manakov=True` sempre (D-2) |
| `Demultiplexer` | gaussiana de ordem 3, BW 70 GHz, centrada no canal |

Aceitação:
- Mux seguido de demux ideal devolve o canal com erro abaixo de uma tolerância
  declarada.
- OSNR medida após o amplificador = OSNR prevista por `NF·h·ν·(G−1)·B_ref`
  (± 0,2 dB).
- O amplificador em modo potência alvo entrega exatamente a potência pedida.
- **D-22 em aberto:** o amplificador precisa aceitar a potência como total do
  campo agregado ou por canal, e o spike testa as duas.

## 5. Receptor coerente (`src/simcopy/components/rx/`)

Parâmetros do VPI (`Rx_DualPol`, `ModelType = Ideal`): LO 1 mW, 500 kHz, offset
500 MHz; híbrida de 90°; fotodiodos balanceados, R = 1,0 A/W, ruído térmico
10 pA/√Hz, shot ligado; LPF **Bessel** de ordem 4 com banda 1,15·R_s; ADC a
2 amostras/símbolo e 12 bits, saída `DigitalSignal`.

Aceitação:
- Filtro com fase linear na banda do sinal (C-21): variação do atraso de grupo
  na banda ocupada < 0,01 T_símbolo.
- Sem ruído e com LO ideal, o front-end devolve o campo do canal (a menos de
  fator de escala) com erro abaixo de uma tolerância declarada.
- Variância do ruído na saída = orçamento de shot + térmico
  (seção 14.16 da arquitetura), ± 5 %.
- O resultado não depende de `SpS` da grade (o defeito P-15 do código legado).

## 6. DSP data-aided e BER

Decidido: no spike, o receptor usa os símbolos conhecidos. Motivo: na janela
convergida o DSP cego do VPI está perto do ideal, então data-aided mede física
sem confundir com convergência de algoritmo.

- Compensação de CD com `L` e `D` verdadeiros (D-7), incluindo β₃ (tarefa 2).
- Equalizador 2×2 por mínimos quadrados com os símbolos conhecidos.
- Recuperação de fase data-aided.
- Decisão, desmapeamento, BER por polarização.
- Aceitação: em AWGN puro, BER contra SNR bate com a curva teórica de 16QAM
  (± 0,2 dB em SNR para BER entre 1e-4 e 1e-2), com e sem Gray.

## 7. Pipeline e execução

- `pipelines/`: `build_tx`, `build_link`, `build_rx`, montados a partir da
  configuração; a escolha IMDD/DCS acontece no builder, nunca dentro de um bloco.
- Runner: valida, resolve, aloca, executa, mede.
- CLI: `simcopy run <config.yaml> -o resultados/`, gravando `.json` com a
  configuração completa, a versão do pacote e o hash do commit, e `.npz` com os
  dados de plotagem.
- Varredura de potência e de distância.

## 8. Spike: reproduzir A, B e C

Alvos (`reference/vpi/ajuste_abc.py`, ajuste conjunto de 7 distâncias, resíduo
RMS 0,083 dB):

| Alvo | Valor | Bloco responsável | Tolerância |
|---|---|---|---|
| `A125` | 0,1186 | front-end do RX | ± 0,5 dB |
| `B` | 0,01224 (piso de 19,12 dB) | impairments do TX (C-25) | ± 0,5 dB |
| `C125` | 5,671e-6 | fibra: γ, L_eff, XPM | ± 1 dB na penalidade a 17 dBm |
| Potência ótima | por distância, ver saída do script | composição | ± 0,5 dB |

- Rodar as duas hipóteses de D-22 e reportar qual fecha.
- Comparar só dentro da janela convergida de cada distância (a janela sobe de
  0–17 dBm a 80 km para 12–16 dBm a 140 km).
- Uma discrepância em `A` aponta para o front-end; em `B`, para o transmissor;
  em `C`, para a fibra. D-21 (FWM) pode explicar parte de uma discrepância em `C`.
- 262.144 símbolos por polarização, como no VPI.

## 9. DSP cego

Marco próprio, validado contra o penhasco de baixa potência dos dados do VPI
(abaixo do limiar a BER trava em ~0,17).

Cadeia do VPI: compensação de CD → MIMO adaptativo MMA (7 taps, 20 iterações,
passo 5e-4) → demultiplexação de polarização no espaço de Stokes → recuperação
de relógio → compensação de offset de frequência → Blind Phase Search
(64 ângulos de teste).

## 10. Desempenho

Medido no código legado (seção 15 da arquitetura). Aplicar no pacote novo:
correlação por FFT em vez de `np.correlate`, mapeamento e desmapeamento
vetorizados, recuperação de fase por broadcasting, `scipy.fft` com `workers=-1`.

## Depois

- IMDD: OOK e PAM4, modulação direta (equação de taxa do laser) e externa,
  detecção direta, shot noise instantâneo (C-9).
- Compensação em linha: bloco de DCF por span.
- Publicação: `CITATION.cff`, `CHANGELOG.md`, `__all__` explícito, documentação
  com mkdocs, CI rodando `pytest`, `mutacoes.py` e `import-linter`.
