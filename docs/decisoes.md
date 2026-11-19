# Decisões de projeto

Fonte única de verdade sobre o que está decidido. Quando este arquivo e
`docs/arquitetura.md` divergirem, vale este. Os identificadores `[D-n]` são os
mesmos usados no documento de arquitetura e nos comentários do código.

Legenda: **Fechada** (decidida, implementar assim) · **Aberta** (não decidir
sozinho: perguntar ao Gabriel) · **Verificar** (depende de conferir o arquivo VPI).

De onde vem cada decisão e cada implementação (livro, artigo, código do professor, VPI ou medida
nossa): `docs/referencias.md`. O que falta fazer: `docs/proximos_passos.md`. As respostas de cada
dúvida, com o raciocínio: `docs/duvidas_abertas.md`.

## Fechadas

| ID | Decisão | Consequência no código |
|---|---|---|
| D-1 | `QAM_mod` existe no notebook original (célula 5); se perdeu na exportação para `draft.py` | Reescrever vetorizado, sem laço por símbolo |
| D-2 | Propagação SP também usa Manakov (fator 8/9) | No pipeline, `FiberKernelParams(manakov=True)` sempre. `manakov=False` só em testes que comparam com a NLSE escalar de forma fechada |
| D-3 | O fator 8/9 entra | Fechada pelo arquivo VPI (`NonlinearAdjustmentFactors = (8./9.) ...`). O Gabriel não respondeu explicitamente; se algum resultado antigo do grupo, calculado com γ puro, precisar ser reproduzido, perguntar |
| D-4 | `fs` é múltiplo inteiro de todos os baud rates | Sem reamostragem fracionária na v1. Validador `IntegerSamplesPerSymbol` |
| D-5 | Sinal elétrico real é guardado como complexo com parte imaginária nula | Um único tipo `Waveform` para elétrico e óptico |
| D-6 | Margem de alargamento espectral no validador de Nyquist é constante configurável | `SimulationConfig.broadening_margin = 1.2` |
| D-7 | O receptor recebe `L` e `D` verdadeiros do enlace | Compensação de CD usa os valores do enlace (oracle). Sem estimação cega de dispersão na v1 |
| D-8 | PMD de 1ª ordem entra nesta versão | Implementada em `kernels/birefringence.py` |
| D-9 | MZM tem razão de extinção como parâmetro; o resto dos dispositivos mantém o nível do rascunho | Sem RIN, corrente de escuro ou TIA detalhado. **Atenção:** os impairments de TX do arquivo VPI (desbalanço IQ, skew, banda S21, ADC de 12 bits) são necessários para reproduzir o termo `B` (C-25); não contam como aumento de fidelidade, são parte da referência |
| D-10 | Sem GPU por enquanto | NumPy/SciPy puros. Não escrever abstração de backend |
| D-11 | Pacote instalável | Layout `src/`, `pyproject.toml` |
| D-12 | Será publicado | Licença permissiva (MIT no `pyproject`), SemVer começando em `0.x`, API pública explícita, `CITATION.cff` e DOI no futuro |
| D-13 | Tipos de sinal da seção 4.2.1 da arquitetura, mais uma quarta representação | `Bits`, `Symbols`, `Waveform` (grade global) e `DigitalSignal` (taxa do DSP, entrada do DAC e saída do ADC, sem `f0`), todos ligados por `SignalGeometry`. `domain` é metadado para pré-condição e visualização, nunca para escolher física |
| D-14 | Canais SP e DP podem coexistir no mesmo enlace | Basta um canal DP para o solver ser vetorial: `SystemConfig.link_is_vectorial` |
| D-15 | IMDD é sempre SP | `IMDDChannel` não tem campo `pol` |
| D-16 | O pacote se chama `simcopy` | Imports `from simcopy...` |
| D-17 | Coeficiente de PMD padrão: 0,1 ps/√km | Valor do arquivo VPI (`PMDCoefficient = 0.1e-12/31.62`). O valor de 0,05 sugerido a partir do BIFROST foi descartado por afastar da referência |
| D-18 | `correlation_length` é o comprimento de seção do modelo convencional de segmentos (não o período de spin); valor 50 m | Valor do arquivo VPI (`CorrelationLength = 50.0`) |
| D-19 | PMD por estatística calibrada; sem modelo por mecanismo físico (temperatura, curvatura, torção) | BIFROST serve apenas offline, como fonte de números. **É GPL-3.0: nunca importar**, senão o pacote inteiro vira GPL |
| D-20 | `f0` (posição na grade WDM) fica no nível comum do canal; `lo_frequency_offset` fica dentro de `CoherentReceiver` | O multiplexador precisa ver a grade inteira |
| — | `TransmitterConfig` e `ReceiverConfig` ficam abaixo da decisão IMDD/DCS | Uniões discriminadas em `config/system.py`; no YAML, `tx` e `rx` dentro de cada canal |
| — | Span pode ter ou não compensação em linha | `SpanConfig.line_compensation`. O bloco de DCF ainda não existe |
| — | BER é a métrica de comparação | O VPI só exporta BER com facilidade. EVM pode ser usada como diagnóstico interno rápido, nunca como base de comparação com a referência |
| — | A referência de correção é o VPI TransmissionMaker, não o código legado | O código legado tem BER que depende de `SpS` (seção 15.3 da arquitetura) e não serve de gabarito |
| — | Plotagem fica fora da física | Blocos não importam matplotlib. Dados de plotagem são exportados (`.npz`/`.json`) e consumidos por scripts separados |
| — | Receptor data-aided no spike; DSP cego depois, como marco próprio | Na janela convergida o DSP do VPI está perto do ideal, então data-aided mede física sem confundir com convergência de algoritmo |
| — | 400ZR é caso de teste, não forma do código | Nenhum módulo em `src/` pode conter a string `400ZR` nem os valores do cenário. **Hoje o código viola isso** nos valores padrão; ver `docs/roadmap.md`, tarefa 1 |
| — | Validação física: analítico contra implementação, um fenômeno por vez, depois empilhando | Figuras com a forma fechada e o SSFM juntos e o erro em painel próprio (`docs/figuras.md`) |
| D-21 | FWM: aceitar a diferença (o VPI tem `FWM = No`, o SSFM de campo agregado inclui FWM) | Fechada em 2026-11-16 (dúvida 36). Registrado em `docs/comparacao_vpi.md`, seção 6 |
| D-22 | `PowerdBm` do VPI é a potência total dos 4 canais | Fechada em 2026-11-16 (dúvida 61). O YAML do 400ZR usa 13 dBm no booster e `scripts/compara_vpi.py` faz a conta por canal |
| D-23 | O DSP padrão do cenário 400ZR é o receptor cego: espaço de Stokes, módulo constante, ajuste dirigido por decisão, desvio de frequência e busca cega de fase conjunta com janela de 43 | `DSPConfig(mode="blind", phase_window=43)` no YAML do 400ZR. O modo com símbolos conhecidos continua disponível. A ordem e os padrões estão em aberto (dúvida 66) |
| D-24 | O YAML do 400ZR reproduz o transmissor do VPI medido nas constelações dele: desbalanço de fase de −5°, `drive_depth` 0,51 e skew de uma amostra (1,044 ps) | `examples/400zr_4ch_75ghz.yaml`. O skew de uma amostra é inferência (o VPI arredondaria o atraso para amostras inteiras): confirmar com o grupo |
| D-25 | Figuras no estilo do guia do TCC: PDF fora do notebook, SVG dentro; o estilo novo vale do notebook 07 em diante | `plotting.apply_style` e `save_figure`. Os notebooks 01 a 06 ficaram na versão anterior até consertar as legendas e a largura (ver `docs/proximos_passos.md`) |
| D-26 | A replicação do VPI está encerrada (professor, 2026-11-19: "não precisa simular o VPI") | Os experimentos ficam em `scripts/experimentos/`, fora de `src/`. O resumo e as hipóteses ficam no diário e no notebook 09 |
| D-27 | O DSP coerente ganha MIMO 2×2 adaptativo como bloco, ortogonalização de Gram-Schmidt (GSOP) e compensação de PMD (professor, 2026-11-19) | A implementar (`docs/proximos_passos.md`, seção 1). O receptor ainda não tem parâmetros de desbalanço IQ para testar a GSOP |
| D-28 | Fechamento do projeto: o Gabriel e o grupo escrevem os testes, a documentação e o resto das implementações, pelas issues | Cada issue de teste segue o modelo de `memory/issues/v2/_comum.md`. As fontes de cada implementação estão em `docs/referencias.md` |

## Abertas

| ID | Questão | Estado da evidência | Impacto |
|---|---|---|---|
| — | Valores padrão da configuração | A dúvida 39 resolveu atenuação, dispersão, área efetiva, n₂, responsividade e ruído térmico. Ainda sobram valores do 400ZR em `src/`: `center_frequency = 193,7 THz` (`config/system.py` e `config/io.py`) e a docstring de `cli.py` | Opções: campos obrigatórios (sem padrão) ou padrões genéricos. Decisão do Gabriel |
| — | Dúvidas sem resposta: 59 (perda do modulador IQ), 65 (fundo de escala do ADC), 66 (padrões e ordem do receptor cego), 67 (PAM4), 68 (β₃ por padrão), 70 (gráfico de símbolos), 71 (notebook 08) | `docs/duvidas_abertas.md` | Listadas em `docs/proximos_passos.md`, seção 4 |
| — | Respostas interpretadas, a confirmar: 4, 9, 17 e 21 | `docs/duvidas_abertas.md` | A 9 ("sim" a G ≥ 0) conflita com a 22 ("sempre o Agrawal"): não aplicada |

## A verificar no arquivo VPI (só se a replicação for retomada: D-26)

| Item | Por que importa |
|---|---|
| `GrayMapping = No` | A conversão BER → SNR do ajuste assume Gray; os SNR absolutos extraídos ficam deslocados. O mapeador precisa suportar os dois |
| `PRBS_Order = 7` | Período de 127 bits repetido ~8.258 vezes na janela: cria correlação artificial entre símbolos e entre canais, o que afeta XPM |
| Exportação do campo óptico | Permitiria o nível de validação V3 (forma de onda após a fibra). Hoje só BER é exportável |
| Se o amplificador é só booster | Os dados indicam que sim (o termo `A` escala com a perda do span); a documentação do resumo sugeria o contrário |
