# simcopy

Simulador modular de sistemas de comunicações ópticas (WDM, IMDD e coerente digital).

## Para que serve o código que existe hoje

**Atualização (2026-11-12):** o sistema IMDD (OOK, modulação direta ou externa, um ou
vários canais WDM, fibra padrão e compensadora, amplificadores) já roda de ponta a ponta
a partir de um arquivo YAML ou JSON, até o fator Q e a BER de cada canal:

```
simcopy run examples/imdd_ook_mzm.yaml -o resultados/mzm
simcopy run examples/imdd_wdm_4ch.yaml -o resultados/wdm
simcopy sweep examples/imdd_ook_dml.yaml --param "link.spans[0].fiber.length" --range 20e3 40e3 1e3
```

Guia: [`docs/sistema_imdd.md`](docs/sistema_imdd.md).

**Atualização (2026-11-14):** primeira rodada do sistema coerente digital: BPSK e QPSK de
uma polarização, sem fibra, com oscilador local, híbrida de 90° e detecção balanceada, até
a BER contada, a EVM e o fator Q (`simcopy run examples/coerente_qpsk.yaml`). Guia:
[`docs/sistema_coerente.md`](docs/sistema_coerente.md).

**Atualização (2026-11-15):** 16QAM, dupla polarização, vários canais, fibra e DSP com os
símbolos conhecidos (compensação de dispersão, desvio de frequência, polarização, fase e
equalizador). O cenário de 4 canais DP-16QAM do VPI (`examples/400zr_4ch_75ghz.yaml`) foi
comparado com as constelações do VPI em 125 e 140 km: a SNR fica de 0,03 a 0,16 dB acima
([`docs/comparacao_vpi.md`](docs/comparacao_vpi.md)).

**Atualização (2026-11-16):** DSP cego (`rx.dsp.mode: blind`: espaço de Stokes, equalizador
adaptativo, busca cega de fase e ajuste por decisão), ADC com quantização, 64QAM, PAM4 no
IMDD (`examples/imdd_pam4_mzm.yaml`), inclinação da dispersão na fibra, e MI e GMI nas
métricas do coerente. Guias: seção 12 de [`docs/sistema_coerente.md`](docs/sistema_coerente.md)
e seção 11 de [`docs/sistema_imdd.md`](docs/sistema_imdd.md). Com o DSP cego e a janela de
fase de 21, a SNR fica de 0,21 a 0,33 dB acima da do VPI em 125 e 140 km, e a queda de baixa
potência do VPI aparece 1 dB abaixo (seção 9 de [`docs/comparacao_vpi.md`](docs/comparacao_vpi.md)).
Ainda sem recuperação de relógio, código corretor de erros e testes automáticos do que entrou
desde o IMDD. A revisão dos módulos novos está na seção 9 de
[`docs/revisao_2026-11-14.md`](docs/revisao_2026-11-14.md).

**Respostas às dúvidas (2026-11-16):** todas as dúvidas de
[`docs/duvidas_abertas.md`](docs/duvidas_abertas.md) têm resposta, e as que pediam código foram
aplicadas: `vbias = vpi` exigido no coerente, campo ao longo de z na fibra (`Fiber.along`), passo
limitado com PMD, avisos de PMD com campo escalar e de potência acima do limiar de Brillouin no
IMDD, potência do laser em dBm, OOK como formato padrão do IMDD, ruído térmico padrão do Agrawal
(5,76 pA/√Hz) e amplificador controlado por potência medindo só a banda dos canais. O sinal do
desbalanço de 5° do cenário do VPI foi corrigido (−5°, excursão de 0,51) e a comparação refeita em
120, 125, 130, 135 e 140 km: seção 10 de [`docs/comparacao_vpi.md`](docs/comparacao_vpi.md). Os
documentos do pacote (`docs/decisoes.md`, `historico.md`, `roadmap.md`, `arquitetura.md`) vieram
para o repositório. Para discutir: [`docs/proposta_ideal_real.md`](docs/proposta_ideal_real.md).
O exemplo `coerente_dp16qam.yaml` e os notebooks 09 e 10 ainda usam +5° e 0,47.

**Planos (2026-11-16):**
[`docs/plano_replicar_vpi.md`](docs/plano_replicar_vpi.md) (hipóteses para fechar as diferenças
que sobram para o VPI) e [`docs/plano_validacao_cruzada.md`](docs/plano_validacao_cruzada.md)
(conferir bloco a bloco contra outros simuladores e validar os dados). Do primeiro, a etapa A
e os testes 1 e 2 estão feitos ([`docs/replicacao_vpi_diario.md`](docs/replicacao_vpi_diario.md)): com
o skew de uma amostra e a janela de fase de 43, o receptor cego fica a −0,04 a +0,14 dB do VPI,
com média de +0,04 dB. Desde 2026-11-17 o cenário (`examples/400zr_4ch_75ghz.yaml`) usa essa
configuração, e o notebook 09 muda de números quando for reexecutado. O segundo plano ainda não
começou.

**Reprodução do VPI (2026-11-19):** o notebook 09 passou a ser só de 16QAM e da reprodução do VPI,
com o DSP cego (que veio do notebook 10) e a comparação em BER e SER com todos os conjuntos de dados
do grupo, lidos por `scripts/catalogo_vpi.py`. Onde o VPI converge, a BER do simcopy fica de 0,85 a
1,07 vez a dele. Nas bordas o DSP do VPI perde a fase, e as hipóteses para isso vêm uma por seção no
fim do notebook (seção 13 de [`docs/comparacao_vpi.md`](docs/comparacao_vpi.md)). O notebook 10
ficou com ADC, 64QAM, MI e GMI, PAM4 e inclinação da dispersão.

O texto abaixo descreve a fundação, anterior aos componentes. Em termos do `draft.py`:

| Hoje o pacote faz | No `draft.py` isso era |
|---|---|
| Lê um cenário em YAML e recusa combinações inválidas (IMDD com 16QAM, `pol` em IMDD, λ e f incoerentes) | Não existia; qualquer combinação passava |
| Calcula tudo que é derivado: fs, SpS, N, janela, β₂, γ, L_eff, λ | ~66 globais, com `SpS` e `ts` livres ao mesmo tempo e `lambda0` escrito 17 vezes |
| Valida antes de alocar memória (Nyquist, SpS inteiro, tamanho de FFT, sobreposição de canais) | O erro aparecia depois de minutos de SSFM, como `IndexError` |
| Propaga um campo pela fibra: SSFM simétrico, Manakov 8/9, PMD de passo grosso | A função `fiber()`, com os defeitos C-6, C-7 e C-8 |
| Gera streams aleatórios independentes por canal e por papel | `np.random.seed` global: os 4 canais WDM tinham os mesmos bits |

O que falta, nesta ordem: transmissor (bits → forma de onda, com os impairments
do C-25), multiplexador e amplificador, receptor coerente, e aí a BER.

## Testes: dois tipos

```
pytest                 # 73 testes, ~25 s
```

Cada módulo tem a sua pasta de testes: `src/simcopy/<módulo>/tests/`. A ideia é
testar função por função, de baixo para cima (unidades → grade e sinais →
configuração → kernels), para que o comportamento do sistema se apoie em funções
já verificadas.

**Unitários** (`src/simcopy/<módulo>/tests/`) perguntam: *o código faz o que diz que faz?*
Não envolvem física. Exemplo: um canal IMDD configurado com 16QAM tem que ser
recusado com mensagem que diga que a causa é `kind = IMDD`.

**De física** (`src/simcopy/kernels/tests/physical_tests/`) perguntam: *o código faz o que a natureza faz?*
Cada um compara a simulação com uma **fórmula fechada** da literatura, isolando
um efeito por vez. Estão numerados na ordem da lista de testes proposta:

| Arquivo | Efeito isolado | Referência | Tolerância |
|---|---|---|---|
| `test_01_atenuacao.py` | só perda (β₂ = γ = 0) | P(z) = P(0)·e^(−αz) | 1e-12 |
| `test_02_dispersao.py` | só dispersão (α = γ = 0) | FWHM gaussiana, sem e com chirp (Agrawal, NLFO, 3.2) | 2e-3 |
| `test_03_nao_linearidade.py` | só SPM (β₂ = 0) | φ_NL = γ·\|A\|²·L_eff; fase com a forma do pulso; 8/9 em DP (Agrawal, NLFO, 4.1) | 1e-9 |
| `test_04_pmd.py` | só birrefringência | DGD rms = D_p·√z; Maxwelliana; rotação cobre os dois hemisférios (Agrawal, FOCS, cap. 2) | 10 % (estatístico) |
| `test_05_soliton.py` | dispersão **e** SPM juntas | sóliton fundamental N = 1 preserva a forma; SSFM de 2ª ordem (Agrawal, NLFO, 5.2) | 1e-3; ordem > 1,7 |

Os itens "fase deve ter a mesma forma que o pulso" e "não lineares com γ e L_eff"
da lista são o mesmo fenômeno (SPM pura) e estão juntos no teste 03.

O caso **com chirp** do teste 02 foi acrescentado porque o caso sem chirp depende
só de |β₂|: um erro de sinal na dispersão passaria. Com C·β₂ < 0 o pulso comprime
antes de alargar, e só o sinal certo reproduz isso.

### Representações de sinal

| Tipo | O que é | Taxa |
|---|---|---|
| `Bits` | sequência binária | — |
| `Symbols` | símbolos da constelação | 1 por símbolo |
| `DigitalSignal` | sinal na taxa do DSP: entrada do DAC, saída do ADC | `sps_DSP · R_s` (ex.: 2 · 59,84 GBd) |
| `Waveform` | campo na grade global de simulação, elétrico ou óptico | `fs` global (ex.: 16 · 59,84 GBd) |

As quatro compartilham uma `SignalGeometry`. `DigitalSignal` fica fora da grade
global e não tem `f0`: após a detecção coerente o sinal está em banda base, e o
offset residual do LO é conteúdo que o DSP remove, não posição na grade WDM.
Passar de `Waveform` para `DigitalSignal` é dizimação inteira, então o SpS do DSP
tem que dividir o da grade (regra `DSPRateDividesSimulationRate`).

### Tolerâncias

Os testes usam `from simcopy.testing import approx` em vez de `pytest.approx`
direto (nem `np.isclose`, cuja tolerância absoluta padrão é 1e-8). O
`pytest.approx` tem tolerância **absoluta** padrão de 1e-12, e em unidades SI
várias grandezas daqui são dessa ordem (DGD ~1e-13 s, β₂ ~2e-26 s²/m). Com o
padrão, um DGD 3× errado passava.

Pelo mesmo motivo, o gabarito de um teste nunca é uma constante importada do
próprio código: o teste do fator 8/9 escreve `8.0 / 9.0` literalmente. Se usasse
`MANAKOV_FACTOR`, mudar o código mudaria o gabarito junto.

### Revisão por mutação

```
python scripts/mutacoes.py
```

Quebra o código de propósito, um defeito por vez (sinal de β₂ trocado, 8/9
removido, SoP de 2 ângulos, λ fixo em 1550 nm, ...), roda os testes e restaura.
Uma mutação que nenhum teste detecta é uma propriedade sem proteção. Hoje:
10/10 detectadas.

### Testes de contrato

Alguns testes unitários existem para pegar **deriva na integração**: mudanças
que parecem inofensivas mas reintroduzem problemas já corrigidos. O exemplo é
`test_resolve_usa_lambda_derivado_da_frequencia`, que falha se qualquer ponto do
`resolve` voltar a usar λ fixo em vez de λ = c/f. **Rode `pytest` depois de
qualquer integração ou refatoração.** Se um teste falhar, o problema está no
código, não no teste.

## Figuras

```
pip install -e ".[figures]"
python figures/generate_figures/figuras_validacao.py figures/
```

Cada figura tem a mesma estrutura: a **forma fechada** da literatura e a saída do
**nosso SSFM** plotadas juntas, variando um parâmetro, mais um painel com o **erro**
entre as duas. A explicação painel a painel está em [`docs/figuras.md`](docs/figuras.md).

| Figura | Fenômeno isolado | Erro típico |
|---|---|---|
| `01_atenuacao` | só α: varia α, a grade, a forma do pulso e o passo | ~1e-15 |
| `02_dispersao` | só β₂: varia C, D e T₀, com colapso em z/L_D | ~1e-13 |
| `03_nao_linearidade` | só γ: varia φ_max, α (via L_eff) e A_eff (via γ) | ~1e-12 |
| `04_pmd` | só birrefringência: varia D_PMD e l_c | 10 % (estatístico) |
| `05_empilhando` | os fenômenos somados, um a um | ver abaixo |

A figura 05 responde até onde dá para empilhar comparando com o analítico:

| Etapa | Referência | Erro |
|---|---|---|
| α | forma fechada | 1e-14 |
| α + β₂ | forma fechada (gaussiana com chirp) | 1e-13 |
| α + γ | forma fechada (SPM pura) | 1e-12 |
| β₂ + γ | forma fechada **só** para o sóliton fundamental | 1e-4 |
| β₂ + γ + PMD | **não existe forma fechada** | ordem de convergência e conservação de energia |

O salto de 1e-12 para 1e-4 é onde os dois efeitos passam a interagir e aparece o
erro de separação do split-step. É ele que o passo controla, e é por isso que a
figura 05(c) verifica que o erro cai com h².

## Estrutura

```
src/simcopy/                          o pacote
src/simcopy/<módulo>/tests/           testes de cada módulo (core, config, validation, kernels)
src/simcopy/kernels/tests/physical_tests/   testes contra fórmula fechada
src/simcopy/testing/                  utilitários de teste (approx, pulsos, medidas)
figures/                              figuras de validação (PNG)
figures/generate_figures/             scripts que geram as figuras
scripts/                              teste de mutação
examples/                             cenários em YAML (o 400ZR mora só aqui)
```

Os testes ficam junto de cada módulo, como no NumPy e no SciPy
(`numpy/linalg/tests/`). Um teste em `core/tests/` só importa `core`: a regra de
camadas do `import-linter` vale também para os testes.

## Regra do cenário de referência

Nenhum módulo sob `src/simcopy/` pode conter a string `400ZR` nem os valores do
cenário. Ele vive em `examples/400zr_4ch_75ghz.yaml` e nos testes de regressão.
