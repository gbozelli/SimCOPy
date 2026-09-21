# simcopy

Simulador modular de sistemas de comunicações ópticas (WDM, IMDD e coerente digital).

## Para que serve o código que existe hoje

Ainda **não simula um enlace de ponta a ponta**: não há transmissor, receptor nem
BER. O que existe é a fundação e a fibra. Em termos do `draft.py`:

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
pytest                 # 70 testes, ~25 s
```

**Unitários** (`tests/unit/`) perguntam: *o código faz o que diz que faz?*
Não envolvem física. Exemplo: um canal IMDD configurado com 16QAM tem que ser
recusado com mensagem que diga que a causa é `kind = IMDD`.

**De física** (`tests/physics/`) perguntam: *o código faz o que a natureza faz?*
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

Os testes usam `tests/_approx.py` em vez de `pytest.approx` direto. O
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
python scripts/figuras_validacao.py figuras/
```

As figuras não substituem as asserções dos testes; servem para inspeção visual e
documentação. São quatro, uma por teste de física.

## Estrutura

```
src/simcopy/        o pacote (é o que vai para o PyPI)
tests/unit/         testes de software
tests/physics/      testes contra fórmula fechada
scripts/            figuras de validação
examples/           cenários em YAML (o 400ZR mora só aqui)
```

Em layout `src/`, os testes ficam na raiz em `tests/`, não dentro de `src/`:
tudo que está em `src/` é empacotado e distribuído.

## Regra do cenário de referência

Nenhum módulo sob `src/simcopy/` pode conter a string `400ZR` nem os valores do
cenário. Ele vive em `examples/400zr_4ch_75ghz.yaml` e nos testes de regressão.
