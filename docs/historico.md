# Histórico: bugs encontrados e lições

Registro do que já deu errado neste projeto e por quê. Existe para que o mesmo
erro não seja cometido de novo. Cada item tem o sintoma, a causa e a proteção
que ficou no repositório.

## Bugs no pacote `simcopy`

### 1. DGD dependia do passo numérico
- **Sintoma:** com passo do SSFM de 5 m e comprimento de correlação de 50 m, o
  DGD rms saía 3,18× maior que o esperado.
- **Causa:** o DGD por passo era `D_PMD·√passo` sempre. Dentro de uma seção de
  birrefringência a orientação não muda, então o DGD cresce linearmente; o
  acúmulo em raiz só vale entre seções.
- **Correção:** `dgd_for_step = D_PMD · passo / √max(passo, L_corr)`.
- **Proteção:** `tests/physics/test_04_pmd.py::test_b` e mutação 4.

### 2. Fase não linear com erro de 5e-4 em passos longos
- **Causa:** o termo não linear usava o comprimento geométrico do passo `dz`
  com a intensidade do ponto médio, ignorando a atenuação dentro do passo.
- **Correção:** comprimento efetivo do passo, `2·sinh(α·dz/2)/α`. Com β₂ = 0 a
  SPM passa a ser exata.
- **Proteção:** `test_03_nao_linearidade.py` com tolerância 1e-9 e mutação 3.

### 3. Último passo ultrapassava o fim da fibra
- **Sintoma:** a curva de convergência do SSFM ficava serrilhada.
- **Causa:** `dz = max(min_step, min(dz, L - z))` forçava o último passo a
  `min_step` mesmo quando faltava menos que isso.
- **Correção:** `dz = min(max(dz, min_step), L - z)`.
- **Proteção:** `test_01_atenuacao.py::test_comprimento_nao_multiplo_do_passo`.

### 4. λ fixo em 1550 nm introduzido na integração
- **Sintoma:** β₂ e γ calculados com λ = 1550 nm independentemente da
  frequência central (a 193,7 THz, λ = 1547,715 nm; erro de 0,3 % em β₂).
- **Causa:** ao integrar o código no repositório via Antigravity/Gemini, o
  agente reescreveu `resolve.py` com `system.center_wavelength if
  hasattr(system, "center_wavelength") else 1550e-9`. O atributo não existia,
  então o `hasattr` era sempre falso.
- **Correção:** `center_wavelength` virou propriedade derivada (`c/f`) de
  `SystemConfig`; o YAML aceita frequência ou comprimento de onda e recusa os
  dois se incoerentes.
- **Proteção:** `tests/unit/test_resolucao.py::test_resolve_usa_lambda_derivado_da_frequencia`
  (falha a 191, 193,7 e 196 THz se λ voltar a ser fixo) e mutação 7.

### 5. Offset de 11 % nas figuras de SPM
- **Causa:** o script de figuras chamava o kernel com `manakov=True` (padrão)
  e comparava com a fórmula escalar, que não tem o 8/9.
- **Lição:** comparações com a NLSE escalar de forma fechada exigem
  `manakov=False` explícito. O painel de erro separado tornou o desvio visível;
  a sobreposição simples teria escondido.

## Bugs nos testes

### 6. `pytest.approx` tornava testes decorativos
- **Causa:** tolerância absoluta padrão de 1e-12. DGD em segundos é ~1e-13 e β₂
  é ~2e-26 s²/m, então qualquer valor passava.
- **Proteção:** `tests/_approx.py`, com tolerância absoluta zero por padrão.
  Sempre usar `from tests._approx import approx`.

### 7. Gabarito importado do próprio código
- **Causa:** o teste do fator 8/9 usava `MANAKOV_FACTOR` do kernel como valor
  esperado. Mudar o código mudava o gabarito junto.
- **Regra:** o valor esperado de um teste vem da física ou da referência,
  escrito literalmente (`8.0 / 9.0`), nunca importado do módulo testado.

### 8. Teste removido em reorganização
- **Causa:** ao reorganizar os testes de física, o teste de uniformidade da
  rotação na esfera de Poincaré foi retirado. Os outros não cobriam o defeito,
  porque com muitas seções até a rotação de 2 ângulos acaba espalhando o estado.
- **Proteção:** `test_04_pmd.py::test_d` e mutação 5.

Os itens 6 a 8 foram encontrados por `scripts/mutacoes.py`: na primeira rodada,
só 7 de 10 mutações eram detectadas. Hoje são 10 de 10.

## Diagnósticos errados já corrigidos

- **F-3 do rascunho:** a versão 0.1 da arquitetura afirmava que o literal `10`
  nos índices de sincronismo deveria ser `N_zeros_init = 100`. Estava errado: o
  `10` compensa a guarda do próprio gabarito de correlação. A causa real dos
  `IndexError` era a ausência de verificação de limites sobre o pico de
  correlação.
- **"O VPI não exercita PMD":** errado. O arquivo tem `PMDCoefficient`,
  `CorrelationLength` e `PolarizationAnalysis = VectorPMD`.
- **Coeficiente de PMD de 0,05 ps/√km:** proposto a partir do BIFROST; o arquivo
  VPI usa 0,1. Vale o arquivo.
- **Regeneração de bits pela seed no receptor:** a versão 0.1 dizia que o
  receptor não deveria regenerar a sequência pela seed. Estava errado. O
  problema do rascunho era o `np.random.seed` global, não a regeneração.

## Lições sobre o código legado (`legacy/`)

- **BER dependente de SpS (P-15):** o fotodiodo do rascunho usa filtro elíptico
  com `gpass = 3 dB` aplicado com `lfilter`; a ordem do filtro depende de `fs` e
  a fase não linear introduz ISI. Trocar SpS de 16 para 4 mudava a BER em três
  ordens de grandeza. O código só funcionava porque a banda elétrica estava em
  80 GHz, quase 6× a taxa de símbolo. Por isso o legado não serve de gabarito.
- **Canais WDM idênticos:** os 4 transmissores da célula WDM do notebook usavam
  as mesmas seeds, e a seed também gerava o payload e reiniciava o RNG global.
  Os quatro canais carregavam a mesma sequência de bits e o mesmo ruído de fase.
- **BER "boa demais":** no cenário de 1 span de 120 km a margem de OSNR era de
  12 a 24 dB; BER zero era o resultado correto. O rascunho também não tinha o
  piso de implementação do TX (`B`), simulava 14 GBd em vez de 59,84 GBd e
  lançava 2 mW.
- **Custo:** 2 milhões de símbolos não são necessários. Para BER de 1e-3 com 10 %
  de erro relativo bastam ~3.300 símbolos; o VPI usa 262.144.

## Lições de processo

- **Agentes de integração reescrevem código.** Ao passar arquivos para outro
  agente, exigir cópia literal e rodar `pytest` e `scripts/mutacoes.py` depois.
- **Agentes consertam o teste em vez do código.** Nunca editar um teste para ele
  passar. Se um teste falha, o problema está no código até prova em contrário.
- **Revisar código que não se escreveu:** a pergunta útil não é "isto está
  certo?", é "se isto estivesse errado, algum teste avisaria?". Para isso existe
  `scripts/mutacoes.py`.
