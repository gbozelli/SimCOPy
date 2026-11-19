# Figuras para relatório, slides e artigo, e em inglês: proposta para discussão (2026-11-18)

Nasceu do pedido de 2026-11-18: as figuras seguem o guia do TCC
(`TCC/codigos/guia_figuras.md`) no tamanho dos notebooks, e fica para o futuro gerar cada figura
no tamanho de cada documento e em outro idioma.

**Proposta em uma frase:** o estilo ganha dois parâmetros, o destino (notebook, relatório, slides,
artigo) e o idioma (português ou inglês), e cada função de figura passa a ler os dois em vez de ter
tamanho, fonte e textos fixos.

Nada aqui foi implementado. A seção 4 lista o que precisa ser decidido.

---

## 1. Como está hoje

- `src/simcopy/plotting.py` segue o guia: fonte Times, ticks para dentro, grade tracejada, sem
  título, rótulo de painel "(a) ...", legenda acima, cores nomeadas, vírgula decimal e constelações
  em mapa de densidade. As medidas de fonte e linha são as do relatório (fonte 11, ticks 10, linha
  1,6), e o tamanho é o dos notebooks (cerca de 10 polegadas de largura).
- O tamanho de cada figura está escrito dentro de cada função, e os textos também, em português.
- As figuras gravadas saem em PDF (`save_figure`, na CLI e nos scripts), com um PNG ao lado só
  quando um documento em Markdown precisa mostrá-las.
- Os notebooks chamam `apply_style()` na primeira célula (desde 2026-11-18), então as figuras
  feitas direto neles também seguem o estilo, menos o que o estilo não alcança: títulos, posição
  da legenda e rótulos de painel continuam como cada célula escreve.
- Os scripts de comparação (`scripts/compara_vpi.py` e `scripts/ber_vpi_simcopy.py`) usam o mesmo
  estilo.
- As figuras feitas direto nos notebooks não usam o estilo, a não ser que o notebook chame
  `apply_style()`.

## 2. Destinos

O guia gera uma versão de cada figura por documento, já no tamanho final, porque reescalar no LaTeX
muda o tamanho da fonte (guia, seção 1):

| Destino | Largura | Fonte base | Linhas |
|---|---|---|---|
| notebook (hoje) | cerca de 10 pol | 11 pt | 1,6 |
| relatório ou TCC | 6,3 pol | 11 pt | 1,6 |
| slides | até 5,3 pol | 9 pt | 1,1 |
| artigo em duas colunas | 3,39 pol (coluna) ou 7,08 pol (página) | 8 pt | 1,1 |

Proposta:

- O dicionário de estilo de hoje vira uma função do destino, com as medidas da tabela.
- Cada função de figura recebe o destino (padrão: notebook, e nada muda para quem não pedir) e tira
  a largura da tabela. A altura fica escrita na função como proporção da largura.
- Uma função de gravar refaz a mesma figura em cada destino pedido e grava o PDF na pasta de cada
  documento, mais um PNG de conferência, como o `salvar` do guia.
- Quando faltar espaço no slide, rótulos mais curtos ("potência (dBm)" no lugar de "potência total
  de lançamento, soma dos 4 canais (dBm)"). Cada rótulo longo ganha a sua versão curta.

Refazer a figura em vez de reescalar é barato: as funções só desenham, a partir do resultado já
calculado.

## 3. Idioma

- Os textos das figuras saem de um dicionário por idioma, com as mesmas chaves: "em fase (I)" e
  "in-phase (I)", "potência total de lançamento (dBm)" e "launched optical power (dBm)".
- Os números seguem o idioma: em português, vírgula decimal e ponto no milhar (o formatador já
  existe); em inglês, ponto decimal.
- Para o TCC, o dicionário em inglês usa os rótulos fixos do guia (seção 8): "Bit error ratio",
  "Launched optical power per channel (dBm)".

## 4. O que precisa ser decidido

1. **Como pedir.** Um parâmetro em cada chamada (`plot_constellation(r, destino="slides",
   idioma="en")`), uma escolha que vale para as figuras seguintes (`usar_figuras("slides", "en")`),
   ou os dois? Recomendo os dois: a escolha geral como padrão e o parâmetro para uma figura só.
2. **Onde gravar.** Pastas fixas por documento (como no guia) ou passadas por quem chama?
   Recomendo passadas por quem chama, porque o simcopy não sabe onde fica o relatório.
3. **Por onde começar.** Pelas figuras da comparação com o VPI, que são as que vão para relatório e
   slides, ou por todas de uma vez?
4. **Títulos nos notebooks.** As cerca de 210 chamadas de título nas células dos notebooks
   continuam. Tirar os títulos e pôr rótulos de painel e legenda acima, como o guia pede, é
   reescrever o código de gráfico de cada célula. Faço, ou ficam os títulos nos notebooks?

outra coisa, eu gostei dos graficos nos notebooks de imdd e dcs, mas o svg fica grande a ponto de n couber na tela e eu ter que rolar lateralmente. acho q n precisa disso, eu queria ver ele inteiro, e os notebooks anteriores ficaram bem bugados, com legendas acima de figura, mal formatados. se for complexo demais, deixamos esse estilo de formatação pros outros notebooks depois, e mantemos essa formatação só de imdd pra frente.