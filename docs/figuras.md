# Figuras de validação: analítico contra implementação

## Método

Cada figura compara a forma fechada da literatura com a saída do integrador
implementado em `simcopy/kernels/ssfm.py`. O parâmetro variado muda de painel
para painel e o erro relativo entre as duas soluções ocupa painel próprio. O erro
é definido como `max|A_sim − A_ana| / max|A_ana|` sobre o campo complexo, de modo
que discrepâncias de amplitude e de fase são capturadas em conjunto.

A separação por fenômeno é necessária devido à equação não linear de Schrödinger
não possuir solução fechada geral. Isolando um efeito por vez existe solução
exata para comparação; combinando dispersão e não linearidade, a solução exata
existe apenas para o sóliton fundamental. Além desse ponto a verificação passa a
depender de ordem de convergência e de invariantes, o que a figura 05 explicita.

Reprodução: `python scripts/figuras_validacao.py figuras/`. As asserções
numéricas correspondentes estão em `tests/physics/`.

---

## 01 — Atenuação

Referência: `P(z) = P(0)·exp(−αz)`, com α em Np/m obtido de α[dB/m]/4,343.

| Painel | Variação | Conteúdo |
|---|---|---|
| (a) | α = 0,15 / 0,20 / 0,25 dB/km | decaimento analítico e simulado |
| (b) | idem | erro relativo em z |
| (c) | grade e forma do pulso | CW e gaussianas de 10 e 50 ps, com N e fs distintos |
| (d) | passo do SSFM, 50 m a 20 km | erro em função do passo |

O erro permanece na ordem de 1e-15 devido ao operador de atenuação ser exato em
cada passo. Dessa forma o resultado não depende da discretização, o que o painel
(d) confirma, e tampouco da grade temporal ou da forma do pulso, o que o painel
(c) confirma. Verificar essa independência importa porque um acoplamento
espúrio entre perda e parâmetros numéricos indicaria erro na separação dos
operadores.

---

## 02 — Dispersão cromática

Referência: forma fechada da gaussiana com chirp (Agrawal, *Nonlinear Fiber
Optics*, seção 3.2), da qual decorre
`T1/T0 = sqrt[(1 + C·β₂z/T0²)² + (β₂z/T0²)²]`.

| Painel | Variação | Conteúdo |
|---|---|---|
| (a) | chirp C = −2 / 0 / +2 | perfil de intensidade em z = 0,5·L_D |
| (b) | D = 8 / 16 / 22 ps/(nm·km) | FWHM em função de z em quilômetros |
| (c) | quatro pares (T0, D), com L_D de 1,2 a 57 km | FWHM em função de z/L_D |
| (d) | idem (c) | erro relativo no campo |

O painel (a) é sensível ao sinal de β₂: com β₂ < 0 e C > 0 o produto C·β₂ é
negativo, de modo que o pulso comprime antes de alargar, enquanto C < 0 produz
alargamento monotônico. Sem chirp os três casos coincidem, o que implica que o
caso C = 0 verifica apenas |β₂| e não distingue o sinal.

O painel (c) mostra que quatro combinações de T0 e D com comprimentos de
dispersão que diferem por fator 47 colapsam sobre uma única curva quando a
distância é normalizada por L_D = T0²/|β₂|. É possível aferir que o integrador
respeita a lei de escala da dispersão devido ao colapso ocorrer sem ajuste, e
isso é evidência mais forte que a concordância em um único conjunto de
parâmetros. O erro permanece em 1e-13, compatível com acúmulo de arredondamento
nas transformadas.

---

## 03 — Automodulação de fase

Referência (Agrawal, seção 4.1), válida para β₂ = 0 e propagação escalar:
`A(z,T) = A(0,T)·exp(−αz/2)·exp[i·γ|A(0,T)|²·L_eff]`, com
`L_eff = (1 − e^{−αz})/α`.

| Painel | Variação | Conteúdo |
|---|---|---|
| (a) | φ_max = 0,5π / 1,5π / 3,5π | perfil da fase não linear |
| (b) | idem | espectro em dB |
| (c) | idem | erro relativo no campo |
| (d) | α = 0 / 0,20 / 0,35 dB/km | φ_max em função de L |
| (e) | A_eff = 60 / 80 / 110 µm² | φ_max em função de P0 |
| (f) | α | erro do painel (d) |

A fase não linear reproduz o perfil de intensidade de entrada, com máximo no
mesmo instante, devido a φ_NL ser proporcional a |A(0,T)|². O número de picos do
espectro segue φ_max/π + 1/2, o que decorre da interferência entre instantes de
mesma frequência instantânea.

O painel (d) separa L de L_eff: com α = 0 a fase cresce linearmente com L,
enquanto com α > 0 a fase satura em γP₀/α, indicada pelas linhas pontilhadas.
Dessa forma é possível verificar que o integrador usa o comprimento efetivo e
não o geométrico. O painel (e) varia γ através de A_eff, com
γ = 2πn₂/(λ·A_eff): as três retas têm inclinação γ·L_eff, o que informa que a
dependência com a área efetiva está correta.

Observação sobre o fator 8/9: estes painéis usam propagação escalar
(`manakov=False`). Com duas polarizações a equação de Manakov introduz o fator
8/9 no termo não linear, e a fase máxima passa a ser (8/9)·γP₀L_eff. O teste
`test_03_nao_linearidade.py::test_d` cobre esse caso.

---

## 04 — Dispersão dos modos de polarização

Referência (Agrawal, *Fiber-Optic Communication Systems*, capítulo 2):
`σ_T² = 2(Δβ₁·l_c)²·[e^{−z/l_c} + z/l_c − 1]`, com `D_p = Δβ₁·sqrt(2·l_c)`.

| Painel | Variação | Conteúdo |
|---|---|---|
| (a) | D_p = 0,05 / 0,1 / 0,2 ps/√km, com l_c = 50 m | DGD rms em função de z |
| (b) | l_c = 25 / 50 / 100 m, com D_p = 0,1 ps/√km | DGD rms em função de z |
| (c) | ambos | razão entre simulado e analítico |
| (d) | — | distribuição do DGD a 20 km |

O painel (b) mostra que, para z muito maior que l_c, as três curvas convergem
para D_p·√z independentemente de l_c, de modo que o comprimento de correlação
determina apenas onde ocorre a transição entre os regimes linear e de raiz. Isso
é relevante para a configuração de referência devido ao 400ZR operar entre 80 e
150 km contra l_c de 50 m, ou seja, três ordens de grandeza acima da transição.

O painel (c) quantifica a diferença entre os dois modelos: o implementado usa
seções discretas de comprimento l_c, enquanto o do Agrawal supõe birrefringência
com correlação exponencial. A razão vale aproximadamente √2 em z ≈ l_c e tende a
1 a partir de cerca de 5·l_c. A discrepância é portanto de modelo e não de
implementação, e restringe-se ao regime que não é usado.

O painel (d) verifica a forma da distribuição pela razão ⟨τ⟩/τ_rms, cujo valor
teórico para a Maxwelliana é sqrt(8/3π) = 0,921.

---

## 05 — Empilhando fenômenos

| Painel | Conteúdo |
|---|---|
| (a) | erro contra a forma fechada em quatro etapas cumulativas |
| (b) | sóliton fundamental após 1, 2 e 4 períodos |
| (c) | erro em função do passo, contra referência de passo fino |
| (d) | conservação de energia com β₂, γ e PMD simultâneos |

As etapas do painel (a) e o erro típico de cada uma:

| Etapa | Referência analítica | Erro |
|---|---|---|
| α | exponencial | 1e-14 |
| α + β₂ | gaussiana com chirp | 1e-13 |
| α + γ | SPM pura | 1e-12 |
| β₂ + γ | sóliton fundamental, N = 1 | 1e-4 |
| β₂ + γ + PMD | não existe | — |

O salto de 1e-12 para 1e-4 ocorre devido à separação dos operadores: enquanto um
único efeito atua, o passo do split-step é exato para aquele operador, e o erro
se restringe ao arredondamento. Quando dispersão e não linearidade atuam
simultaneamente, os operadores não comutam, dessa forma surge o erro de
separação, que é o termo controlado pelo passo. É por esse motivo que o painel
(c) verifica a ordem de convergência: o erro cai com h², o que confirma a
implementação simétrica do split-step. Uma implementação assimétrica cairia com
h, e essa distinção não seria visível em nenhum dos painéis anteriores.

Na última etapa não existe solução fechada, de modo que a verificação passa a
usar dois critérios indiretos. O primeiro é a convergência do próprio método,
no painel (c). O segundo é a conservação de energia com α = 0, no painel (d),
que vale porque tanto o termo não linear quanto a rotação de polarização são
unitários, e portanto não podem alterar a energia total do campo.
