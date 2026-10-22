# Laser: guia de implementação passo a passo

Documento de trabalho para implementar o laser no simcopy, seguindo a ordem da aula
(`docs/04-class/context.md`). Cada afirmação está marcada:

- **[Fato: fonte]**: está escrito na fonte indicada (equação, seção ou linha de código).
- **[Dedução]**: obtida por mim a partir de equações citadas; a conta está no texto.
- **[Proposta]**: escolha de implementação minha; precisa de aprovação.
- **[Dúvida]**: não está nas fontes; precisa de resposta do Gabriel ou do orientador.
- **[Decidido: Gabriel, data]**: resposta do Gabriel a uma dúvida.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014), cap. 3:
   seção 3.3 (características de lasers) e 3.4 (geração do sinal óptico), exercícios
   3.5 a 3.16. Citado como "Agrawal" com o número da equação, que é estável.
2. Código do professor:
   - `legacy/draft.py` (idêntico a `draft.py` na raiz do repo): laser CW.
   - Notebook `2_Laser_sources` (export do Colab; cópia em
     `pacote/simcopy/docs/2_laser_sources.py`): equações de taxa. É rascunho (seção 3.6).
3. Notas da aula: `docs/04-class/context.md`.
4. Projeto: `docs/arquitetura.md` (v0.9), `docs/decisoes.md`,
   `reference/vpi/parametros_vtmu.md` (no pacote `simcopy_1`, ainda fora do repo).

---

## 1. O que a aula pede

**[Fato: context.md]** "Lasers: equações descritas matematicamente; modelos de rate
equations (mod. direta) e o modelo quase-monocromático (mod. externa); para diferentes
correntes de injeção existem diferentes fluxos de fótons. Gráfico N_c × t para um
I_bias, embaixo um fluxo de fótons φ × t, e embaixo P_out × I, e um de saída pela
entrada (entrada um pulso quadrado I(t), sai uma potência(t))."

Daí saem dois modelos e quatro gráficos:

| Modelo | Usado em | Gráficos |
|---|---|---|
| A. Equações de taxa | modulação direta (IMDD) | N(t) para degrau de corrente; fótons P(t); P_out × I; pulso quadrado I(t) → P_out(t) |
| B. Quase-monocromático (CW com ruído de fase) | modulação externa (MZM, IQ) e oscilador local | (não pedido na aula; ver passo 2) |

**[Decidido: Gabriel, 2026-10-08]** O "modelo quase-monocromático" da aula é o laser CW
com ruído de fase da seção 2.

**[Decidido: Gabriel, 2026-10-08]** O laser devolve **duas polarizações**, forma (2, N),
com todo o campo alinhado em uma (X) e a outra (Y) igual a zero. Vale para os dois modelos.

**[Decidido: Gabriel, 2026-10-08]** O "fluxo de fótons φ" da aula é **S(t)**, a
densidade de fótons do notebook do professor. No código, `phi` fica reservado para a
**fase** do campo (Agrawal, Eq. 3.3.16), como no notebook do professor (`phi_array`).

---

## 2. Modelo B: laser quase-monocromático (CW)

### 2.1 O que as fontes dizem

- **[Fato: Agrawal, 3.4.2]** Na modulação externa o laser opera em corrente constante e
  produz saída CW; o modulador converte a luz CW nos dados.
- **[Fato: Agrawal, Eq. 3.3.33–3.3.34]** O campo é $E(t) = \sqrt{P}\,e^{i\phi(t)}$. Se as
  flutuações de intensidade forem desprezadas, a autocorrelação do campo é
  $\Gamma_{EE}(\tau) = \exp[-\langle\Delta\phi^2(\tau)\rangle/2]$, com
  $\Delta\phi(\tau) = \phi(t+\tau) - \phi(t)$ gaussiano.
- **[Fato: Agrawal, Eq. 3.3.35–3.3.37]** Desprezando o termo das oscilações de relaxação,
  $\Gamma_{EE}$ decai exponencialmente, o espectro é lorentziano e a largura total à meia
  altura é $\Delta\nu = R_{sp}(1+\beta_c^2)/(4\pi P)$.
- **[Fato: draft.py, linhas 924–941]** `single_frequency_laser(t, P, Delta_nu, Freq_offset,
  X_fraction, phi_pol)`: amplitude $\sqrt{P}$; incrementos de fase gaussianos de desvio
  $\sqrt{2\pi\Delta\nu\,t_s}$ acumulados com `cumsum`; fase total
  $2\pi\,\Delta f\,t + \phi_{ruído}$; divide o campo entre X e Y com `X_fraction` e
  defasagem $\pm\phi_{pol}/2$; devolve forma (2, N).
- **[Fato: parametros_vtmu.md, linhas 146 e 467]** VPI: `Linewidth = 500e3` Hz nos
  transmissores e no oscilador local.
- **[Fato: arquitetura.md, seções 4.6.1 e 8]** O projeto prevê `CWLaserConfig` para
  modulação externa, `RateEquationLaserConfig` para direta (4.6.1), e o sorteio da fase
  por `ctx.rng_for(channel.index, "tx_laser_phase")` (8).

### 2.2 Dedução que liga as fontes

**[Dedução: Agrawal 3.3.34 + 3.3.37]** Um espectro lorentziano de largura $\Delta\nu$ tem
autocorrelação $\exp(-\pi\Delta\nu|\tau|)$. Igualando a $\exp[-\langle\Delta\phi^2\rangle/2]$:

$$\langle\Delta\phi^2(\tau)\rangle = 2\pi\,\Delta\nu\,|\tau|$$

A variância cresce linearmente com τ: é um passeio aleatório (processo de Wiener). Na
grade de passo $t_s$, cada incremento é gaussiano com variância $2\pi\Delta\nu\,t_s$, que é
exatamente o que `single_frequency_laser` faz. O código do professor e o livro concordam.

### 2.3 Implementação proposta

**[Proposta]** `src/simcopy/kernels/laser/cw.py`:

```python
def cw_field(n_samples, dt, power, linewidth, freq_offset, rng) -> np.ndarray:
    """Campo (N,) complexo: sqrt(P)·exp[i(2π·Δf·t + φ(t))], φ passeio aleatório com
    incrementos N(0, 2π·Δν·dt). Agrawal Eq. 3.3.33–3.3.37."""
```

Diferenças em relação ao `draft.py`, cada uma por uma regra do projeto:

| `draft.py` | Proposta | Por quê |
|---|---|---|
| `np.random.randn` (gerador global) | `rng` recebido como argumento, vindo de `ctx.rng_for(...)` | Regra 5 do `CLAUDE.md`; C-1 da arquitetura (seção 9) e o diagnóstico de canais WDM com o mesmo ruído de fase |
| a primeira amostra já tem ruído (`cumsum` começa com um incremento) | fase inicial sorteada uniforme em [0, 2π) | **[Decidido: Gabriel, 2026-10-08]** sorteia uniformemente |
| divide X/Y dentro do laser com `X_fraction` e `phi_pol` | saída (2, N): campo em X, Y = 0 | **[Decidido: Gabriel, 2026-10-08]** "todas as entradas vetoriais para PMD"; alinhar com uma polarização. **[Proposta]** o kernel `cw_field` continua devolvendo (N,); o bloco do componente monta o (2, N) |
| nome `freq_noise` para o que é incremento de fase | `dphi` | o nome deve dizer o que a variável é |

### 2.4 Testes (contra forma fechada) e mutações

| Teste | Valor esperado | Mutação que ele pega |
|---|---|---|
| potência constante: $\lvert E\rvert^2 = P$ em toda amostra | P (literal) | amplitude `P` em vez de `sqrt(P)` |
| variância dos incrementos de fase = $2\pi\Delta\nu\,t_s$ | tolerância estatística escrita no teste (qui-quadrado com N amostras) | `np.pi` em vez de `2*np.pi` |
| autocorrelação média $\lvert\Gamma_{EE}(\tau)\rvert/P = e^{-\pi\Delta\nu\tau}$ | Agrawal 3.3.34 | incremento com `dt**2` |
| $\Delta\nu = 0$: tom puro na frequência `freq_offset` (espectro em um bin) | bin de $\Delta f$ | sinal do offset trocado |
| mesma chave de `rng_for` → mesmo campo; chaves diferentes → campos diferentes | igualdade exata | gerador global no lugar do `rng` |

**Estado (2026-10-08):** `cw_field` implementado pelo Gabriel em
`src/simcopy/kernels/laser/cw.py`; verificado que $\lvert E\rvert^2 = P$ em todas as
amostras. Falta a fase inicial uniforme (decidida acima). O teste
`kernels/laser/tests/cw_test.py` ainda não vale:
- o pytest não o coleta, porque a função se chama `cw_field_test` e precisa começar com
  `test`; faltam também `__init__.py` em `kernels/laser/` e `kernels/laser/tests/`;
- a asserção `cw_wave[-1] < 1e-3` contraria a física: o laser CW não atenua no tempo,
  e o esperado é $\lvert E\rvert^2 = P$;
- o offset de 1 GHz com `dt` = 1 ns é igual à taxa de amostragem e aparece em 0 Hz
  (verificado: pico no bin 0); é preciso $f_s > 2\lvert\Delta f\rvert$ e frequências por
  `np.fft.fftfreq`.

---

## 3. Modelo A: equações de taxa (modulação direta)

### 3.1 Equações

**[Fato: Agrawal, Eq. 3.3.1–3.3.3, 3.3.15, 3.3.16]** Laser monomodo, com P = número de
fótons e N = número de elétrons na região ativa (não densidade):

$$\frac{dP}{dt} = GP + R_{sp} - \frac{P}{\tau_p},\qquad
\frac{dN}{dt} = \frac{I}{q} - \frac{N}{\tau_c} - GP$$

$$G = G_N(N - N_0)(1 - \varepsilon_{NL}P),\qquad
\frac{d\phi}{dt} = \frac{1}{2}\beta_c\left[G_N(N - N_0) - \frac{1}{\tau_p}\right]$$

| Símbolo | Significado | Fonte |
|---|---|---|
| $G_N$, $N_0$ | ganho diferencial [1/s] e número de transparência | Eq. 3.3.3; $G_N = \Gamma v_g\sigma_g/V$, $N_0 = N_T V$ |
| $\tau_p$ | tempo de vida de fótons, $\tau_p^{-1} = v_g(\alpha_{mir} + \alpha_{int})$ | Eq. 3.3.4 |
| $\tau_c$ | tempo de vida de portadores | Eq. 3.3.2 e 3.1.3 |
| $R_{sp}$ | emissão espontânea no modo, $R_{sp} = n_{sp}G$, $n_{sp} \approx 2$ | texto após Eq. 3.3.3 |
| $\varepsilon_{NL}$ | compressão de ganho, típico $10^{-7}$, válido para $\varepsilon_{NL}P \ll 1$ | Eq. 3.3.15 |
| $\beta_c$ | fator de aumento de largura de linha, típico 4–8 (InGaAsP) | Eq. 3.3.16 |

**[Fato: Agrawal, Eq. 3.3.8]** Potência emitida por faceta:
$P_e = \tfrac{1}{2}(v_g\alpha_{mir})\,\hbar\omega\,P$.

**[Dedução: Eq. 3.3.4 + 3.3.8]** Como $v_g\alpha_{mir} = (\alpha_{mir}/\alpha_{cav})/\tau_p$:

$$P_e = \frac{\hbar\omega}{2\tau_p}\,\frac{\alpha_{mir}}{\alpha_{mir}+\alpha_{int}}\,P$$

Dessa forma a conversão fótons → watts depende só de $\tau_p$ e da fração de perda pelos
espelhos, sem precisar de $v_g$ e L separadamente. **[Proposta]** usar esse par como
parâmetros de configuração.

**[Dedução: Agrawal, Eq. 3.3.33]** O campo que entra no resto do simulador é
$A(t) = \sqrt{P_e(t)}\,e^{i\phi(t)}$, em √W como o resto do pacote.

### 3.2 Formas fechadas disponíveis (são a base dos testes)

| Grandeza | Expressão | Fonte |
|---|---|---|
| abaixo do limiar ($R_{sp} = 0$) | $P = 0$, $N = \tau_c I/q$ | texto antes da Eq. 3.3.6 |
| população de limiar | $N_{th} = N_0 + (G_N\tau_p)^{-1}$ | idem |
| corrente de limiar | $I_{th} = qN_{th}/\tau_c$ | Eq. 3.3.6 |
| fótons acima do limiar | $P = (\tau_p/q)(I - I_{th})$ | Eq. 3.3.7 |
| eficiência diferencial | $dP_e/dI = (\hbar\omega/2q)\,\eta_d$, $\eta_d = \eta_{int}\alpha_{mir}/(\alpha_{mir}+\alpha_{int})$ | Eq. 3.3.10 |
| oscilações de relaxação | $\Omega_R = [GG_NP_b - (\Gamma_P-\Gamma_N)^2/4]^{1/2}$, $\Gamma_R = (\Gamma_P+\Gamma_N)/2$ | Eq. 3.3.20–3.3.21 |
| resposta de pequeno sinal | $H(\omega_m) = (\Omega_R^2+\Gamma_R^2)/[(\Omega_R+\omega_m-i\Gamma_R)(\Omega_R-\omega_m+i\Gamma_R)]$ | Eq. 3.3.22 |
| banda de 3 dB | $f_{3dB} \approx \sqrt{3}\,\Omega_R/2\pi$ | Eq. 3.3.23–3.3.24 |
| chirp | $\delta\nu(t) = \frac{\beta_c}{4\pi}[G_N(N-N_0) - 1/\tau_p]$ (**[Decidido: Gabriel, 2026-10-15]** δν > 0 é deslocamento para o azul, como na Fig. 3.16) | Eq. 3.4.1 |
| largura de linha | $\Delta\nu = R_{sp}(1+\beta_c^2)/(4\pi P)$ | Eq. 3.3.37 (exige o ruído de Langevin) |
| atraso de ligação, de $I_b < I_{th}$ para $I_{on}$ | $t_d = \tau_c\ln[(I_{on}-I_b)/(I_{on}-I_{th})]$ | **[Dedução]** Eq. 3.3.2 com $P \approx 0$ até N atingir $N_{th}$; não está no livro |

### 3.3 Valores de referência (exercícios do Agrawal)

**[Fato: Agrawal, Exercícios 3.7, 3.8, 3.10]** L = 250 µm, $\alpha_{int}$ = 40 cm⁻¹, índice
modal 3,3, índice de grupo 3,4, λ = 1,55 µm, $G_N = 6\times10^3$ s⁻¹, $N_0 = 10^8$,
$\tau_c$ = 2 ns, $G_P = -4\times10^4$ s⁻¹, $R_{sp} = 2/\tau_p$ (o livro escreve "2/πp";
**[Decidido: Gabriel, 2026-10-08]** é 2/τp).

**[Dedução: contas minhas com Eq. 3.1.5, 3.1.7, 3.3.4–3.3.10, 3.3.20–3.3.24; conferir antes
de virar teste]**

| Grandeza | Valor |
|---|---|
| refletividade $R = [(n-1)/(n+1)]^2$ | 0,2861 |
| $\alpha_{mir}$ | 50,06 cm⁻¹ |
| $\tau_p$ | 1,2593 ps |
| $N_{th}$ | 2,3234 × 10⁸ |
| $I_{th}$ | 18,613 mA |
| P em $I = 2I_{th}$ | 1,4630 × 10⁵ fótons |
| $P_e$ em $I = 2I_{th}$ | 4,138 mW por faceta |
| $dP_e/dI$ ($\eta_{int}$ = 1) | 0,2223 W/A |
| $\Omega_R/2\pi$ (Ex. 3.10, completo) | 4,187 GHz; $\Gamma_R^{-1}$ = 276 ps |
| $f_{3dB}$ (aprox. 3.3.24) | 7,28 GHz |

**[Fato: Agrawal, Exercício 3.14]** Pulso de corrente retangular de 200 ps,
$I_b/I_{th}$ = 0,8, $I_m/I_{th}$ = 3, $\tau_p$ = 3 ps, $\tau_c$ = 2 ns, $R_{sp} = 2/\tau_p$,
$G_N = 10^4$ s⁻¹, $N_0 = 10^8$, $\varepsilon_{NL} = 10^{-7}$. Pede a forma do pulso óptico e o
chirp, e pergunta por que o pulso óptico é mais curto que o de corrente. Esse é o gráfico
"pulso quadrado I(t) → P(t)" da aula.

**[Dedução]** Com esses valores, $I_{th}$ = 10,68 mA e o atraso de ligação é
$t_d \approx$ 138 ps: dos 200 ps do pulso de corrente, ~138 ps são gastos levando N até
$N_{th}$, e só o restante vira luz. Isso responde à pergunta do exercício.

**[Fato: Agrawal, Fig. 3.16]** Referência qualitativa: pulso de 100 ps, $I_b = I_{th}$,
$I_m = 3I_{th}$, $\beta_c = 5$: o pulso óptico sai atrasado, sem bordas abruptas, e o chirp é
primeiro para o azul (frente) e depois para o vermelho (cauda).

### 3.4 Método numérico

**[Decidido: Gabriel, 2026-10-08]** Euler explícito, como no notebook do professor
(passo de 0,1 ps lá).

**[Dedução]** Consequências para os testes:
- No estado estacionário as derivadas são nulas, e o ponto fixo do Euler coincide com o
  da equação contínua. Dessa forma os testes de estado estacionário ($I_{th}$, S(I),
  $dP_e/dI$) podem ter tolerância apertada, independente do passo.
- A dinâmica (oscilações de relaxação, atraso de ligação, resposta de pequeno sinal) tem
  erro proporcional ao passo, porque o Euler é de 1ª ordem. Os testes dinâmicos precisam
  declarar o passo usado e uma tolerância compatível com ele.

O passo da grade do simulador é ~1 ps (ex.: 16 × 59,84 GBd → 1,04 ps), dez vezes o do
notebook (0,1 ps).

**[Decidido: Gabriel, 2026-10-08]** Subpassos de Euler: cada amostra da grade é dividida em
k passos internos, com $k = \lceil dt_{grade}/dt_{máx}\rceil$, e a saída é guardada só nos
instantes da grade. **[Proposta]** $dt_{máx}$ = 0,1 ps (o do notebook) como parâmetro
numérico, ao lado de `ssfm_max_step` em `SimulationConfig`; I(t) interpolada linearmente
entre amostras da grade nos subpassos; com ruído, cada subpasso sorteia o próprio
incremento com variância $2D\,dt_{sub}$.

**Solução alternativa (documentada, não adotada): grade fina.** Exigir que o passo da
grade global já seja ≤ $dt_{máx}$, sem subpassos. Vantagem: o laser integra direto na
grade, sem interpolar I(t) dentro da amostra. Custo: o número de amostras N cresce na
mesma proporção (~10× no exemplo acima), e com ele o custo de todas as FFTs da fibra e dos
filtros, por causa de um único bloco. Reamostrar para uma grade fina só no laser não é
opção na v1: D-4 exclui reamostragem fracionária, e subpassos inteiros já são a forma
inteira dessa ideia.

**[Proposta]** Condição inicial: estado estacionário em $I(0)$ (resolve as equações com
derivadas nulas), para que a simulação não comece com um transiente artificial. O gráfico
"N × t para um I_bias" da aula é justamente o transiente de ligar o laser; nesse caso a
condição inicial é N = P = 0 e o degrau é explícito.

### 3.5 Ruído (Langevin)

**[Fato: Agrawal, Eq. 3.3.25–3.3.28]** Forças de Langevin gaussianas com
$\langle F_i F_j\rangle = 2D_{ij}\delta(t-t')$; dominantes $D_{PP} = R_{sp}P$ e
$D_{\phi\phi} = R_{sp}/(4P)$.

**[Decidido: Gabriel, 2026-10-08]** O ruído entra. Isso revisa D-9 ("sem RIN") para o
laser de modulação direta.

**[Dedução]** Discretização (Euler–Maruyama, coerente com o Euler decidido): de
$\langle F F\rangle = 2D\,\delta(t-t')$, a integral da força num passo dt é gaussiana com
variância $2D\,dt$. Logo o incremento de ruído por passo é $\sqrt{2D\,dt}\cdot\mathcal{N}(0,1)$:
proporcional a $\sqrt{dt}$, não a dt.

**[Dedução]** Na formulação em densidade (S = P/V, $R_{sp}^{dens} = R_{sp}/V$ em
m⁻³s⁻¹, Γ = 1): como $F_S = F_P/V$, $D_{SS} = D_{PP}/V^2 = R_{sp}^{dens}S$; e
$D_{\phi\phi} = R_{sp}/(4P) = R_{sp}^{dens}/(4S)$. Os dois coeficientes coincidem com os do
notebook do professor (linha 648, comentada: `2*R_sp*S`; linha 654: `2*R_sp/(4*S)`); a
diferença está no passo: lá o incremento multiplica dt, aqui $\sqrt{dt}$.

**[Dedução: Agrawal 3.3.37 em densidade]** $\Delta\nu = R_{sp}(1+\alpha_H^2)/(4\pi S)$: o V
cancela. Esse é o teste do ruído de fase.

~~**[Decidido: Gabriel, 2026-10-08]** Manter `R_sp` do notebook,
`R_sp = Γ·n_sp/(τ_c·V)/1e40` com n_sp = 1,5.~~ Revisto em 2026-10-16 (abaixo).

**[Fato: Agrawal, texto após a Eq. 3.3.3]** $R_{sp} = n_{sp}G$, com $n_{sp}$ da ordem de 2
em lasers de semicondutor. **[Fato: Agrawal, Ex. 3.10 e 3.14]** $R_{sp} = 2/\tau_p$.

**[Decidido: Gabriel, 2026-10-16]** $R_{sp}$ do Agrawal. **[Interpretação, avisada]**
$R_{sp} = n_{sp}G$ avaliado no limiar ($G = 1/\tau_p$), constante, como nos exercícios:
$R_{sp} = n_{sp}/(\tau_p V)$ em densidade (2,88e27 m⁻³s⁻¹ no laser do professor). A
alternativa $R_{sp} = n_{sp}G(N, S)$ dinâmica fica negativa abaixo da transparência
(G < 0). Implementado em `resolve_rate_equation_laser`; o teste do resolve mudou de valor
(mudança de especificação).

**[Fato: verificado rodando]** Efeitos da troca no notebook: o joelho da curva P × I
arredonda (erro de 6e-4 contra a forma fechada perto do limiar); ao ligar o laser a luz
aparece 109 ps depois de $t_d$ em 30 mA (antes, 406 ps); no pulso de I_th a 3·I_th o pico sai
em 76 ps e aparece uma segunda oscilação de relaxação dentro do pulso.

**Implementação (2026-10-16, passo 7).** `integrate(..., rng=None)`: com `rng`, Euler–Maruyama
com $\sqrt{2D_{SS}h}\,\xi$ em S e $\sqrt{2D_{\phi\phi}h}\,\xi'$ em φ, normais independentes;
$D_{NN}$ e cruzados desprezados (Agrawal). Sem `rng`, ou com $R_{sp} = 0$, o resultado é o
determinista.

**[Proposta, a confirmar]** Com poucos fótons a força gaussiana deixa de valer (abaixo do
limiar há menos de um fóton na cavidade, e um passo de ruído deixaria S < 0). O integrador
limita S ≥ 0 depois de cada passo e calcula $D_{\phi\phi}$ com S ≥ $R_{sp}\tau_p$
($n_{sp}$ fótons), para não dividir por zero partindo de S = 0.

**[Dedução, verificada rodando]** A largura de linha medida segue a Eq. 3.3.37 entre 1,5 e
6·I_th (notebook, seção 7). Linearizando em frequência zero, o fator de aumento exato é
$1 + \alpha_H^2\kappa$ com $\kappa = [g_cS\,a_{NS}/(a_{SN}a_{NS} + \Gamma_P\Gamma_N)]^2$; com
ε = 0, κ = 1 − 1e-5 para este laser, e o $b$ da Eq. 3.3.35 (que o Agrawal faz igual a 1 na
3.3.37) vale 0,995 em 3·I_th. Com ε > 0 a fase usa o ganho sem compressão (Eq. 3.3.16) e
κ ≈ 1/(1 − εS)²: a largura passa um pouco da Eq. 3.3.37. Não testado com ε > 0.

**[Proposta]** n0 e s0 podem ser arrays (M,): o integrador avança M lasers independentes
juntos, e a saída fica (M, T). Serve para medir estatísticas de conjunto (o teste da largura
de linha usa 4000 lasers em ~1 s).

### 3.6 Notebook do professor (`2_Laser_sources`): formulação em densidade

**[Fato: notebook, linhas 172–185]** Equações implementadas, com N e S em densidade
(m⁻³) e V = L·W·H o volume da região ativa:

$$\frac{dS}{dt} = (G - \tau_p^{-1})\,S + R_{sp},\qquad
\frac{dN}{dt} = \frac{I}{qV} - G\,S - \frac{N}{\tau_c},\qquad
G = \Gamma\,\frac{c}{n_{gr}}\,G_N\,(N - N_{tr})$$

**[Fato: Agrawal, texto após Eq. 3.3.3]** $G_N^{Agrawal} = \Gamma v_g\sigma_g/V$ e
$N_0 = N_T V$: as duas formulações são a mesma equação, dividida ou não por V.

**[Decidido: Gabriel, 2026-10-08]** O simcopy usa a formulação em **densidade**, como o
professor. Onde for preciso converter, funções próprias para isso.

**[Dedução, verificada rodando]** Os parâmetros dos exercícios do Agrawal vêm em número e
sem volume. Para usá-los em densidade basta escolher um V qualquer: $\Gamma v_g\sigma_g =
G_N V$, $N_{tr} = N_0/V$, $R_{sp}^{dens} = R_{sp}/V$. A solução em número ($P = SV$,
$N = N'V$) não depende do V escolhido: com o Ex. 3.14 em I = 2I_th, V = 1e-16, 1,8e-16 e
1 m³ dão o mesmo P = 1,833265e5 fótons.

**[Proposta]** `src/simcopy/kernels/laser/conversions.py`: número ↔ densidade
(`number_to_density`, `density_to_number`), $G_N^{Agrawal}$ ↔ $\Gamma v_g\sigma_g$, e
densidade de fótons → potência por faceta ($P_e = \tfrac{1}{2}v_g\alpha_{mir}\hbar\omega\,SV$,
Eq. 3.3.8). Cada uma com teste de ida e volta.

**[Decidido: Gabriel, 2026-10-08]**
- O notebook é rascunho; `G_N = 1` (linha 112) e `kaixo` (linha 682) são descartados.
  **[Dedução]** Pela equação de limiar do próprio notebook (linha 107), $G_N$ é o ganho
  diferencial $\sigma_g$ = `g_m` = 4e-20 m². Verificado rodando: com `G_N = 1` o laser não
  liga a 20 mA; com `G_N = g_m`, S(20 mA) = 1,46e20 m⁻³ contra 1,49e20 da forma fechada.
- `alpha_mir` = 2410 m⁻¹ e `alpha_int` = 1510 m⁻¹ são os coeficientes de perda (espelhos
  e interna) e ficam como estão.
- Usar **os dois** conjuntos de parâmetros (professor e exercícios do Agrawal): o kernel
  recebe números, e cada conjunto vira um caso de teste com forma fechada.

**Parâmetros do professor** **[Fato: notebook, linhas 40–52]** e valores derivados
**[Dedução: contas minhas, Eq. 3.3.4, 3.3.6, 3.3.7, 3.3.10; conferir antes de virar teste]**:

| Grandeza | Valor |
|---|---|
| λ; L × W × H; V | 1550 nm; 300 × 3 × 0,2 µm; 1,8e-16 m³ |
| $\tau_c$; $N_{tr}$; $\sigma_g$ (`g_m`) | 2,1 ns; 1,25e24 m⁻³; 4e-20 m² |
| Γ; $\alpha_H$; $n_{gr}$ | 1; 3; 3,4 |
| $\alpha_{mir}$; $\alpha_{int}$ | 2410 m⁻¹; 1510 m⁻¹ |
| $\tau_p$ | 2,8932 ps |
| $N_{th}$ | 1,34800e24 m⁻³ |
| $I_{th}$ | 18,512 mA |
| S em 20 mA | 1,4928e20 m⁻³ (2,687e4 fótons) |
| $\eta_d = \alpha_{mir}/(\alpha_{mir}+\alpha_{int})$; $dP_e/dI$ | 0,6148; 0,2459 W/A por faceta |
| no formato do Agrawal: $G_N$; $N_0$ | 1,9594e4 s⁻¹; 2,25e8 (mesma ordem dos exercícios 3.7–3.14) |

**Achados no rascunho, para não copiar** **[Fato: verificado no código; os de execução,
rodando]**:

| Linha | Achado | Correção |
|---|---|---|
| 112 | `G_N = 1` | $G_N = \sigma_g$ (decidido acima) |
| 180, 654 | fase sem `*deltaT`: φ chega a 1e16–1e27 rad | $d\phi = \tfrac{1}{2}\alpha_H(G - \tau_p^{-1})\,dt$ (Agrawal 3.3.16) |
| 648–654 | ruído multiplicado por dt e por um fator 1e25 | $\sqrt{2D\,dt}$ (seção 3.5) |
| 664, 682 | `S[...]` indexa um escalar; `kaixo` | a célula de ruído não roda |
| 371–379 | célula de compressão: atualiza N com `g` da célula anterior, percorre `I_bias_array` em vez de `I_bias_array_GC`, usa `g_m` e não `G_N` | usar o mesmo G nas duas equações |
| 438 | `S_half_m = np.max(...)` | `np.min` (já corrigido na linha 495) |
| 176–177, 184–185 | `max(S, 0)` só no array salvo; S e N internos podem ficar negativos sem aparecer | testar que não ficam negativos |
| — | não há potência em watts, só S | **[Dedução]** $P = SV$; $P_e = \tfrac{1}{2}v_g\alpha_{mir}\hbar\omega P$ (Eq. 3.3.8) |
| 9–20 | a introdução anuncia a seção 2.2 (modelo CW), que não está no arquivo | — |

---

## 4. Onde o código mora

**[Fato: CLAUDE.md, regra 7]** Camadas `core < config < validation < kernels < components <
pipelines`. **[Fato: arquitetura.md, seção 11.1]** Blocos `Laser` e `DirectModulator` em
`components/tx/`.

**[Proposta]** Seguindo a sua ideia de "pasta laser" nos kernels:

```
src/simcopy/kernels/laser/
    __init__.py
    cw.py                 # modelo B: cw_field(...)
    rate_equations.py     # modelo A (densidade): derivadas, Euler com subpassos, estado estacionário
    conversions.py        # número <-> densidade; fótons -> potência (seção 3.6)
    tests/
        test_cw.py
        test_rate_equations.py
src/simcopy/components/tx/laser.py   # depois: blocos que leem a config e chamam os kernels
notebooks/02_laser.py                # gráficos da aula, chamando os kernels
```

Os kernels recebem números (não configurações) e devolvem arrays, como `ssfm.propagate`.
O bloco do componente lê a configuração, pede o `rng` ao contexto e monta o `Waveform`.

**[Proposta]** Configuração (em `config/system.py`), conforme a arquitetura:
`CWLaserConfig(power, linewidth)` e `RateEquationLaserConfig`.

**[Decidido: Gabriel, 2026-10-08]** Grandezas como $\tau_p$ são **derivadas**, nunca
entrada (mesma regra do resto do projeto: derivadas em `config/resolve.py`).

**[Proposta]** Divisão, seguindo os parâmetros do notebook do professor:

| Entrada (`RateEquationLaserConfig`) | Derivada (`resolve`) | Fórmula |
|---|---|---|
| λ, L, W, H | V | L·W·H |
| $n_{gr}$, $\alpha_{mir}$, $\alpha_{int}$ | $v_g$, $\tau_p$ | $c/n_{gr}$; $[v_g(\alpha_{mir}+\alpha_{int})]^{-1}$ (Eq. 3.3.4) |
| $\sigma_g$, $N_{tr}$, Γ, $\tau_c$ | $N_{th}$, $I_{th}$ | $N_{tr} + (\Gamma v_g\sigma_g\tau_p)^{-1}$; $qN_{th}V/\tau_c$ (Eq. 3.3.6) |
| $n_{sp}$ | $R_{sp}$ | $n_{sp}/(\tau_p V)$, Agrawal no limiar (seção 3.5) |
| $\varepsilon$ (compressão), $\alpha_H$ | $\eta_d$, $dP_e/dI$ | $\alpha_{mir}/(\alpha_{mir}+\alpha_{int})$; Eq. 3.3.10 |

---

## 5. Passo a passo

Cada passo termina com `pytest` passando, uma mutação nova detectada por
`scripts/mutacoes.py` e a seção correspondente do notebook `02_laser`.

| Passo | Entrega | Validação |
|---|---|---|
| 1 | Notebook, só conceito e equações (sem código do pacote) | revisão do Gabriel e do orientador |
| 2 | `kernels/laser/cw.py` (modelo B) | testes da seção 2.4 |
| 3 | Estado estacionário das equações de taxa e curva $P_e \times I$ | $I_{th}$, $P(I)$ e $dP_e/dI$ da seção 3.3 |
| 4 | Integrador no tempo; degrau de corrente: N(t), P(t) | estado final = estacionário; $\Omega_R$, $\Gamma_R$ e $t_d$ medidos contra as formas fechadas |
| 5 | Resposta de pequeno sinal | $\lvert H(\omega_m)\rvert$ medido contra a Eq. 3.3.22 |
| 6 | Pulso quadrado e chirp (Ex. 3.14) | forma qualitativa da Fig. 3.16; $\delta\nu(t)$ calculado pela Eq. 3.4.1 a partir do N(t) simulado |
| 7 | Forças de Langevin (aprovado) | $\Delta\nu$ contra a Eq. 3.3.37 em densidade |
| 8 | Componentes `CWLaser` e `DirectModulatedLaser` + configuração | integração com o resto do TX |

Mutações sugeridas: sinal do termo $-P/\tau_p$; $\tau_p$ e $\tau_c$ trocados; fator ½ da
Eq. 3.3.8 removido; sinal de $\beta_c$ no chirp; $\varepsilon_{NL}$ ignorado.

---

## 6. Dúvidas consolidadas

Resolvidas: 1 (fluxo de fótons = S(t)), 2 (fase inicial uniforme), 3 (saída (2, N),
campo em X), 4 (2/τp), 5 (Euler explícito), 6 (ruído entra), 7 (τp e afins derivados),
8 (quase-monocromático = CW), 9 (subpassos; grade fina documentada), 10 (`R_sp` do
notebook), 11 (densidade, com funções de conversão), 12 (Γ é parâmetro de entrada).
Todas em 2026-10-08.

**[Dedução]** Sobre a posição de Γ (a preocupação da dúvida 12): a formulação do
professor é a do Agrawal dividida por V, com S = P/V. Nessa convenção o mesmo
$G = \Gamma v_g\sigma_g(N - N_{tr})$ aparece nas duas equações, e o número de fótons é
P = S·V (sem 1/Γ). A forma com Γ só na equação de S corresponde a outra convenção, com S
medido no volume do modo (V/Γ); o simcopy não a usa.

---

## 7. Estado da implementação

**Atualização (2026-10-22): testes com o grupo.** **[Decidido: Gabriel, 2026-10-22]** Aqui só
se implementa; os testes são escritos pelo grupo (issue #14). Os testes do laser que eu tinha
escrito (`kernels/laser/tests/test_*.py`, `components/tx/tests/test_laser.py`,
`config/tests/test_resolve_laser.py`) e as mutações 10–51 saíram do repo; há uma cópia em
`rascunhos/testes_claude/` (ignorada pelo git). As menções a testes nas atualizações abaixo
são históricas.

**Atualização (2026-10-16): passo 8 (blocos e configuração).**
- `config/system.py`: `LaserConfig` virou `CWLaserConfig` (sem apelido); `chirp_alpha` saiu de
  `DirectModulation` e do `io.py` (**[Decidido: Gabriel, 2026-10-16]**).
- `components/tx/laser.py`, os primeiros componentes do pacote:
  - `CWLaser(laser, frequency, frequency_offset=0, rng_role="tx_laser_phase")`: a entrada só
    dá a grade e o canal; `rng_role` separa o stream do LO (`"rx_lo_phase"`, por exemplo).
  - `DirectModulatedLaser(laser, frequency, bias_ratio, modulation_ratio, noise=True,
    max_step=0.1e-12)`: $I(t) = I_{th}[I_b/I_{th} + (I_m/I_{th})\,x(t)]$ (Eq. 3.3.14,
    **[Decidido: Gabriel, 2026-10-16]** correntes em múltiplos de $I_{th}$), estacionário em
    I(0), saída $\sqrt{P_e}\,e^{i\phi}$ por faceta (**[Decidido]**), stream
    `ctx.rng_for(canal, "tx_laser_noise")`.
  - `rate_equation_params(cfg, resolvido)`: o mapeamento configuração → kernel, que o
    notebook também passou a usar.
  - Os dois devolvem (2, N) com Y = 0, domínio óptico, f0 = frequência do canal, e exigem
    `channel_id` (sem ele não há stream do contexto).
- **[Proposta, a confirmar]** Com `noise=True` a fase inicial do `DirectModulatedLaser` é
  uniforme, como no CW (decisão 2); com `noise=False`, φ(0) = 0.
- Testes em `components/tx/tests/test_laser.py` (10), com entrada elétrica sintética
  (**[Decidido]**); 9 mutações novas. Notebook: seção 9 (OOK sintético a 10 Gb/s).
- **[Pendente, pergunta ao Gabriel]** `DirectModulation.laser` ainda é `CWLaserConfig`. Trocar
  para `RateEquationLaserConfig` (que não tem padrões) quebra `IMDDChannel()` sem `tx`, usado
  em 2 testes de `config/tests/test_system.py`. Opção A: `tx` obrigatório no IMDD e os 2
  testes passam a informar o transmissor. Opção B: o padrão do IMDD vira modulação externa.
- **[Pendente, pergunta ao Gabriel]** Frequência no estacionário. **[Fato: Agrawal, texto
  da Eq. 3.4.1]** o chirp é a mudança da frequência "em relação ao valor de estado
  estacionário ν₀". **[Fato: verificado rodando]** Com $R_{sp}$ e ε, a Eq. 3.4.1 não zera no
  estacionário. Com o ε do Agrawal: −0,64 GHz em $I_{th}$, +0,55 GHz em 1,2·I_th, +2,8 GHz
  em 2·I_th, +5,9 GHz em 3·I_th e +12,7 GHz em 5·I_th (comparar com 75 GHz entre canais). Hoje o bloco devolve a fase crua, ou seja, f0 é a frequência em que
  $g_c(N - N_{tr}) = 1/\tau_p$, e a portadora sai deslocada de f0 por um valor que depende
  de $I_b$. Alternativa: descontar a frequência do estacionário em $I_b$, para a portadora
  ficar em f0.

**Atualização (2026-10-16): passo 7 (ruído de Langevin).** Detalhes na seção 3.5.
- `resolve`: $R_{sp} = n_{sp}/(\tau_p V)$ (Agrawal); `test_resolve_laser` atualizado para o
  valor novo, 2,8803589e27 m⁻³s⁻¹.
- `integrate(..., rng=None)` com Euler–Maruyama; n0 e s0 aceitam arrays (conjunto).
- Testes novos (8): variância por passo do ruído de intensidade e de fase, independência
  das duas forças, semente, R_sp = 0, S ≥ 0 abaixo do limiar, conjunto = realizações
  separadas, largura de linha contra a Eq. 3.3.37 (4000 lasers, tolerância 15 %).
- Mutações novas: 7 (R_sp com τ_c; ruído ∝ h; $D_{\phi\phi}$ sem 1/4; forças correlacionadas;
  S negativo; sem piso em $D_{\phi\phi}$; ruído ligado com R_sp = 0).
- Notebook: seção 7 nova (ruído), CW passa a ser a seção 8; textos das seções 2, 3 e 6
  atualizados para o R_sp novo.
- Os testes do kernel continuam usando o R_sp antigo (3,97e-16) como entrada nos casos
  deterministas: lá ele é só um número de entrada e os testes não dependem do resolve.

**Atualização (2026-10-15): passo 5 (pequeno sinal).** `kernels/laser/small_signal.py`:
`linearize`, `modulation_response` (Eq. 3.3.19), `normalized_modulation_response`
(Eq. 3.3.22), `bandwidth_3db` (Eq. 3.3.23, |H| = 1/2).
- **[Decidido: Gabriel, 2026-10-15]** Linearização completa, sem descartar os termos de
  ordem ε·S; coincide com o Agrawal quando ε = 0.
- **[Decidido: Gabriel, 2026-10-15]** Compressão de ganho do Agrawal: ε_NL = 1e-7 por fóton,
  ε = ε_NL·V (1,8e-23 m³ no laser do professor).
- **[Fato: verificado rodando]** Com o ε do notebook do professor (1e-22 m³), ε·S passa de 0,2
  a partir de ~40 mA; a fórmula literal do livro erra Ω_R em 25 % já em 31 mA e prevê
  ressonância onde a simulação é superamortecida (60 mA: zero cruzamentos).
- Testes em `tests/test_small_signal.py` (5), inclusive a medição por lock-in com o
  integrador (amplitude a 1,5 %, fase a 0,02 rad); 4 mutações novas.
- Notebook: seção 5 nova (pequeno sinal) e seção 6 refeita (pulso de I_th a 3·I_th,
  potência em watts, **[Decidido]** δν > 0 = azul).

**Atualização (2026-10-08, segunda rodada).** CW e equações de taxa agora estão juntos em
`kernels/laser/` no branch `laser-estado-estacionario` (pasta `Área de Trabalho\simcopy\laser`):
- **Passo 2 (CW):** `cw.py` do Gabriel trazido para o branch, com a fase inicial uniforme
  (decidida) e a fase em t = 0 sem incremento. Testes novos em `tests/test_cw.py`
  (potência constante, variância dos incrementos, autocorrelação lorentziana, tom puro,
  fase inicial uniforme, determinismo). O `cw_test.py` antigo ficou só no `repo/`.
- **Passo 4 (integrador):** `integrate(...)`, Euler explícito com subpassos e I(t)
  interpolada. Testes: ponto fixo, convergência após degrau, oscilações de relaxação contra
  a Eq. 3.3.20 (Ω_R a 1 %, Γ_R a 5 %), subpassos equivalentes à grade fina.
- **Notebook `notebooks/02_laser.py`** (+ `.ipynb` executado): limiar, P × I, ligar o laser
  (N(t), S(t)), oscilações de relaxação, pulso quadrado (Agrawal Ex. 3.14) e CW.
- Total: 101 testes; 25/25 mutações (7 novas: 4 do CW, 3 do integrador).

**Passo 3 (estado estacionário e P_e × I): implementado em 2026-10-08**, no branch
`laser-estado-estacionario` (ainda sem commit):
- `config/system.py`: `RateEquationLaserConfig` (só entradas, sem padrões).
- `config/resolve.py`: `resolve_rate_equation_laser(cfg, frequency)` → V, v_g, τ_p,
  g_c = Γv_gσ_g, N_th, I_th, R_sp, h·f, η_d, dP_e/dI.
- `kernels/laser/conversions.py`: número ↔ densidade, G_N → g_c, potência por faceta.
- `kernels/laser/rate_equations.py`: `RateEquationParams`, `gain`, `derivatives`
  (dN/dt, dS/dt, dφ/dt) e `steady_state` (raiz de uma equação em S pelo método de Brent).
- Testes: 18, com os dois conjuntos de parâmetros (professor e Agrawal 3.7, 3.8, 3.10),
  inclusive a independência do volume escolhido na conversão.
- Mutações novas: 8 (de 10 a 17 em `scripts/mutacoes.py`).

**[Proposta, a confirmar]** O ganho não é limitado a G ≥ 0. O notebook do professor
usa `max(G, 0)`; o Agrawal (Eq. 3.3.3) não limita. A diferença só aparece abaixo da
transparência (N < N_tr), onde G < 0 representa absorção.
