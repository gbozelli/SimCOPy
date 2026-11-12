# Amplificadores (EDFA): guia de implementação passo a passo

Documento de trabalho para o amplificador óptico, seguindo a aula
(`docs/04-class/context.md`). Marcação igual à de `docs/laser.md` e `docs/moduladores.md`:
**[Fato: fonte]**, **[Dedução]**, **[Proposta]**, **[Dúvida]**, **[Decidido: Gabriel, data]**.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014): seção 7.2.3
   (ruído do amplificador, Eq. 7.2.8–7.2.15), 7.4 (SNR óptica, Eq. 7.4.1–7.4.3), 7.5 (SNR
   elétrica, Eq. 7.5.1–7.5.14).
2. Código do professor: notebook `4_EDFAsAndReceivers` (cópia em `pacote/`), seções 0 e 1;
   `legacy/draft.py`, `EDFA` (linha 1162).
3. Notas da aula: `docs/04-class/context.md`.
4. Projeto: `reference/vpi/parametros_vtmu.md` (módulo `AmpSysOpt`), `docs/arquitetura.md`
   (defeito P-1, correção C-2), `config/system.py` (`AmplifierConfig`).

Testes: com o grupo (**[Decidido: Gabriel, 2026-10-22]** aqui só se implementa).

---

## 1. O que a aula pede

**[Fato: context.md]** "Amplificadores (EDFA, não usaremos Raman): como se agrega ruído no
EDFA, como se calcula SNR na saída e na entrada."

---

## 2. Física (Agrawal)

**[Fato: Eq. 7.2.8–7.2.10]** A figura de ruído é a razão das SNR elétricas na entrada e na
saída, com detector ideal limitado por ruído de disparo: $F_n = \mathrm{SNR}_{in}/\mathrm{SNR}_{out}$,
$\mathrm{SNR}_{in} = P_{in}/(2h\nu\Delta f)$.

**[Fato: Eq. 7.2.11–7.2.12]** Densidade espectral da ASE: $S_{ASE} = n_{sp}h\nu_0(G - 1)$, com
$n_{sp} = \sigma_eN_2/(\sigma_eN_2 - \sigma_aN_1)$ o fator de emissão espontânea.

**[Fato: Eq. 7.2.13–7.2.15]** O ruído dominante é o batimento sinal-ASE;
$\mathrm{SNR}_{out} \approx GP_{in}/[(4S_{ASE} + 2h\nu)\Delta f]$ e
$F_n = 2n_{sp}(1 - 1/G) + 1/G \approx 2n_{sp}$. Mínimo de 3 dB ($n_{sp} = 1$); na prática,
maior que 3 dB, podendo chegar a 6–8 dB.

**[Fato: Eq. 7.4.1–7.4.2]** Numa cadeia de $N_A$ amplificadores com $G = e^{\alpha L_A}$,
$P_{ASE}^{tot} = 2N_AS_{ASE}\Delta\nu_o$, "o fator 2 leva em conta a natureza despolarizada da
ASE"; $\mathrm{SNR}_o = P_{in}/P_{ASE}^{tot}$.

**[Fato: Eq. 7.5.1–7.5.11]** Fotocorrente com ASE copolarizada e ortogonal; variâncias
$\sigma^2_{sig-sp} = 4R_d^2GP_sS_{ASE}\Delta f$, $\sigma^2_{sp-sp} = 4R_d^2S_{ASE}^2\Delta f(\Delta\nu_o - \Delta f/2)$,
disparo $2q[R_d(GP_s + P_{ASE})]\Delta f$ com $P_{ASE} = 2S_{ASE}\Delta\nu_o$; SNR elétrica
(Eq. 7.5.11). A pré-amplificação melhora a SNR quando o ruído térmico domina (Eq. 7.5.14).

**[Dedução]** $S_{ASE}$ é a densidade **por polarização** (o fator 2 das Eq. 7.4.1 e 7.5.10).
Conferido pela conta da Eq. 7.5.8: o batimento ASE-ASE de uma polarização com densidade S numa
banda $\Delta\nu_o$ dá $2R_d^2S^2\Delta f(\Delta\nu_o - \Delta f/2)$; as duas polarizações dão o
dobro, que é a Eq. 7.5.8.

---

## 3. Código do professor

**[Fato: notebook 4, seção 1]** Modelo $s_{out} = \sqrt{G}\,s_{in} + n_{ASE}$, ruído branco com
$\eta_{ASE} = n_{sp}h\nu(G - 1)$ e $n_{sp} \approx N_F/2$; potência
$P_{ASE} = \eta_{ASE}\cdot BW_{sim}$, com $BW_{sim} = f_s$ (banda inteira da simulação).

**[Fato: notebook 4 e `draft.py`]** Problemas no código:
- o ruído é sorteado com `np.random.rand` (uniforme em [0, 1), média não nula); a arquitetura
  registra como defeito P-1 e propõe a correção C-2 (gaussiano circular);
- o `EDFA` vetorial divide $\sqrt{P_{ASE}}$ por 2 em cada quadratura de cada polarização: com
  ruído gaussiano isso daria $P_{ASE}/2$ por polarização, ou seja, $S_{ASE}$ total nas duas, a
  metade do Agrawal;
- no notebook, o `EDFA` vetorial lê a variável global `s` em vez de `s_in`, o `set_OSNR` usa as
  globais `Popt`, `BW_ref` e `N`, e o `SignalGenerator` faz `exp(2πft)` sem o `1j`.

**[Fato: notebook 4, seção 0]** OSNR medida pela FFT: potência na banda do sinal menos o ruído
estimado em duas bandas vizinhas, com banda de referência de 12,5 GHz. Exercícios: OSNR na saída
em função da OSNR na entrada para várias NF (Ex. 3) e para vários ganhos (Ex. 4).

---

## 4. VPI e configuração

**[Fato: `parametros_vtmu.md`, `AmpSysOpt`]** `AmplifierType = PowerControlled`,
`OutputPower = 10^((PowerdBm-30)/10)`, `NoiseFigure = 5` dB, `SaturationPower = 10e-3` W,
`GainMax = 100` dB, `NoiseBins = On`, `NoiseCenterFrequency = TxCarrierFrequency`,
`PolarizFilter = None`. Resumo do arquivo: "um amplificador óptico controlado por potência,
compensando exatamente a perda da fibra".

**[Fato: `config/system.py`]** `AmplifierConfig(mode="target_output_power" | "fixed_gain",
gain_db=None, output_power=None, noise_figure_db=5.0)`; `SpanConfig.amplifier` e
`LinkConfig.booster`.

---

## 5. Implementado (2026-10-22): kernel `kernels/amplifier.py`

Só física, sem configuração:
- `ase_psd(n_sp, gain, frequency)`: Eq. 7.2.11, por polarização.
- `nsp_from_noise_figure(F, G)`: Eq. 7.2.15 invertida, $n_{sp} = (GF - 1)/[2(G - 1)]$
  (**[Dedução]**).
- `amplify(campo, G, S_ASE, f_s, rng)`: $\sqrt{G}\cdot$campo + ruído gaussiano circular de
  densidade $S_{ASE}$ em cada polarização, na banda inteira da grade; sem `rng`, só o ganho.

Notebook `notebooks/04_amplificadores.py`: ASE (espectro, estatística, despolarização), OSNR na
entrada e na saída de um EDFA (com entrada limpa e ruidosa), cadeia de EDFA (Eq. 7.4.1–7.4.2),
SNR elétrica, batimentos e figura de ruído (Eq. 7.5.7–7.5.11, 7.2.15). Tudo bate com o Agrawal
dentro do erro estatístico (~0,5 %).

---

## 6. Proposta de bloco (aguarda as dúvidas)

**[Proposta]** `components/optical/amplifier.py`: `Amplifier(cfg: AmplifierConfig, frequency)`,
com a interface de todo bloco, `__call__(x, ctx)`, entrada e saída ópticas (2, N):
- `fixed_gain`: G de `gain_db`;
- `target_output_power`: G tal que a potência de saída seja `output_power` (como o
  `PowerControlled` do VPI);
- ruído de `ctx.rng_for(<índice do trecho>, "edfa_ase")` (a arquitetura usa essa chave).

---

## 7. Dúvidas

**A1. NF → n_sp.** Usar a Eq. 7.2.15 invertida (a NF medida fica igual à configurada em
qualquer ganho) ou a aproximação $n_{sp} = F_n/2$ do professor (a NF real cai em ganho baixo:
1,8 dB abaixo em G = 3 dB; igual acima de ~15 dB)? Como o VPI define a NF não está no arquivo.

**A2. ASE por polarização.** Seguir o Agrawal ($S_{ASE}$ em cada polarização, total
$2S_{ASE}$) e não o código do professor (total $S_{ASE}$)? A diferença é 3 dB na OSNR.

**A3. Campo de uma polarização (1, N).** O laser do pacote já devolve (2, N). Se chegar (1, N),
somar ASE só nessa polarização (falta a ortogonal, que entra no ASE-ASE e na potência) ou exigir
(2, N)?

**A4. Ganho no modo `target_output_power`.** Calcular G pela potência total de entrada (sinal +
ASE acumulada dos amplificadores anteriores), como um amplificador real controlado por potência?
E a saturação (`SaturationPower` de 10 mW no VPI) entra agora ou fica para depois?

**A5. hν.** Usar a frequência central da simulação (`SystemConfig.center_frequency`) para toda
a banda WDM? O VPI usa `NoiseCenterFrequency = TxCarrierFrequency`.

**A6. Banda do ruído.** ASE na banda inteira da grade, como o professor (o VPI representa o ruído
fora da banda amostrada por "noise bins")?

---

## 8. Decisões de 2026-10-22 (`docs/duvidas_abertas.md`)

- **Potência [Decidido: Gabriel]** Por canal: `AmplifierConfig.output_power` é por canal; o total é
  N vezes.
- **A1 [Decidido: Gabriel]** Agrawal: $n_{sp}$ pela Eq. 7.2.15 invertida (já é o kernel).
- **A2 [Decidido: Gabriel]** Agrawal: $S_{ASE}$ em cada polarização (já é o kernel).
- **A3 [Decidido: Gabriel; minha leitura, a confirmar]** Campo (1, N) vira (2, N) com Y = 0 e a ASE
  entra nas duas polarizações. Aplicado no `amplify`.
- **A4 [Decidido: Gabriel]** "Considere sempre o Agrawal; se houver conflito com o que definimos
  para o sistema, vale o sistema." **[Minha leitura, a confirmar]** Modos e potência como na
  configuração, ganho aplicado ao campo total (sinal + ASE), sem saturação (a configuração não tem).
- A5 e A6 seguem abertas.

## 9. Decisões de 2026-11-12 (`docs/sistema_imdd.md`, seção 7)

- **Saturação [Decidido: Gabriel]** "Pode colocar". Aplicado no bloco `Amplifier`: `saturation_power`,
  g = g0/(1 + P/Ps) (Agrawal, Exercício 7.2), ln G + (G − 1)·P_in/Ps = ln G0 com a potência média de
  entrada, só no modo `fixed_gain`. **[Dúvida]** Saturação no modo controlado por potência (o VPI usa
  `SaturationPowerType = OutputPower`): recusada até decidir.
- Bloco `Amplifier` em `components/optical/amplifier.py` (seção 6 deste guia): feito.
