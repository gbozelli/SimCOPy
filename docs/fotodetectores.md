# Fotodetectores: guia de implementação passo a passo

Documento de trabalho para o fotodetector, seguindo a aula (`docs/04-class/context.md`).
Marcação igual à dos outros guias: **[Fato: fonte]**, **[Dedução]**, **[Proposta]**,
**[Dúvida]**, **[Decidido: Gabriel, data]**.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014): seção 4.1
   (Eq. 4.1.1–4.1.9), 4.4 (Eq. 4.4.2–4.4.14), 4.6.1 (Eq. 4.6.3–4.6.11).
2. Código do professor: notebook `4_EDFAsAndReceivers` (cópia em `pacote/`), seção 2;
   `legacy/draft.py`, `photodiode` (linha 962).
3. Projeto: `reference/vpi/parametros_vtmu.md` (`Rx_DualPol`), `docs/arquitetura.md` (P-3,
   P-15, C-21, C-22), `config/system.py` (`PhotodiodeConfig`), `config/simulation.py`.

Testes: com o grupo (**[Decidido: Gabriel, 2026-10-22]** aqui só se implementa).

---

## 1. Física (Agrawal)

**[Fato: Eq. 4.1.1–4.1.3]** $I_p = R_dP_{in}$, $R_d = \eta q/h\nu$.

**[Fato: Eq. 4.1.6–4.1.9]** RC: resposta ao degrau $1 - e^{-t/RC}$, tempo de subida
$T_r = \ln 9\,RC$, banda $\Delta f = 1/[2\pi(\tau_{tr} + \tau_{RC})]$.

**[Fato: Eq. 4.4.2–4.4.5]** Ruído de disparo: densidade **bilateral** $qI_p$ (unilateral
$2qI_p$); $\sigma_s^2 = 2q(I_p + I_d)\Delta f$, com $\Delta f = \int_0^\infty|H_T(f)|^2df$.

**[Fato: Eq. 4.4.7–4.4.10]** Ruído térmico: densidade bilateral $2k_BT/R_L$;
$\sigma_T^2 = (4k_BT/R_L)F_n\Delta f$; as variâncias se somam.

**[Fato: Eq. 4.4.11–4.4.14]** $\mathrm{SNR} = R_d^2P_{in}^2/\sigma^2$; limite térmico
(SNR ∝ $P^2$), NEP.

**[Fato: Eq. 4.6.6–4.6.11]** Limiar ótimo $I_D = (\sigma_0I_1 + \sigma_1I_0)/(\sigma_0 + \sigma_1)$;
$\mathrm{BER} = \tfrac{1}{2}\mathrm{erfc}(Q/\sqrt{2})$, $Q = (I_1 - I_0)/(\sigma_1 + \sigma_0)$.

---

## 2. Código do professor

**[Fato: notebook 4, seção 2.0]** $i_{pd}(t) = h_{pd}(t)*[R(|A_x|^2 + |A_y|^2) + n(t)]$; disparo com
$S_{shot} = 2q\,i_{pd}(t)$ e térmico com $S_{thermal} = 4k_BT/R_L\cdot F_n$ (no código, sem
$F_n$), ambos com variância $S\cdot BW_{sim}$, $BW_{sim} = f_s$; T = 300 K, $R_L$ = 50 Ω.

**[Fato: notebook 4, seção 2.4, e `draft.py`]** `photodiode`: filtro elíptico (passa em BW,
corta em 1,5·BW, 3 dB / 20 dB) aplicado com `lfilter` depois do ruído. Problemas: ignora as
flags `shot_noise` e `thermal_noise` (soma sempre os dois) e lê as globais `BW_sim`, `N` e `fs`.
A arquitetura registra o formato de saída que muda com a flag (P-3) e o elíptico com ordem
dependente de $f_s$ e fase não linear (P-15).

**[Dedução]** Numa corrente real amostrada a $f_s$, o ruído branco ocupa $[-f_s/2, f_s/2]$: a
variância por amostra é a densidade **bilateral** vezes $f_s$ ($qI f_s$; $2k_BTf_s/R_L$). O
professor usa a unilateral vezes $f_s$, o que dobra os dois ruídos (3 dB).

---

## 3. VPI e configuração

**[Fato: `Rx_DualPol`]** `Responsivity = 1.0` A/W, `DarkCurrent = 0`, `ThermalNoise = 10.0e-12`
A/√Hz, `ShotNoise = On`, filtro `Bessel` de ordem 4 com `Bandwidth = 1.15*SymbolRate`.

**[Fato: `config/system.py`]** `PhotodiodeConfig(responsivity=1.0, thermal_noise_density=10e-12,
dark_current=0.0, shot_noise=True, bandwidth_factor=1.15, filter_type="bessel", filter_order=4)`.
**[Fato: `config/simulation.py`]** `SimulationConfig` também tem `shot_noise` e `thermal_noise`.

**[Fato: arquitetura, linha 1443]** $4k_BT/R_L$ a 300 K e 50 Ω dá 18,2 pA/√Hz, "próximo dos
21 pA/√Hz do VPI" (o arquivo do VPI no pacote diz 10 pA/√Hz).

---

## 4. Implementado (2026-10-22): kernel `kernels/photodetector.py`

- `photocurrent(campo, R, I_d)`: $R(|A_x|^2 + |A_y|^2) + I_d$.
- `thermal_density(T, R_L, F_n)`: $\sqrt{4k_BTF_n/R_L}$ [A/√Hz].
- `receiver_noise(I, f_s, rng, thermal, shot)`: disparo com densidade bilateral $qI(t)$ e térmico
  com densidade unilateral `thermal`², ambos vezes a banda da grade (dedução da seção 2).
- `lowpass_response(f, B, "bessel"|"butterworth", ordem)`: filtro analógico com −3 dB em B,
  aplicado na FFT por `apply_response`; `noise_bandwidth` calcula o Δf da Eq. 4.4.4.

Notebook `notebooks/05_fotodetectores.py`: detecção quadrática (batimentos, polarizações
ortogonais, cancelamento por fase), filtros (módulo, atraso de grupo, degrau, Δf), ruído de
disparo e térmico (Eq. 4.4.5, 4.4.9, 4.4.12) e **BER de OOK back-to-back contra
½erfc(Q/√2)**: o primeiro teste de ponta a ponta com forma fechada.

---

## 5. Proposta de bloco (aguarda as dúvidas)

**[Proposta]** `components/rx/photodiode.py`: `Photodiode(cfg: PhotodiodeConfig, symbol_rate)`,
entrada óptica (2, N), saída elétrica (1, N) com a corrente em A; ruído de
`ctx.rng_for(canal, "rx_photodiode")`; filtro com banda `bandwidth_factor·symbol_rate`.

---

## 6. Dúvidas

**F1. Densidade do ruído.** Seguir o Agrawal (bilateral × $f_s$, como o kernel) e não o
professor (unilateral × $f_s$, 3 dB a mais)?

**F2. `thermal_noise_density`.** É a densidade unilateral de corrente ($\sigma_T^2 = i_n^2\Delta f$),
como assumi para o `ThermalNoise` do VPI? E o padrão fica em 10 pA/√Hz (VPI) ou sai de
$4k_BT/R_L$ (professor, 18,2 pA/√Hz)?

**F3. Filtro.** Bessel analógico aplicado na FFT (inclui o atraso de grupo de ~34 ps em 10 GHz),
ou outro? E no IMDD, a banda também é `bandwidth_factor·symbol_rate`?

**F4. Chaves de ruído com dois donos.** `shot_noise` está em `PhotodiodeConfig` e em
`SimulationConfig` (e `thermal_noise` só no segundo). Qual fica?

**F5. Saída.** Corrente em A num `Waveform` elétrico (1, N) com parte imaginária nula (D-5)?
