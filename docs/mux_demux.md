# Mux e demux WDM: guia de implementação passo a passo

Documento de trabalho para o multiplexador e o demultiplexador, seguindo a aula
(`docs/04-class/context.md`). Marcação igual à dos outros guias: **[Fato: fonte]**,
**[Dedução]**, **[Proposta]**, **[Dúvida]**, **[Decidido: Gabriel, data]**.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014): seção 6.2.2
   (multiplexadores e demultiplexadores), 6.3.1 (interferência linear fora da banda,
   Eq. 6.3.1–6.3.4).
2. Código do professor, `legacy/draft.py`: montagem WDM (linhas 2016–2052) e `optical_filter`
   (linha 2058).
3. Projeto: `reference/vpi/parametros_vtmu.md` (`WDM_MUX_N_1`, `FilterOpt`), `docs/arquitetura.md`
   (seção 4.2, P-8, C-11, C-20), `config/system.py` (`SystemConfig.channel_offsets`).

Testes: com o grupo (**[Decidido: Gabriel, 2026-10-22]** aqui só se implementa).

---

## 1. Fatos

**[Fato: Agrawal, seção 6.2.2]** Demultiplexadores usam grades de difração ou filtros e
acopladores; o mesmo dispositivo serve de mux e de demux, pela reciprocidade.

**[Fato: Agrawal, Eq. 6.3.1 e 6.3.4]** Com o filtro no canal m, $P = P_m + \sum_{n\ne m}T_{mn}P_n$;
a penalidade é $\delta_X \approx 10\log_{10}(1 + X)$, $X = \sum_{n\ne m}T_{mn}$.

**[Fato: `draft.py`]** O mux é a soma dos campos dos canais, cada um já deslocado pelo próprio
laser (`Freq_offset_TXk` em $\pm\delta f/2$, $\pm3\delta f/2$, δf = 75 GHz); o demux é
`optical_filter`: máscara retangular na FFT, centrada em `f0`, de largura `BW` (75 GHz no
script).

**[Fato: VPI]** `WDM_MUX_N_1`: perda de 3 dB, filtro `BandPass` `Gaussian` de `FilterOrder = 3`
(`GaussianOrder = 3`) e 70 GHz, 4 canais a 75 GHz. `FilterOpt` (seleção do canal no receptor):
`Gaussian`, ordem 3, 70 GHz, centrado no canal.

**[Fato: arquitetura]** `f0` do `Waveform` é o desvio em relação ao centro da banda de
simulação (seção 4.2). A máscara retangular dá *ringing* tipo sinc e ISI (P-8); correção: filtro
supergaussiano ou Butterworth (C-11). O mux deve declarar a normalização de potência e
registrar a grade (C-20). `SystemConfig.channel_offsets()` dá a grade simétrica.

---

## 2. Implementado (2026-10-22): kernel `kernels/wdm.py`

- `frequency_shift(campo, f, t)`: campo·$e^{+i2\pi ft}$ (eixo da FFT).
- `multiplex(campos, desvios, t)`: soma dos canais deslocados (o mux do professor).
- `bandpass_response(f, centro, B, "rectangular"|"supergaussian", ordem)`: máscara retangular
  do professor, ou $\exp[-\tfrac{\ln 2}{2}(2(f - f_c)/B)^{2n}]$, com $|H|^2 = 1/2$ em $\pm B/2$.
  **[Proposta]** é a minha leitura do "Gaussian de ordem n" do VPI.
- `apply_optical_filter(campo, H)`: filtra cada polarização na FFT.

Notebook `notebooks/06_mux_demux.py`: grade de 4 canais de 50 GBd a 75 GHz (de
`channel_offsets`), filtros e respostas ao impulso, mux + demux (a máscara retangular devolve o
canal no arredondamento; os filtros suaves distorcem), interferência fora da banda (Eq. 6.3.1 e
6.3.4: medido = previsto).

---

## 3. Proposta de blocos (aguarda as dúvidas)

**[Proposta]** `components/optical/wdm.py`: `Multiplexer(grade)` recebe os campos (2, N) dos
canais e devolve o campo agregado; `Demultiplexer(canal, filtro)` devolve o campo do canal.

---

## 4. Dúvidas

**X1. Quem desloca o canal para a sua frequência?** O laser (como o professor e o VPI, com
`EmissionFrequency` por canal) ou o mux (a arquitetura: "soma N canais na grade")? Ligado a isso:
os blocos de laser que escrevi põem em `f0` a frequência absoluta do canal; a arquitetura diz
que `f0` é o desvio em relação ao centro. Corrijo os lasers para o desvio?

**X2. Convenção de frequência** (dúvida M1 de `docs/moduladores.md`): o sinal do desvio no eixo
da FFT.

**X3. Filtros no mux.** O mux do VPI filtra cada entrada (gaussiana de ordem 3, 70 GHz) e tem
3 dB de perda; o do professor é uma soma ideal. Qual?

**X4. Filtro do demux.** Retangular (professor) ou supergaussiano (arquitetura C-11, VPI)? E o
"Gaussian de ordem 3" do VPI é a supergaussiana do kernel, com B a banda de 3 dB?

**X5. Saída do demux.** Devolver o canal na posição dele na grade ou trazê-lo para a banda base?
(No IMDD não importa; no coerente quem faz isso é o oscilador local.)

**X6. Potência.** Potência de lançamento total ou por canal (decisão D-22, ainda aberta), e
quem normaliza (C-20).

---

## 5. Decisões de 2026-10-22 (`docs/duvidas_abertas.md`)

- **X2 [Decidido: Gabriel]** Convenção do Agrawal. Aplicado: o kernel recebe desvios e centros em
  frequência óptica relativa (positivo = azul), leva um canal para +Δ com $e^{-i2\pi\Delta t}$, e os
  filtros são avaliados em `TimeGrid.optical_frequency`. O notebook 06 usa esse eixo.
- **X6 [Decidido: Gabriel]** Potência do amplificador por canal.
- X1, X3, X4 e X5 seguem abertas.

## 6. Decisões de 2026-11-12 (`docs/sistema_imdd.md`, seção 7)

- **X1 [Decidido: Gabriel]** O mux leva o canal para a sua frequência; `Waveform.f0` é o desvio do canal
  na saída do transmissor e 0 depois do mux.
- **X3 [Decidido: Gabriel]** As duas: soma ideal por padrão (`SystemConfig.mux`) ou filtro por canal e
  perda, como o VPI.
- **X4 ["não sei"]** Filtro do demux configurável e sem padrão (`rx.optical_filter`); estudo na seção 8
  do `docs/sistema_imdd.md`. Exemplo: supergaussiano de ordem 3 com 30 GHz para 10 Gb/s a 50 GHz.
- **X5 [Decidido: Gabriel]** O demux deixa o canal na posição dele.
- Blocos `Multiplexer` e `Demultiplexer` em `components/optical/wdm.py`: feitos.
