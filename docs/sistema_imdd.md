# Sistema IMDD: guia

Documento de trabalho da cadeia IMDD montada a partir da configuração (2026-11-12). Mesma convenção
de marcas dos outros guias: **[Fato: fonte]**, **[Dedução]**, **[Proposta]**, **[Dúvida]**,
**[Decidido: Gabriel, data]**; **[Medido]** são números obtidos com o ambiente do repositório.

Fontes, em ordem de prioridade:

1. Agrawal, *Sistemas de Comunicação por Fibra Óptica*, 4ª ed. (trad. 2014): seção 4.6.1
   (Eq. 4.6.6–4.6.11, fator Q e BER).
2. Notebook 07 (`notebooks/07_sistema_imdd.ipynb`): a cadeia montada à mão, que este pacote
   reproduz; notebook 05, seção 4 (a mesma conta de Q, conferida contra a BER contada).
3. Projeto: `docs/mapa_arquitetura.md` (lacunas e molde dos blocos), `docs/proposta_bloco_fibra.md`,
   os guias de cada componente e `docs/duvidas_abertas.md` (dúvidas citadas por número).

Testes: com o grupo (**[Decidido: Gabriel, 2026-10-22]** aqui só se implementa). A seção 9 é o
insumo para as issues.

---

## 1. O que a aula pede

**[Fato: `docs/04-class/context.md`]** Etapas do IMDD:
- sequência pseudoaleatória convertida em sinal elétrico com filtro de Bessel (não Nyquist);
- laser, fibra e fotodetector, monitorando a saída de cada um;
- laser CW, MZM, fibra e fotodetector;
- duas fibras (SMF e DCF), com um amplificador depois da DCF;
- múltiplos spans;
- "MZM";
- múltiplos canais (mux e demux, enlace compartilhado).

O notebook 07 percorre todas, na ordem. Até o MZM, ele monta a cadeia bloco a bloco; daí em diante,
o sistema sai de um arquivo de configuração e roda pelo executor. A primeira seção do notebook
mapeia cada etapa à sua seção.

**[Fato: notebook 07, célula do fotodetector]** Anotação: "20 - 40 km passo 1, estimar BER por fator
Q (IMDD não dá, pois teríamos que simular com 10^10 símbolos)".

**[Fato: `docs/mapa_arquitetura.md`, comentário do Gabriel]** Pedido de um mecanismo simples para
plotar entrada e saída de interesse de cada componente, no tempo e na frequência.

---

## 2. O que roda

Um ou mais canais IMDD OOK, de um arquivo YAML ou JSON até o fator Q de cada canal:

| Estágio (id) | Bloco | Arquivo | Campo do arquivo de cenário |
|---|---|---|---|
| `chC.tx.bits` | `BitSource` | `components/tx/bits.py` | `seed` do canal |
| `chC.tx.driver` | `NRZDriver` | `components/tx/driver.py` | `tx.driver` (filtro) |
| `chC.tx.laser` | `DirectModulatedLaser` | `components/tx/laser.py` | `tx.frontend`, `scheme: direct` |
| `chC.tx.modulator` | `ExternalModulator` | `components/tx/modulator.py` | `tx.frontend`, `scheme: external` |
| `chC.tx.power` | `PowerControl` | `components/optical/power.py` | `launch_power` ou `launch_power_dbm` (opcional) |
| `mux` | `Multiplexer` | `components/optical/wdm.py` | só com mais de um canal; grade de `channel_spacing` |
| `link.booster` | `Amplifier` | `components/optical/amplifier.py` | `link.booster` (opcional) |
| `link.spanK.fiber` | `Fiber` | `components/optical/fiber.py` | `link.spans[K].fiber` |
| `link.spanK.amplifier` | `Amplifier` | `components/optical/amplifier.py` | `link.spans[K].amplifier` (opcional) |
| `chC.rx.demux` | `Demultiplexer` | `components/optical/wdm.py` | `rx.optical_filter` (obrigatório com mais de um canal) |
| `chC.rx.photodiode` | `Photodiode` | `components/rx/photodiode.py` | `rx.photodiode` |

Cada estágio registra as suas entradas (o mux recebe a saída do transmissor de cada canal; o demux
de cada canal recebe a saída do enlace), e o executor segue esse registro. Depois dos receptores,
`kernels/decision.py` mede o fator Q de cada canal com os bits transmitidos (seção 3).

**Fibra compensadora.** Cada elemento de `link.spans` é uma fibra com um amplificador opcional depois
dela. SMF + DCF + amplificador são então dois elementos, `[SMF]` e `[DCF + amplificador]`; vários
trechos são essa dupla repetida. Valores típicos de DCF (Agrawal, seção 8.2): D = −100 ps/(km·nm),
~0,5 dB/km, A_eff ~ 20 µm²; comprimento por D₁L₁ + D₂L₂ = 0 (Eq. 8.2.3).

Todos os blocos seguem o molde dos lasers (`docs/mapa_arquitetura.md`, 4.5): dataclass congelado,
configuração do dispositivo e números de cima (frequência, passo do SSFM) no construtor, processo em
`__call__(x, ctx)`. Os blocos novos têm `id`, que é a chave das sementes (`ctx.rng_for(id, papel)`)
e o nome do estágio. A montagem fica em `pipelines/imdd.py`, e a execução em `pipelines/runner.py`.

Kernels novos: `kernels/filters.py` (os filtros elétricos, que saíram de
`kernels/photodetector.py` e continuam importáveis de lá), `kernels/electrical.py` (NRZ) e
`kernels/decision.py` (Q e BER).

### 2.1 Como rodar

Linha de comando, na pasta `repo/`:

```
simcopy run examples/imdd_ook_mzm.yaml -o resultados/mzm
simcopy run examples/imdd_ook_dml.json
simcopy sweep examples/imdd_ook_dml.yaml --param "link.spans[0].fiber.length" --range 20e3 40e3 1e3 -o resultados/dml
```

`run` imprime o tempo de cada estágio e as métricas; com `-o`, grava `metricas.json`, o olho, a
cadeia e um gráfico por estágio. `sweep` grava `varredura.json` e a curva de Q e BER.

Em Python (é o que o notebook 07 faz no fim):

```python
from simcopy.config.io import load
from simcopy.pipelines.runner import run, sweep
from simcopy.plotting import plot_chain, plot_eye, plot_stage, plot_sweep

system, simulation = load("examples/imdd_ook_mzm.yaml")      # ou .json
r = run(system, simulation)
r.metrics[0].q, r.metrics[0].ber_q                           # canal 0
r.output("link.span0.fiber")                                 # sinal na saída de um estágio
s = sweep(system, simulation, "link.spans[0].fiber.length", np.arange(20e3, 41e3, 1e3))
```

`sweep` aceita qualquer campo do sistema (`channels[0].launch_power`, `link.spans[0].fiber.dispersion`)
ou da simulação (`simulation.root_seed`). Os estágios cujo bloco e cujo contexto não mudam são
reaproveitados do ponto anterior: variando a fibra, o transmissor roda uma vez só.

### 2.2 Exemplos

- `examples/imdd_ook_mzm.yaml`: laser CW de 10 mW, MZM ideal na quadratura (o da última célula do
  notebook 07), 20 km.
- `examples/imdd_ook_dml.yaml` e `.json` (o mesmo cenário): laser de modulação direta com os
  parâmetros do notebook 07, I = I_th·(1,2 + 1,8·x), 20 km.
- Os dois: 10 Gb/s, 2048 bits, 16 amostras por símbolo, Bessel de ordem 4 com −3 dB em 10 GHz no
  acionamento, fotodetector com os padrões da configuração.
- `examples/imdd_wdm_4ch.yaml`: 4 canais (laser CW de 1 mW + MZM, 0 dBm lançados por canal),
  espaçamento de 50 GHz, 80 km de SMF + 12,8 km de DCF + EDFA de 22,4 dB, demux supergaussiano de
  ordem 3 com 30 GHz (seção 8), 32 amostras por símbolo.

---

## 3. A medida: fator Q, instante de decisão e guardas de borda

**[Fato: Agrawal, Eq. 4.6.6–4.6.11]**
$Q = (I_1 - I_0)/(\sigma_1 + \sigma_0)$, $\mathrm{BER} = \tfrac{1}{2}\mathrm{erfc}(Q/\sqrt{2})$,
limiar ótimo $I_D = (\sigma_0I_1 + \sigma_1I_0)/(\sigma_0 + \sigma_1)$.

`measure_ook` (em `kernels/decision.py`) devolve Q, a BER por Q, a BER contada com o limiar ótimo,
os níveis, os desvios, o instante e o trecho medido:

- **Instante de decisão.** É o atraso, de 0 a 2 símbolos, que maximiza Q, com os bits transmitidos
  conhecidos (a busca do notebook 05). Cobre o atraso de grupo dos filtros e a resposta do laser.
- **Guardas de borda.** A FFT trata a janela como periódica, e o início e o fim do sinal não se
  emendam. Por isso a medida descarta, em cada borda,
  $\lceil R_s(\pi|\beta_2|Lf_s + 1\,\mathrm{ns})\rceil + 1$ símbolos (`config/resolve.py`,
  `edge_guard_symbols`; `docs/mapa_arquitetura.md`, 5.7). O primeiro termo é o atraso da componente
  na borda da banda de simulação, e o segundo cobre a memória do transmissor. **[Medido,
  2026-11-12]** Contra o mesmo trecho numa janela 4 vezes maior, sem ruído: o laser direto difere
  em até 6 bits no início, mesmo sem fibra (ele parte do estacionário de I(0), sem a história dos
  bits anteriores); com a fibra, até 6 bits no início e 4 no fim, de 0 a 80 km. A guarda (11 a 20
  bits) cobre todos os casos medidos.
- **Regra de validação** `EdgeGuardFitsWindow`: erro se as guardas não deixam bits para medir; aviso
  com menos de 100 bits medidos (desvio de cada nível com erro relativo acima de ~10 %).
- **Limite da fórmula.** Com interferência entre símbolos, as distribuições não são gaussianas, e
  a BER por Q é uma estimativa. Com o olho fechado (Q perto de 1) há erros suficientes para contar,
  e a BER contada é a referência; `plot_sweep` mostra as duas.

---

## 4. Configuração: o que mudou

**[Fato: `config/system.py`, `config/io.py`]**

| Campo | Antes | Agora | Dúvida |
|---|---|---|---|
| `DirectModulation.laser` | `CWLaserConfig` | `RateEquationLaserConfig`, mais `bias_ratio` e `modulation_ratio`, sem padrão | 5a |
| frontend padrão do IMDD | direta | externa (CW + MZM), que tem padrões | 5b |
| `ExternalModulation.drive_amplitude` | não existia | pico por braço / Vπ, padrão 0,25 (nulo a pico) | 11 |
| `IMDDTransmitter.driver` | não existia | `ElectricalFilterConfig(filter_order, bandwidth_factor, filter_type="bessel")`, sem padrão; ausente = NRZ sem filtro | 50 |
| `MZMConfig.s21_bandwidth` | 40 GHz | `None` = banda infinita, como no VPI | M5 |
| `launch_power` no YAML | lido e descartado | lido (W, ou `launch_power_dbm`) | 3 |
| `amplifier` do trecho no YAML | descartado | lido; `gain_db` do booster também | — |
| `birefringence.enabled` no YAML | não lido | lido | — |
| `shot_noise`, `thermal_noise` da simulação | não lidos | lidos | 6 |
| arquivo JSON | não | pela extensão `.json` | — |
| canal IMDD do YAML | `NameError` | carrega | 43 |
| chave desconhecida | descartada em silêncio | `ConfigError` com o caminho, nas partes novas (tx e rx do IMDD, frontend, filtro do demux, trecho, amplificador, birrefringência) | issue #9 |
| `IMDDReceiver.optical_filter` | não existia | `OpticalFilterConfig(bandwidth, filter_type, filter_order)`, sem padrão; obrigatório com mais de um canal | 30 |
| regra `ChannelSeedsAreDistinct` | não existia | aviso se dois canais têm a mesma semente (mesmos bits, defeito C-1 do rascunho) | — |
| regra `SSFMStepResolvesNonlinearity` | divisão por zero com γ = 0 | ignora a fibra sem Kerr | — |
| `AmplifierConfig.saturation_power` | não existia | W; `None` = sem saturação; só no modo `fixed_gain` (também `saturation_power_dbm` no YAML) | 22 |
| `SystemConfig.mux` | não existia | `MuxConfig(insertion_loss_db=0, optical_filter=None)`: soma ideal por padrão | 29 |
| `Waveform.f0` dos transmissores | frequência absoluta | desvio do canal em relação ao centro; 0 depois do mux | 2 |
| bits do canal | `ctx.rng_for("payload", seed)` | `rng_from_seed(seed, "payload")`: não dependem de `root_seed` | A5 do mapa |

---

## 5. Gráficos de interesse

**[Proposta, implementada]** Resposta ao comentário do Gabriel no mapa
(`docs/mapa_arquitetura.md`, 5.12). Os blocos não plotam (regra 6 do `CLAUDE.md`). O executor
guarda a saída de todos os estágios, e `simcopy.plotting` lê essas saídas:

- `plot_stage(resultado, id)`: entrada e saída de um estágio, no tempo e na frequência. Os painéis
  dependem do domínio do sinal, não do componente: sinal elétrico mostra forma de onda e espectro;
  sinal óptico mostra potência, chirp e espectro óptico; bits aparecem em degraus. Por isso a mesma
  função serve a todos os componentes.
- `plot_io(entrada, saída)`: o mesmo, para dois sinais quaisquer.
- `plot_chain(resultado, channel=c)`: o caminho de um canal no tempo (transmissor dele, mux,
  enlace, receptor dele), um painel por estágio.
- `plot_eye(resultado, c)`: olho no fotodetector do canal, com I₁, I₀, o limiar e os histogramas
  no instante de decisão.
- `plot_sweep(varreduras)`: Q e BER (por Q e contada) ao longo de uma varredura.

As contas ficam em `simcopy.analysis` (potência, chirp com o sinal da convenção do Agrawal,
espectros com média móvel, olho) e devolvem arrays, para os notebooks e uma futura interface
(`docs/proposta_api_interface.md`). A unidade de cada saída vem do atributo `output_unit` do bloco
(`"A"` no fotodetector). O contrato do `import-linter` põe `analysis` logo acima de `core`, e
`plotting` no topo, onde nenhum bloco o importa.

Ainda não entra o que é interno ao componente e não está na saída: N(t) e S(t) do laser, o campo ao
longo de z na fibra. Para isso falta o gravador no contexto (`docs/mapa_arquitetura.md`, 5.2, e
dúvida A10).

---

## 6. O que a montagem ainda recusa

Recusa com `ConfigError` e a causa, em vez de ignorar:

- mais de um canal sem `rx.optical_filter` em cada canal;
- canal coerente (falta o receptor);
- PAM4 (dúvida 51);
- `line_compensation: true` (dúvida 52). A fibra compensadora já entra como mais um elemento de
  `link.spans` (seção 2);
- amplificador controlado por potência que precisaria de ganho menor que 1 (`SignalError`);
- `saturation_power` no amplificador controlado por potência (seção 7, dúvida no fim).

---

## 7. Decisões

**[Decidido: Gabriel, 2026-11-12]**, a partir das propostas desta implementação (respostas também em
`docs/duvidas_abertas.md`):

1. **Dúvida 5: manter.** Modulação direta com o laser de equações de taxa e as razões de corrente
   obrigatórias (5a); padrão do IMDD = modulação externa (5b, opção B), e os testes de
   `test_system.py` passam sem mudança; a união `direct | external` fica (5c).
2. **Dúvida 49: manter.** Bits aleatórios equiprováveis pela semente do canal, não PRBS.
3. **Dúvida 50: "pode revisar e ver".** Revisão das fontes: a aula pede Bessel sem número; o
   Agrawal não especifica o filtro do transmissor; o professor não tem IMDD (só o DAC de Nyquist do
   coerente); o VPI tem no transmissor `LPF_CutoffFrequency = BaudRate` e `LPF_FilterOrder = 4` (o
   tipo do filtro não está no arquivo). Adotado nos exemplos e no notebook 07: Bessel de ordem 4 com
   −3 dB na taxa de símbolos. Os campos continuam obrigatórios no YAML.
4. **Dúvida 11: "não sei".** Mantido `drive_amplitude` = 0,25 como padrão: com V_b = Vπ/2 é a
   excursão de nulo a pico, a de maior extinção (Eq. 10.1.5). O professor só usa sinal pequeno no
   coerente (`ind_mod` = 0,1), e o VPI não tem IMDD.
5. **Dúvidas M9 e 3: pode ser.** Saída física do modulador; o `PowerControl` no fim do transmissor
   fixa a potência quando o canal tem `launch_power`.
6. **Dúvida 6: manter as duas chaves separadas.** O ruído do fotodetector entra se o dispositivo e a
   simulação o pedem.
7. **Dúvida 22 (A4): saturação, "pode colocar".** `AmplifierConfig.saturation_power`: o coeficiente
   de ganho satura como g = g0/(1 + P/Ps) (Agrawal, Exercício 7.2); integrando,
   ln G + (G − 1)·P_in/Ps = ln G0, com G0 o ganho de pequeno sinal (`gain_db`) e P_in a potência
   média de entrada (o EDFA responde à potência média; seção 7.2.2). Vale no modo `fixed_gain`. No
   controlado por potência, a montagem recusa (ver a dúvida abaixo). O resto da A4 fica como estava:
   ganho pela potência total de entrada; hν no centro e ASE na banda inteira (A5 e A6).
8. **Semente (A5 do mapa): a alternativa.** A semente do canal define os bits sozinha
   (`core.context.rng_from_seed`), como o `RandomNumberSeed` do VPI; `root_seed` define os ruídos.
   Trocar `root_seed` mantém os bits.
9. **Dúvida 2: recomendação.** O mux leva o canal para a sua frequência (opção (b)). `Waveform.f0`
   passa a ser o desvio, em relação ao centro da banda, da frequência em que o envelope está
   referenciado: na saída do transmissor, o desvio do canal (campo em banda base do canal); depois do
   mux, 0. Os lasers e o modulador gravam o desvio (`channel_offset`), e o mux confere esse `f0`
   antes de deslocar. hν continua com a frequência absoluta do canal.
10. **Dúvida 29: as duas.** `SystemConfig.mux` (`MuxConfig`): soma ideal por padrão (professor); com
    `optical_filter` e `insertion_loss_db`, filtro por canal e perda, como o WDM_MUX_N_1 do VPI.
11. **Dúvida 30: "não sei".** Estudo na seção 8. Mantido sem padrão; o exemplo usa supergaussiano de
    ordem 3 com 30 GHz (3× a taxa de bits).
12. **Dúvida 31: manter.** O demux deixa o canal na posição dele.
13. **Linha "MZM" da aula: isso.** Repetir os múltiplos trechos com o MZM e comparar com a modulação
    direta no mesmo enlace (notebook 07, seção "Múltiplos spans").

**[Dúvida]** Saturação no amplificador controlado por potência: o VPI usa `PowerControlled` com
`SaturationPower` = 10 mW e `SaturationPowerType = OutputPower`. O que a saturação faz nesse modo
(limitar a potência de saída? o ganho máximo, `GainMax` = 100 dB?) não está no arquivo, e a definição
do VPI (potência de saída) pode não ser o Ps do Agrawal (potência interna). Hoje a montagem recusa
essa combinação em vez de supor.

---

## 8. Verificações feitas

**[Medido, 2026-11-12, com os bits pela semente do canal e o Bessel de ordem 4]**

- **Pacote contra o notebook.** Com os mesmos bits e as mesmas sementes, as saídas de acionamento,
  MZM, fibra e fotodetector do pacote são iguais às da cadeia montada com os kernels: diferença
  máxima 0,0. A perda em 20 km é −4,0000 dB. No notebook 07, o `NRZDriver` dá o mesmo sinal da
  geração manual (diferença 0,0).
- **YAML e JSON** do laser direto dão o mesmo Q (1,511).
- **Back-to-back.** Externo: Q = 988, limitado pelo ruído do receptor. Direto: Q = 5,8, limitado
  pelas oscilações de relaxação.
- **20 km.** Externo: Q = 26,6, sem erros em 2020 bits. Direto: Q = 1,51, com 42 erros (BER contada
  2,1e-2; por Q, 6,5e-2).
- **20 a 40 km.** Externo: Q cai de 26,6 para 9,9. Direto: Q perto de 1 a partir de 24 km, com 153 a
  257 erros contados (0,08 a 0,13), perto da estimativa por Q (0,12 a 0,14).
- **Fibra padrão + compensadora** (60 km de SMF + 9,6 km de DCF, valores do Agrawal, seção 8.2).
  Só SMF: Q = 6,7 com o MZM e 1,2 com a modulação direta (198 erros). Com a DCF: Q = 36,0 e 5,9,
  com a potência recebida em −9,8 dBm. Com um EDFA de 16,8 dB depois da DCF: Q = 43,5 com o MZM,
  e a potência volta a 7,00 dBm.
- **Saturação** (Ps = 10 mW, G0 = 16,8 dB). O ganho medido satisfaz ln G + (G − 1)·P_in/Ps = ln G0
  (resíduo de ~1e-5, da ASE na saída). Com a entrada do enlace (−9,8 dBm) o ganho cai para 15,3 dB,
  e Q quase não muda (43,5). A compressão de 3 dB fica em −5,2 dBm na entrada.
- **Múltiplos trechos** (SMF + DCF + EDFA, 0 dBm lançados). MZM: Q cai de 59,5 (1 trecho) para 8,8
  (20 trechos); até 2 trechos, aproximadamente 1/√N (batimento sinal-ASE). Sem Kerr segue 1/√N até
  20 trechos (13,0 contra 14,5); a queda a mais vem do efeito Kerr acumulado. Modulação direta: Q de
  5,8 a 4,6, limitado pelo próprio laser.
- **4 canais WDM** (`examples/imdd_wdm_4ch.yaml`, demux de 30 GHz). Com 50 GHz entre canais: Q de
  32,8 a 34,9, igual ao canal sozinho no enlace (34,0); sem demux, o canal sozinho tem 30,8 (o demux
  também corta ASE). Com 25 GHz: Q de 22 a 29, pior nos canais do meio; sem Kerr fica quase igual
  (28,2; 22,8; 22,6; 26,5): diafonia linear no demux. Com o mux como o do VPI (filtro de 45 GHz por
  canal e 3 dB de perda): Q ~24, e a potência na saída do mux cai 3 dB.
- **Filtro do demux** (dúvida 30; 4 canais a 50 GHz; Q do pior canal, com os bits e o Bessel de
  ordem 3 da época do estudo): supergaussiano de ordem 3 com 12, 15, 20, 25, 30, 40 e 50 GHz: 18,3;
  17,7; 27,0; 30,9; 32,8; 31,4; 30,9. Retangular de 50 GHz (professor): 30,5. Retangular de 25 GHz:
  20,6 (*ringing*). Gaussiano de 25 GHz: 34,5.
- **Semente.** Trocar `root_seed` mantém os bits e muda os ruídos.
- **`f0`.** No exemplo de 4 canais, o canal 1 sai do modulador com f0 = −25 GHz, e o campo tem
  f0 = 0 do mux em diante.
- **Reaproveitamento.** Numa varredura de 4 canais, o segundo ponto reaproveita os 17 estágios dos
  transmissores e do mux, e todas as saídas são idênticas às de rodar do zero.
- **Amostras por símbolo.** O laser direto com 16 ou 32 amostras por símbolo dá Q com diferença
  menor que 4 % em 0, 10 e 20 km.
- **Suíte e contrato.** Os 73 testes da suíte passam (sem `test_birefringence.py`, dúvida 45). O
  contrato de camadas está mantido, e as 10 mutações de `scripts/mutacoes.py` continuam achando o
  trecho que alteram (o script não foi rodado).

---

## 9. Insumo para as issues de teste (com o grupo)

- `measure_ook`: sinal OOK sem ruído e sem interferência dá Q infinito e zero erro; com ruído
  gaussiano de desvio conhecido, Q bate com (I₁ − I₀)/(2σ); a BER contada segue ½erfc(Q/√2), como no
  notebook 05.
- `arm_amplitudes`: IL = 0 e ER infinita dão ½ e ½; a extinção medida na saída do
  `ExternalModulator` é a ER configurada; a perda medida no nível 1 é a IL.
- `NRZDriver`: −3 dB do filtro em `bandwidth_factor`·R_s; sem filtro, níveis exatos 0 e 1.
- `Fiber`: β₂ e γ iguais aos do `resolve` do sistema; PMD ligada só com campo (2, N); mesma saída
  com a mesma semente.
- `Amplifier`: potência de saída igual ao alvo; ASE com densidade da Eq. 7.2.11; recusa de ganho
  menor que 1; `ConfigError` sem `gain_db` ou sem `output_power`.
- `Photodiode`: corrente média R·P + I_d; banda do filtro.
- `load`: canal IMDD de YAML e de JSON; chave desconhecida recusada com o caminho; campos
  obrigatórios da modulação direta.
- `run`: mesma configuração, mesmo resultado; `sweep` reaproveita o transmissor e dá o mesmo
  resultado que rodar cada ponto do zero.
- `edge_guard_symbols` e a regra `EdgeGuardFitsWindow`.
- `Multiplexer` e `Demultiplexer`: mux + demux retangular de um canal devolve o campo do canal (como
  o notebook 06); diafonia medida igual à da Eq. 6.3.1; potência de cada canal preservada no mux.
- Montagem com N canais: um transmissor e um receptor por canal, entradas do mux e do demux certas;
  `ConfigError` sem `rx.optical_filter`; aviso com sementes repetidas.

---

## 10. Próximos passos

- Gravador no contexto, para N(t) e S(t) do laser e o campo ao longo de z na fibra (A10).
- Receptor coerente (a sequência coerente da aula).
- Fibra compensadora dentro do trecho na configuração (dúvida 52) e PAM4 (dúvida 51).
- Regra de banda ocupada que enxergue o chirp do laser direto (`docs/mapa_arquitetura.md`, 5.10).
- Sensibilidade (BER contra potência recebida), com `launch_power` ou um atenuador.
