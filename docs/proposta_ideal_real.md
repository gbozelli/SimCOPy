# Componentes ideais e reais: proposta para discussão (2026-11-16)

Nasceu da resposta à dúvida 54: a híbrida, o oscilador local e a detecção devem poder ser ideais ou
reais conforme os parâmetros do usuário, e isso vale para os outros componentes, porque para testar
é bom ter a versão ideal de cada um.

**Proposta em uma frase:** cada defeito de cada componente é um parâmetro com um valor neutro que o
desliga, o padrão de `src/` é sempre o valor neutro, e os valores reais ficam nos exemplos em YAML.

Nada aqui foi implementado. A seção 5 lista o que precisa ser decidido.

---

## 1. Como está hoje

A maior parte já funciona assim: há um valor do parâmetro que desliga o defeito. O que não é
uniforme é o padrão: alguns componentes nascem ideais e outros nascem com os valores do VPI.

| Componente | Defeito | Parâmetro | Valor que desliga | Padrão hoje |
|---|---|---|---|---|
| Laser CW | ruído de fase | `laser.linewidth` | 0 | 500 kHz (real) |
| MZM | perda de inserção | `mzm.insertion_loss_db` | 0 | 6 dB (real) |
| MZM | extinção finita | `mzm.extinction_ratio_db` | `.inf` | 30 dB (real) |
| MZM | banda elétrica | `mzm.s21_bandwidth` | `null` | `null` (ideal) |
| Modulador IQ | desbalanço de fase | `iq.phase_imbalance_deg` | 0 | 0 (ideal) |
| Modulador IQ | skew entre I e Q | `iq.skew` | 0 | 0 (ideal) |
| Modulador IQ | desbalanço de ganho | `iq.gain_imbalance_db` | 0 | 0 (ideal) |
| Fibra | atenuação, dispersão, Kerr | `attenuation`, `dispersion`, `n2` | 0 | valores de fibra padrão (real) |
| Fibra | inclinação da dispersão | `dispersion_slope` | `null` | `null` (ideal) |
| Fibra | PMD | `birefringence.enabled` | `false` | `true` (real) |
| Amplificador | saturação | `saturation_power` | `null` | `null` (ideal) |
| Mux | filtro e perda | `mux.optical_filter`, `insertion_loss_db` | `null`, 0 | ideal |
| Fotodetector | ruído térmico | `thermal_noise_density` | 0 | 5,76 pA/√Hz (real) |
| Fotodetector | ruído de disparo | `shot_noise` | `false` | `true` (real) |
| Fotodetector | corrente de escuro | `dark_current` | 0 | 0 (ideal) |
| Oscilador local | laser próprio | `rx.lo_source` | `tx_laser` | `free` (real) |
| Oscilador local | desvio de frequência | `rx.lo_frequency_offset` | 0 | 500 MHz (real) |
| ADC | quantização | `adc.resolution_bits` | 0 | 12 (real) |
| DSP | cego ou com símbolos conhecidos | `rx.dsp.mode` | `data_aided` | `data_aided` |
| Simulação | ruídos do receptor | `simulation.shot_noise`, `thermal_noise` | `false` | `true` |

---

## 2. O que não tem versão real, ou não tem versão ideal

| Componente | O que falta | Parâmetro proposto | Fonte |
|---|---|---|---|
| Híbrida de 90° | erro de fase entre os ramos I e Q | `rx.hybrid.phase_error_deg` | VPI: o transmissor tem `PhaseShift`; o do receptor não está no documento |
| Híbrida de 90° | desbalanço de perda entre as saídas | `rx.hybrid.loss_imbalance_db` | VPI: `InsertionLossImbalanceI/Q` = 0 no receptor (dúvida 18) |
| Detecção balanceada | fotodiodos com responsividades diferentes | `rx.balanced.responsivity_mismatch` | Agrawal, Eq. 10.2.4 e texto sobre ruído de intensidade no capítulo 10: o termo DC e o ruído de intensidade do oscilador local só se cancelam "se os dois ramos forem balanceados" |
| Detecção balanceada | atraso entre os dois fotodiodos | `rx.balanced.skew` | sem fonte ainda |
| Oscilador local | ruído de intensidade (RIN) | `lo.rin_db_hz` | só importa com detecção desbalanceada |
| Oscilador local | divisão entre X e Y | `rx.lo_x_fraction` | professor: `X_fraction = 0.5`, fixo no bloco hoje |
| Laser CW | ruído de intensidade (RIN) | `laser.rin_db_hz` | o laser de equações de taxa já tem, pelo ruído de Langevin |
| Fotodetector | não há como tirar o filtro elétrico | `bandwidth_factor: null` | hoje o filtro sempre entra |
| Amplificador | não há como tirar a ASE | `simulation.ase` (como `shot_noise`) | hoje só sai com ganho 1 |
| ADC | jitter do relógio e banda analógica | `adc.jitter_rms`, `adc.bandwidth_factor` | VPI: não tem no documento |

Os quatro primeiros são os da dúvida 54. Com eles desligados, o receptor é o de hoje.

---

## 3. Três formas de pedir o componente ideal

**A. Só parâmetros, como hoje.** O ideal é escrever os valores neutros.
É o que já existe. O usuário precisa saber quais são os valores neutros, e hoje os padrões misturam
ideal e real: `MZMConfig()` tem 6 dB de perda, e `IQImpairments()` não tem defeito nenhum.

**B. Padrão sempre neutro.** `MZMConfig()` é o MZM ideal, `CWLaserConfig()` é o laser sem ruído de
fase, `FiberConfig(length=...)` continua sendo a fibra padrão (a fibra não tem "ideal" útil). Os
valores reais ficam nos exemplos, como já acontece com o cenário do VPI.
É a regra que a tarefa 1 já tinha: "padrão neutro (impairment desligado)". Muda padrões que os
testes do grupo usam, então precisa do aval de vocês e de issues para os testes.

**C. Atalho por nome.** No YAML, `mzm: ideal`, ou `mzm: {preset: ideal, vpi: 3.0}`. No código,
`MZMConfig.ideal(vpi=3.0)`.
É curto e explícito, mas cria uma sintaxe nova no `io.py`.

**Minha recomendação:** B, mais o construtor `ideal()` de C só no código (para notebooks e testes).
Sem sintaxe nova no YAML. Assim "ideal" e "padrão" passam a ser a mesma coisa, e um teste de bloco
começa de `Config()` e liga um defeito de cada vez, como o notebook 09 faz à mão.

Exemplo de como ficaria o receptor coerente:

```yaml
rx:
  lo: {power: 1.0e-3, linewidth: 0.0}          # padrão: sem ruído de fase
  lo_source: tx_laser                          # padrão: homódino ideal
  hybrid: {phase_error_deg: 0.0, loss_imbalance_db: 0.0}       # padrão: híbrida ideal
  balanced: {responsivity_mismatch: 0.0, skew: 0.0}            # padrão: subtração perfeita
  photodiode: {thermal_noise_density: 0.0, shot_noise: false, bandwidth_factor: null}
```

E o real, no exemplo do VPI, com cada valor escrito.

---

## 4. Laser e modulador (comentário da dúvida 58)

Na montagem eles já são blocos separados: `tx.laser` gera a portadora e `tx.modulator` a modula. O
que está junto é a configuração (`tx:` tem `laser`, `mzm`, `iq`, `rolloff`, `drive_depth`) e a
escolha do modulador, que é fixa: todo canal coerente passa pelo `IQModulator`, e por isso o BPSK
usa um modulador IQ com o braço Q parado.

Se a preocupação é essa, a mudança seria `tx.modulator` virar uma escolha, como o `frontend` do
IMDD já é (`scheme: external` ou `direct`):

```yaml
tx:
  laser: {power: 1.0e-3, linewidth: 500.0e3}
  modulator: {kind: mzm, vpi: 5.0, vbias: 5.0}                       # BPSK com um MZM só
  # modulator: {kind: iq, mzm: {...}, phase_imbalance_deg: 0.0}      # QPSK, QAM
```

---

## 5. O que precisa ser decidido

1. A forma: A, B ou C da seção 3.
2. Se for B: quais padrões mudam primeiro. Cada mudança de padrão pode quebrar um teste do grupo, e
   a regra é não editar teste para passar.
3. Quais dos defeitos da seção 2 entram, e em que ordem. Os da híbrida e da detecção balanceada só
   fazem diferença visível com desbalanço grande, ou com RIN no oscilador local.
4. O que você vê de problema na mistura de laser e modulador: é a configuração (seção 4), ou outra
   coisa?
5. Valores de referência para os defeitos novos. O documento do VPI só tem
   `InsertionLossImbalanceI/Q` = 0 no receptor.
