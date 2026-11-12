# Proposta: bloco `Fiber` (fibra como componente)

Rascunho de 2026-11-05, para revisão do grupo antes de implementar. Mesma convenção de marcas do
`docs/laser.md`: **[Fato: fonte]**, **[Dedução]**, **[Proposta]**, **[Dúvida]**.

## 1. Situação hoje

**[Fato: código em `src/simcopy`]**

- O kernel existe: `kernels/ssfm.py`, `propagate(A, omega, p, rng=None, max_step=5e3, ...)`, com
  `A` de forma `(n_pol, N)`. Devolve só o campo no fim da fibra.
- A configuração existe: `FiberConfig` (comprimento, atenuação, dispersão, n₂, A_eff,
  `birefringence`) e `resolve_fiber`, que calcula β₂, γ, α em unidades SI.
- **O bloco não existe**: `components/optical/__init__.py` está vazio. O laser já tem o seu
  (`components/tx/laser.py`), com a interface de todo bloco, `bloco(x, ctx) -> Waveform`.

Consequência: quem quer ligar o laser à fibra (notebook 07, E3 do `docs/planejamento.md`) precisa
montar `FiberKernelParams` à mão a cada uso (seis linhas) e embalar o resultado com
`with_samples`. O mapeamento config → kernel fica duplicado em cada notebook, e a PMD (que está na
`FiberConfig` mas **não** passa por `resolve_fiber`) é esquecida sem aviso.

## 2. Proposta

**[Proposta]** Um bloco no mesmo molde do laser, em `components/optical/fiber.py`:

```python
def fiber_kernel_params(cfg: FiberConfig, r: ResolvedFiber) -> FiberKernelParams:
    """Numeros que o kernel recebe: entradas da configuracao + grandezas do resolve."""
    b = cfg.birefringence
    return FiberKernelParams(
        length=cfg.length, beta2=r.beta2, gamma=r.gamma, alpha_np=r.alpha_np,
        pmd_coefficient=b.pmd_coefficient if b.enabled else 0.0,
        correlation_length=b.correlation_length)

@dataclass(frozen=True)
class Fiber:
    fiber: FiberConfig
    index: int = 0                # posicao no enlace, so para a semente da PMD
    max_step: float = 5e3         # passo maximo do SSFM [m] (SimulationConfig.ssfm_max_step)

    def __call__(self, x: Waveform, ctx: SimContext) -> Waveform:
        x.require(Domain.OPTICAL, "Fiber")
        if x.f0 <= 0:
            raise SignalError("Fiber precisa de x.f0 (frequencia do canal) para calcular beta2.")
        p = fiber_kernel_params(self.fiber, resolve_fiber(self.fiber, x.f0))
        rng = (ctx.rng_for("fiber", self.index, "pmd")
               if x.n_pol == 2 and p.pmd_coefficient > 0 else None)
        return x.with_samples(propagate(x.samples, x.grid.omega, p, rng, self.max_step))
```

Uso no notebook 07 (substitui as seis linhas atuais):

```python
from simcopy.components.optical.fiber import Fiber

fiber_waveform = Fiber(FiberConfig(length=20e3))(optical_waveform, ctx)
```

O resto do bloco segue o laser: classe imutável, `require(Domain.OPTICAL)`, ruído/sorteio vindo de
`ctx.rng_for(...)`, saída com a mesma grade, geometria, canal e f0 da entrada.

## 3. Pontos a decidir

| # | Ponto | Por que importa | Recomendação |
|---|---|---|---|
| 1 | De onde vem a frequência central para β₂ e γ | **[Fato: `signals.py`]** `Waveform.f0` existe e o laser o preenche. Usar `x.f0` evita um segundo `frequency` no bloco que poderia divergir do sinal. Mas `f0` vale 0 por padrão. | Usar `x.f0` e dar erro se for 0 (como no esboço). |
| 2 | O kernel sorteia rotações mesmo com PMD zero | **[Fato: `ssfm.py:74`]** a condição é `rng is not None and pmd_coefficient >= 0`: com `rng` dado, a polarização é embaralhada ainda que `pmd_coefficient = 0`. | O bloco só passa `rng` com PMD > 0 (como no esboço). Alternativa: trocar `>=` por `>` no kernel. |
| 3 | PMD ligada por padrão | **[Fato: `system.py`]** `BirefringenceConfig.enabled = True`, 0,1 ps/√km. O laser direto emite só em X, e a PMD joga parte da potência em Y. | Manter o padrão da config; o notebook 07 pode usar `enabled=False` para isolar dispersão. |
| 4 | Passo máximo | `ctx` não carrega `SimulationConfig`; o `ssfm_max_step` não chega ao bloco. **[Dúvida 34]** com PMD, fibras ≤ 5 km caem em um passo só e o DGD sai sempre igual. | Campo `max_step` no bloco, 5 km por padrão; limitar a 1/10 do comprimento quando a PMD estiver ligada. |
| 5 | Atraso de grupo | **[Fato: `01_fibra.py`, seção 0]** o kernel usa o tempo T de um referencial que anda com o pulso; a saída não é atrasada por L/v_g. **[Dedução]** não precisa de sincronismo extra no receptor, mas não dá para ler o atraso real. | Documentar no docstring. |
| 6 | Janela periódica | **[Dedução]** a FFT trata a janela como periódica: o alargamento de dispersão que passa da borda reaparece no outro lado. Ordem de grandeza (D·L·Δλ, com Δλ ≈ 0,4 nm do espectro do notebook 07): ~130 ps em 20 km contra 10 ns de janela; cresce linearmente com L e a janela do caso real pode ser bem menor. | Regra de validação (aviso) quando o alargamento estimado passar de uma fração da janela. |
| 7 | Evolução ao longo de z | **[Dúvida 32]** o kernel devolve só o fim. Gráficos tipo "potência × (t, z)" exigem repropagar a cada z. | Fora desta proposta; se vier, o ponto de entrada é um parâmetro `positions` no kernel. |
| 8 | SMF + DCF + EDFA | **[Dúvida 52]** `SpanConfig` tem uma fibra só. | Fora do escopo: o bloco é uma fibra; o trecho SMF → DCF → EDFA é composição de blocos. |

## 4. O que não muda

- O kernel `propagate` (e a física que o notebook 01 valida) fica como está. Único ajuste opcional:
  o ponto 2.
- `FiberConfig` e `resolve_fiber` ficam como estão.
- `Probe` (E4): o bloco não conhece o monitor, como o laser.

## 5. O que verificar (insumo para a issue de testes, que fica com o grupo)

Os critérios abaixo já estão no `docs/planejamento.md` (semana 1) ou saem do notebook 01:

- Potência média na saída = $P\,e^{-\alpha L}$ (com 20 km e 0,2 dB/km, −4,00 dB). **[Fato: teste manual, 2026-11-05]**
  o laser direto do notebook 07 deu −4,00 dB na média e −4,00 dB na energia do espectro.
- Com β₂ = γ = 0 e PMD desligada, a saída é a entrada multiplicada por $e^{-\alpha L/2}$.
- A dispersão sozinha não muda o módulo do espectro óptico (só a fase); o Kerr muda.
- O bloco devolve o mesmo `grid`, `geometry`, `channel_id` e `f0` da entrada.
- Erros: entrada elétrica, `f0 = 0`.
- Mesma semente → mesma saída (PMD ligada, duas polarizações).

## 6. Arquivos

| Arquivo | Mudança |
|---|---|
| `src/simcopy/components/optical/fiber.py` | novo: `fiber_kernel_params` e `Fiber` |
| `src/simcopy/components/optical/__init__.py` | exportar `Fiber` |
| `notebooks/07_sistema_imdd.ipynb` | usar o bloco no lugar das seis linhas manuais |
| issue "Testes de Fibra" | formato das #5–#10, #14, #15; itens da seção 5 |

---

## 7. Nota (2026-11-09): o MZM está na mesma situação

**[Fato: código em `src/simcopy`]** O modulador também não é um componente:

- Só existe o kernel `kernels/modulator/mzm.py` (`mzm_field`, `push_pull`, `iq_field`,
  `arcsin_predistortion`). `components/tx/` só tem o `laser.py`.
- O kernel trabalha com **arrays**, não com `Waveform`: recebe o campo e as tensões e devolve um
  array. Quem usa precisa embalar o resultado (`carrier.with_samples(...)`) e levar grade, domínio,
  canal e f0 junto.
- A `MZMConfig` tem `vpi`, `vbias`, `insertion_loss_db`, `extinction_ratio_db` e `s21_bandwidth`,
  mas o kernel só usa `vpi` e `vbias` (o docstring diz: "nada de configuração, perda de inserção
  ou filtro S21"). Perda, extinção finita e S21 estão decididos em `docs/moduladores.md` (§8, M4,
  M5 e M9) e ficam para o bloco.
- **[Dúvida M3, aberta]** Como o sinal elétrico (adimensional) vira tensão. Hoje cada notebook
  escolhe a sua escala.
- `iq_field` é o modulador em quadratura (coerente), com assinatura
  `iq_field(field, v_i, v_q, vpi, v_bias, quadrature_phase, ...)`. Para IMDD com modulação externa
  é `mzm_field` com `push_pull`.

A proposta do bloco já está em `docs/moduladores.md` §6 (`ExternalModulator(laser, mzm, frequency)`
em `components/tx/modulator.py`, entrada elétrica e saída óptica `(2, N)`, com o `CWLaser` dentro).
Como o `Fiber`, ele deveria seguir o molde do laser: `__call__(x, ctx) -> Waveform`, `require(...)`
na entrada e saída com a mesma grade, geometria, canal e f0.

**Uso manual enquanto não existe o bloco** (laser CW + MZM em quadratura, OOK de nulo a pico):

```python
carrier = CWLaser(CWLaserConfig(power=10e-3, linewidth=500e3), f0)(electrical_waveform, ctx)
drive_voltage = 0.5 * mzm_config.vpi * (electrical_waveform.samples[0].real - 0.5)   # x em [0, 1] -> ±Vπ/4
v1, v2 = push_pull(drive_voltage, mzm_config.vbias)                                   # vbias = Vπ/2
mzm_waveform = carrier.with_samples(mzm_field(carrier.samples, v1, v2, mzm_config.vpi))
```

Verificado em 2026-11-09: `P_saída/P_portadora = sin²(πx/2)` (erro de 5·10⁻¹⁶), polarização Y
zerada, e a fibra seguinte dá −4,00 dB em 20 km.
