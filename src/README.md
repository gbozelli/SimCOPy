# simcopy

Simulador modular de sistemas de comunicações ópticas (WDM, IMDD e coerente digital).

Estado: **fatia vertical das camadas 1 e 2**. O que está aqui roda e é testado; o
resto ainda não existe. A arquitetura completa está em `arquitetura-simulador-optico.md`.

## O que já funciona

| Camada | Módulo | Estado |
|---|---|---|
| Núcleo | `core/units.py` | Constantes e conversões, dono único |
| Núcleo | `core/grid.py` | `TimeGrid` compartilhada |
| Núcleo | `core/signals.py` | `SignalGeometry`, `Bits`, `Symbols`, `Waveform` |
| Núcleo | `core/context.py` | `SimContext.rng_for()`, streams independentes por canal e papel |
| Config | `config/system.py` | Uniões discriminadas: IMDD/DCS, direta/externa, tx e rx sob o `kind` |
| Config | `config/simulation.py` | O *como*, separado do *o quê* |
| Config | `config/resolve.py` | Grandezas derivadas (SpS, N, n_símbolos, β₂, γ, L_eff) |
| Config | `config/io.py` | Carregamento YAML |
| Validação | `validation/rules.py` | 6 invariantes, todos rodam e devolvem todas as violações |
| Kernels | `kernels/ssfm.py` | Split-step **simétrico**, Manakov 8/9, passo adaptativo por fase |
| Kernels | `kernels/birefringence.py` | SU(2) Haar (quatérnios), passo grosso, DGD de 1ª ordem |
| CLI | `cli.py` | `simcopy validate <config.yaml>` |

## Testes

```
python tests/test_physics.py
```

Seis testes contra solução analítica, não contra outro simulador:

- alargamento de pulso gaussiano bate com `T0·√(1+(z/L_D)²)` com erro < 0,1 %
- energia conservada a 1e-10 sem atenuação
- atenuação bate com Beer-Lambert exatamente
- `random_su2` uniforme na esfera de Poincaré (8 bins, desvio < 6 %)
- `random_su2` unitária
- DGD escala com √L

## Uso

```
python -m simcopy.cli validate examples/400zr_4ch_75ghz.yaml
```

## O que falta

Nesta ordem:

1. `components/tx/` — `BitSource`, `SymbolMapper`, `PulseShaper` (RRC), `IQModulator`
   com os impairments do arquivo de referência (**C-25**: desbalanço IQ, skew, S21, ER)
2. `components/optical/` — `Multiplexer`, `Fiber` (envelope fino do kernel),
   `Amplifier`, `Demultiplexer` (gaussiana de ordem 3)
3. `components/rx/` — front-end coerente, `Photodiode` com filtro **Bessel** (C-21)
   e banda derivada de `k·R_s` (C-22), ADC
4. Receptor **data-aided** (equalizador 2×2 por mínimos quadrados + fase data-aided)
   para o spike: mede física sem confundir com convergência de DSP
5. `pipelines/` + spike: reproduzir `A`, `B`, `C` do cenário de referência
6. Só depois: DSP cego (MMA, Stokes, clock recovery, BPS) como marco próprio,
   validado contra o penhasco de baixa potência

## Regra do cenário de referência

Nenhum módulo sob `simcopy/` pode conter a string `400ZR` nem os valores do
cenário. Ele vive em `examples/400zr_4ch_75ghz.yaml` e nos testes de regressão.
