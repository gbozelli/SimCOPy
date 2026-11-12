"""Execucao: carregar, resolver, validar, montar, rodar e medir (arquitetura, secao 6.3).

    from simcopy.pipelines.runner import run_file, sweep
    r = run_file("examples/imdd_ook_mzm.yaml")
    r.metrics[0].q                      # fator Q do canal 0
    r.output("link.span0.fiber")        # sinal na saida de qualquer estagio
    s = sweep(r.system, r.simulation, "link.spans[0].fiber.length", np.arange(20e3, 41e3, 1e3))

Cada estagio guarda a saida (as "sondas": docs/mapa_arquitetura.md, 5.12), para os graficos de
simcopy.plotting. Numa varredura, os estagios cujo bloco e cujo contexto nao mudaram sao
reaproveitados do ponto anterior: variando o comprimento da fibra, o transmissor roda uma vez so.
"""
from __future__ import annotations
import dataclasses
import re
import time
from dataclasses import dataclass, field
from typing import Any, Callable
import numpy as np
from ..config.io import load
from ..config.resolve import Resolved, resolve
from ..config.simulation import SimulationConfig
from ..config.system import SystemConfig
from ..core.context import SimContext
from ..core.units import w_to_dbm
from ..kernels.decision import QMeasurement, measure_ook
from ..validation.base import Severity, ValidationFailed, Violation
from ..validation.rules import validate
from .imdd import build

@dataclass
class Stage:
    id: str
    block: Any
    output: Any                       # Bits ou Waveform
    seconds: float                    # tempo de execucao (0 se reaproveitado)
    reused: bool = False
    inputs: tuple[str, ...] = ()      # ids dos estagios que alimentam este (varios no mux)

@dataclass
class Result:
    system: SystemConfig
    simulation: SimulationConfig
    resolved: Resolved
    context: SimContext
    stages: list[Stage]
    metrics: dict[int, QMeasurement]  # por indice de canal
    violations: list[Violation] = field(default_factory=list)   # avisos da validacao

    @property
    def ids(self) -> list[str]:
        return [s.id for s in self.stages]

    def stage(self, stage_id: str) -> Stage:
        for s in self.stages:
            if s.id == stage_id:
                return s
        raise KeyError(f"estagio {stage_id!r} nao existe; estagios: {self.ids}")

    def output(self, stage_id: str):
        """Sinal na saida do estagio."""
        return self.stage(stage_id).output

    def input(self, stage_id: str):
        """Sinal na entrada do estagio (saida do estagio que o alimenta; None na fonte de bits).
        No mux, a entrada do primeiro canal."""
        st = self.stage(stage_id)
        return self.stage(st.inputs[0]).output if st.inputs else None

    def path(self, channel: int = 0) -> list[Stage]:
        """Estagios no caminho de um canal: o transmissor dele, o mux, o enlace e o receptor dele."""
        return [s for s in self.stages
                if s.id.startswith(f"ch{channel}.") or not s.id.startswith("ch")]

    def received_power(self, channel: int = 0) -> float:
        """Potencia optica media na entrada do fotodetector [W]."""
        return self.input(f"ch{channel}.rx.photodiode").power

    def summary(self) -> str:
        lines = []
        for ch, m in self.metrics.items():
            q_db = 20 * np.log10(m.q) if np.isfinite(m.q) and m.q > 0 else float("inf")
            lines.append(
                f"canal {ch}: Q = {m.q:.3f} ({q_db:.2f} dB) | BER por Q = {m.ber_q:.3e} | "
                f"erros contados = {m.errors}/{m.n_bits} | potencia recebida = "
                f"{w_to_dbm(self.received_power(ch)):.2f} dBm | instante = {m.delay} amostras | "
                f"simbolos {m.first_symbol} a {m.last_symbol - 1}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Metricas, estagios e tempos, para gravar em JSON."""
        return {
            "metrics": {str(ch): {**dataclasses.asdict(m),
                                  "received_power_dbm": w_to_dbm(self.received_power(ch))}
                        for ch, m in self.metrics.items()},
            "stages": [{"id": s.id, "block": type(s.block).__name__, "inputs": list(s.inputs),
                        "seconds": s.seconds, "reused": s.reused} for s in self.stages],
            "warnings": [str(v) for v in self.violations],
            "grid": {"sampling_rate": self.resolved.grid.sampling_rate,
                     "n_samples": self.resolved.n_samples, "window": self.resolved.window},
        }

def _measure(system: SystemConfig, resolved: Resolved, stages: list[Stage]) -> dict[int, QMeasurement]:
    by_id = {s.id: s.output for s in stages}
    out = {}
    for c, rc in zip(system.channels, resolved.channels):
        bits = by_id[f"ch{c.index}.tx.bits"]
        current = by_id[f"ch{c.index}.rx.photodiode"]
        out[c.index] = measure_ook(current.samples[0].real, bits.data[0], rc.samples_per_symbol,
                                   guard=rc.edge_guard_symbols)
    return out

def run(system: SystemConfig, simulation: SimulationConfig,
        previous: Result | None = None) -> Result:
    """Roda o sistema. Erros de validacao interrompem antes de alocar; avisos vao no resultado.

    `previous`: resultado de uma execucao anterior. Um estagio e reaproveitado se o contexto e o
    mesmo, e se ele tem o mesmo id, o mesmo bloco, as mesmas entradas e nenhuma entrada recalculada.
    """
    violations = validate(system, simulation)
    if any(v.severity is Severity.ERROR for v in violations):
        raise ValidationFailed(violations)
    resolved = resolve(system, simulation)
    pipeline = build(system, simulation, resolved)
    ctx = SimContext(grid=resolved.grid, root_seed=simulation.root_seed)
    before = ({s.id: s for s in previous.stages}
              if previous is not None and previous.context == ctx else {})
    outputs: dict[str, Any] = {}
    recalculated: set[str] = set()
    stages: list[Stage] = []
    for sid, block, inputs in pipeline:
        old = before.get(sid)
        if old is not None and old.block == block and old.inputs == inputs \
                and not recalculated.intersection(inputs):
            out, dt, reused = old.output, 0.0, True
        else:
            if not inputs:
                x = None
            elif len(inputs) == 1:
                x = outputs[inputs[0]]
            else:
                x = [outputs[i] for i in inputs]
            t0 = time.perf_counter()
            out = block(x, ctx)
            dt, reused = time.perf_counter() - t0, False
            recalculated.add(sid)
        outputs[sid] = out
        stages.append(Stage(sid, block, out, dt, reused, inputs))
    metrics = _measure(system, resolved, stages)
    return Result(system, simulation, resolved, ctx, stages, metrics, violations)

def run_file(path: str, previous: Result | None = None) -> Result:
    """Carrega um cenario em YAML ou JSON e roda."""
    system, simulation = load(path)
    return run(system, simulation, previous)

# --------------------------------------------------------------------------- varredura ---------
_TOKEN = re.compile(r"^(\w+)(?:\[(\d+)\])?$")

def _replace_path(obj, tokens: list[str], value):
    if not tokens:
        return value
    m = _TOKEN.match(tokens[0])
    if m is None:
        raise ValueError(f"trecho de caminho invalido: {tokens[0]!r}")
    name, idx = m.group(1), m.group(2)
    if not hasattr(obj, name):
        raise ValueError(f"{type(obj).__name__} nao tem o campo {name!r}")
    child = getattr(obj, name)
    if idx is None:
        new = _replace_path(child, tokens[1:], value)
    else:
        items = list(child)
        items[int(idx)] = _replace_path(items[int(idx)], tokens[1:], value)
        new = tuple(items)
    return dataclasses.replace(obj, **{name: new})

def with_value(system: SystemConfig, simulation: SimulationConfig, parameter: str, value):
    """(system, simulation) com o campo `parameter` trocado por `value`.

    Caminhos relativos ao sistema ("link.spans[0].fiber.length", "channels[0].launch_power") ou,
    com o prefixo "simulation.", a simulacao ("simulation.root_seed"). As regras da configuracao
    (__post_init__) voltam a valer para o valor novo.
    """
    tokens = parameter.split(".")
    if tokens[0] == "simulation":
        return system, _replace_path(simulation, tokens[1:], value)
    if tokens[0] == "system":
        tokens = tokens[1:]
    return _replace_path(system, tokens, value), simulation

@dataclass
class Sweep:
    parameter: str
    values: np.ndarray
    metrics: list[dict[int, QMeasurement]]
    received_power: list[dict[int, float]]          # W, na entrada do fotodetector
    results: list[Result] | None = None              # so com keep_results=True

    def metric(self, name: str, channel: int = 0) -> np.ndarray:
        """Uma metrica (campo de QMeasurement, por exemplo "q" ou "ber_q") ao longo da varredura."""
        return np.array([getattr(m[channel], name) for m in self.metrics])

    def to_dict(self) -> dict:
        return {"parameter": self.parameter, "values": self.values.tolist(),
                "metrics": [{str(ch): {**dataclasses.asdict(q),
                                       "received_power_dbm": w_to_dbm(p[ch])}
                             for ch, q in m.items()}
                            for m, p in zip(self.metrics, self.received_power)]}

def sweep(system: SystemConfig, simulation: SimulationConfig, parameter: str, values,
          keep_results: bool = False,
          on_point: Callable[[int, float, Result], None] | None = None) -> Sweep:
    """Roda o sistema para cada valor de `parameter`, reaproveitando os estagios que nao mudam.

    Os ruidos e os bits sao os mesmos em todos os pontos (mesmas sementes): a varredura compara
    pontos com a mesma realizacao, nao faz media de Monte Carlo.
    """
    values = np.asarray(values, dtype=float)
    metrics, powers, results = [], [], []
    previous = None
    for i, v in enumerate(values):
        s, sim = with_value(system, simulation, parameter, float(v))
        r = run(s, sim, previous)
        metrics.append(r.metrics)
        powers.append({ch: r.received_power(ch) for ch in r.metrics})
        if keep_results:
            results.append(r)
        if on_point is not None:
            on_point(i, float(v), r)
        previous = r
    return Sweep(parameter, values, metrics, powers, results if keep_results else None)
