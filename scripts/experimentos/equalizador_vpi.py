"""Hipotese 2 das bordas (notebook 09): o equalizador do VPI muda o que a busca de fase recebe?

O MIMO do VPI (`TDE_MIMO`, `AlgTDE = MMA`, `Iterations = 20`, `StepSize = 1e-3`) se adapta por
gradiente, simbolo a simbolo, pelo raio dos pontos (criterio dirigido por raio). O receptor cego do
simcopy usa o modulo constante e depois um ajuste por minimos quadrados com as decisoes, em bloco.
Um equalizador por gradiente deixa os coeficientes oscilando em torno da solucao, e em potencia
alta o criterio do raio ve o ruido de fase nao linear de outro jeito. As duas coisas podem mudar
onde a busca de fase se perde.

Cada equalizador e combinado com tres buscas de fase: a de hoje (conjunta, decisao no 16QAM) e as
duas melhores da hipotese 1 (dois estagios por polarizacao, so os aneis de dentro e de fora no
estagio grosso, janela grossa de 32 e de 21). Equalizadores:
  simcopy             o receptor cego de hoje (modulo constante e 4 ajustes por decisao em bloco)
  raio_P_corrente     Stokes, modulo constante (2 passagens) e o criterio por raio por 20 passagens,
                      passo P, saida com os coeficientes de cada instante (como o VPI processa)
  raio_P_congelada    o mesmo, com os coeficientes do fim aplicados ao bloco inteiro
  raio_lms_corrente   o criterio por raio seguido de 4 passagens dirigidas por decisao (LMS, 1e-3)
O passo vale para sinais de potencia unitaria (o do VPI nao se compara: a normalizacao dele nao esta
no documento).

Uso, na pasta repo:
  python scripts/experimentos/equalizador_vpi.py PASTA_SINAIS SAIDA.json [--files 125km_6dBm ...]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import busca_fase_vpi as bf          # noqa: E402  (troca dsp.carrier_recovery pela versao com variantes)
import equalizador_gradiente as eg   # noqa: E402
from simcopy.kernels import dsp      # noqa: E402
from simcopy.kernels.demodulation import measure_blind  # noqa: E402
from simcopy.kernels.mapping import map_bits  # noqa: E402

MODULATION = "16QAM"
WINDOW = 43
SEARCHES = {
    "conjunta": dict(),
    "dois_estagios_32": dict(joint=False, fourth_power_span_deg=22.5, partition=True, coarse_window=32),
    "dois_estagios_21": dict(joint=False, fourth_power_span_deg=22.5, partition=True, coarse_window=21),
}
DEFAULTS = {"grid_scale": None, "joint": True, "fourth_power_span_deg": None, "coarse_window": None,
            "partition": False}


def use_search(name: str) -> None:
    bf.STATE.update({**DEFAULTS, **SEARCHES[name]})


def radius_outputs(x: np.ndarray, step: float, passes: int = 20) -> dict[str, np.ndarray]:
    """Saidas do equalizador por raio, antes da recuperacao de portadora: corrente e congelada."""
    u, w, radii2 = eg.cma_start(x)
    w = dsp._adapt(u, w, step, radii2, passes - 1)
    running, w = eg.running_pass(u, w, step, radii2)
    frozen = (u @ w.T).T
    norm = lambda y: y / np.sqrt(np.mean(np.abs(y) ** 2, axis=1, keepdims=True))
    return {"corrente": norm(running), "congelada": norm(frozen)}


def record(meta: dict, name: str, z: np.ndarray, bits, guard: int, a: np.ndarray, seconds: float) -> dict:
    m, _ = measure_blind(z, bits, MODULATION, guard)
    return {**meta, "variant": name, "ber": float(m.ber_counted), "ser": float(m.ser),
            "snr_db": float(m.snr_db), "phase_slips": int(m.phase_slips),
            "rotation": bf.rotation(z, a, guard), "seconds": seconds}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--files", nargs="*", default=None)
    ap.add_argument("--steps", type=float, nargs="*", default=[1e-3, 3e-4])
    ap.add_argument("--searches", nargs="*", default=list(SEARCHES))
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    results = json.loads(args.out.read_text(encoding="utf-8")) if args.out.exists() else []
    done = {(r["file"], r["variant"]) for r in results}
    files = sorted(args.folder.glob("*.npz"))
    if args.files:
        files = [f for f in files if any(f.name.startswith(p + "_") for p in args.files)]

    def save(r):
        results.append(r)
        args.out.write_text(json.dumps(results), encoding="utf-8")
        g = r["rotation"]
        print(f"{r['file'][:-4]:22s} {r['variant']:42s}: BER {r['ber']:.2e} | SER {r['ser']:.2e} | "
              f"SNR {r['snr_db']:6.2f} | girado {g['rotated_fraction']:.3f} | {r['seconds']:.0f} s", flush=True)

    for f in files:
        d = np.load(f)
        x = d["adc"].astype(complex)
        bits = d["bits"]
        a = np.array([map_bits(row, MODULATION) for row in bits])
        guard = int(d["guard"])
        meta = {"file": f.name, "length_km": float(d["length_km"]), "power_dbm": float(d["power_dbm"])}
        for search in args.searches:
            name = f"simcopy+{search}"
            if (f.name, name) not in done:
                use_search(search)
                t0 = time.perf_counter()
                z, _ = dsp.blind_receiver(x, MODULATION, window=WINDOW)
                save(record(meta, name, z, bits, guard, a, time.perf_counter() - t0))
        for step in args.steps:
            names = {(form, search): f"raio_{step:g}_{form}+{search}"
                     for form in ("corrente", "congelada") for search in args.searches}
            if all((f.name, n) in done for n in names.values()):
                continue
            t0 = time.perf_counter()
            outs = radius_outputs(x, step)
            for (form, search), name in names.items():
                if (f.name, name) in done:
                    continue
                use_search(search)
                z = dsp.carrier_recovery(outs[form], MODULATION, 64, WINDOW)[0]
                save(record(meta, name, z, bits, guard, a, time.perf_counter() - t0))
        for search in args.searches:
            name = f"raio_lms_corrente+{search}"
            if (f.name, name) in done:
                continue
            use_search(search)
            t0 = time.perf_counter()
            z = eg.dd_lms_variants(x, 1e-3, 1e-3, 20, 4)["corrente"]
            save(record(meta, name, z, bits, guard, a, time.perf_counter() - t0))


if __name__ == "__main__":
    main()
