"""Experimento da hipotese C3 (docs/plano_replicar_vpi.md): equalizador adaptado por gradiente,
simbolo a simbolo, como o `TDE_MIMO` do VPI (`AlgTDE = MMA`, `Iterations = 20`), sobre os sinais
guardados na saida do ADC. Fica fora de src/: so entra no simulador se a hipotese se confirmar.

A pergunta: em potencia baixa o VPI tem cerca de 5 % mais ruido do receptor que o receptor cego
do simcopy, e um pouco menos interferencia entre simbolos (docs/replicacao_vpi_diario.md, linha
A). Um equalizador por gradiente oscila em torno da solucao de minimo erro, e esse desajuste soma
ruido. O passo do VPI (`StepSize = 1e-3`) nao se compara ao do simcopy, porque a normalizacao do
sinal no VPI nao esta no documento: aqui o passo e varrido.

Cadeia: separacao pelo espaco de Stokes, modulo constante (2 passagens, passo 1e-3), criterio
dirigido por raio (o MMA com os tres raios do 16QAM) por `passes` passagens com o passo `step`, e a
recuperacao de portadora do simcopy (busca cega de fase com janela de 43). Sem o ajuste por decisao
em bloco do simcopy. Duas formas de saida, da mesma adaptacao:
  congelada  os coeficientes do fim da ultima passagem aplicados ao bloco inteiro (como o simcopy);
  corrente   cada simbolo sai com os coeficientes daquele instante, na ultima passagem (como um
             equalizador por gradiente que processa o sinal enquanto se adapta).
A referencia e o receptor cego do simcopy com janela de 43 (`base`).

Com --dd-steps, mais uma variante: depois do criterio por raio (passo --mma-step), --dd-passes
passagens dirigidas por decisao por gradiente (LMS), o analogo do ajuste em bloco do simcopy (o
bloco do VPI tambem lista `Mode = Decision-Directed`).

O `step` vale para sinais de potencia unitaria (como em kernels/dsp.py, `adaptive_equalizer`).

Uso, na pasta repo:
  python scripts/experimentos/equalizador_gradiente.py PASTA_SINAIS SAIDA.json
      [--files 125km_8dBm 125km_13dBm] [--steps 3e-4 1e-3 3e-3] [--passes 20]
      [--vpi reference/vpi/replicacao/etapa_a_vpi.json]
Os sinais sao os de `scripts/compara_vpi.py simular --save-signals`. O SAIDA.json e retomado: o
que ja esta nele nao roda de novo.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from medidas_constelacao import measure, summary  # noqa: E402
from simcopy.kernels import dsp  # noqa: E402
from simcopy.kernels.demodulation import align_to_transmitted  # noqa: E402
from simcopy.kernels.mapping import constellation, map_bits  # noqa: E402

MODULATION = "16QAM"
TAPS = 7
WINDOW = 43
CMA_PASSES = 2
CMA_STEP = 1e-3


def cma_start(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stokes, potencia unitaria, regressores e coeficientes depois do modulo constante (com a
    mesma protecao de `adaptive_equalizer` contra as duas saidas na mesma polarizacao)."""
    x = dsp.stokes_polarization_demux(x)
    x = x / np.sqrt(np.mean(np.abs(x) ** 2, axis=1, keepdims=True))
    u = dsp._regressors(x, TAPS)
    r2, radii2 = dsp.modulus_radii(constellation(MODULATION))
    w = np.zeros((2, 2 * TAPS), dtype=complex)
    w[0, TAPS // 2] = w[1, TAPS + TAPS // 2] = 1.0
    w = dsp._adapt(u, w, CMA_STEP, np.array([r2]), CMA_PASSES)
    jones = w.reshape(2, 2, TAPS).sum(axis=2)
    if abs(np.linalg.det(jones)) < 0.3 * np.linalg.norm(jones[0]) * np.linalg.norm(jones[1]):
        wx = w[0].reshape(2, TAPS)
        w[1] = np.concatenate([-np.conj(wx[1, ::-1]), np.conj(wx[0, ::-1])])
        w = dsp._adapt(u, w, CMA_STEP, np.array([r2]), CMA_PASSES)
    return u, w, radii2


def running_pass(u: np.ndarray, w: np.ndarray, step: float, radii2: np.ndarray):
    """Uma passagem dirigida por raio que guarda a saida de cada simbolo antes de atualizar os
    coeficientes. Devolve (saidas (2, n), coeficientes no fim)."""
    w = w.copy()
    out = np.empty((u.shape[0], 2), dtype=complex)
    for k in range(u.shape[0]):
        uk = u[k]
        y = w @ uk
        out[k] = y
        p = y.real ** 2 + y.imag ** 2
        r = radii2[np.argmin(np.abs(radii2[:, None] - p[None, :]), axis=0)]
        w += (step * (r - p) * y)[:, None] * np.conj(uk)[None, :]
    return out.T, w


def gradient_variants(x: np.ndarray, step: float, passes: int) -> dict[str, np.ndarray]:
    """Simbolos depois da recuperacao de portadora, nas duas formas de saida."""
    u, w, radii2 = cma_start(x)
    w = dsp._adapt(u, w, step, radii2, passes - 1)
    running, w = running_pass(u, w, step, radii2)
    frozen = (u @ w.T).T
    out = {}
    for name, y in (("congelada", frozen), ("corrente", running)):
        y = y / np.sqrt(np.mean(np.abs(y) ** 2, axis=1, keepdims=True))
        out[name] = dsp.carrier_recovery(y, MODULATION, 64, WINDOW)[0]
    return out


def dd_lms_variants(x: np.ndarray, step: float, dd_step: float, passes: int,
                    dd_passes: int) -> dict[str, np.ndarray]:
    """O mesmo inicio de `gradient_variants` (Stokes, modulo constante, `passes` passagens por raio
    com o passo `step`), seguido de `dd_passes` passagens dirigidas por decisao, por gradiente
    (LMS, passo `dd_step`): o alvo de cada simbolo e a decisao levada de volta a fase recebida, como
    no ajuste em bloco de `blind_receiver`, e a portadora e recuperada de novo depois de cada
    passagem. Saidas: congelada e corrente, da ultima passagem."""
    from simcopy.kernels.mapping import nearest_points
    u, w, radii2 = cma_start(x)
    w = dsp._adapt(u, w, step, radii2, passes)
    y = (u @ w.T).T
    z, phases, _ = dsp.carrier_recovery(y / np.sqrt(np.mean(np.abs(y) ** 2, axis=1, keepdims=True)),
                                        MODULATION, 64, WINDOW)
    running = y
    for _ in range(dd_passes):
        target = (nearest_points(z, MODULATION) * np.exp(1j * phases)).T
        out = np.empty((u.shape[0], 2), dtype=complex)
        for k in range(u.shape[0]):
            uk = u[k]
            yk = w @ uk
            out[k] = yk
            w += (dd_step * (target[k] - yk))[:, None] * np.conj(uk)[None, :]
        running = out.T
        z, phases, _ = dsp.carrier_recovery(running, MODULATION, 64, WINDOW)
    frozen = (u @ w.T).T
    return {"congelada": dsp.carrier_recovery(frozen, MODULATION, 64, WINDOW)[0], "corrente": z}


def evaluate(z: np.ndarray, a: np.ndarray, guard: int) -> dict:
    y, slips, _ = align_to_transmitted(z, a, MODULATION)
    y, at = y[:, guard:a.shape[1] - guard], a[:, guard:a.shape[1] - guard]
    s = summary(measure(y, at, False, edge=8))
    return {"snr_db": float(-10 * np.log10(s["inv_snr_total"])), "slips": int(slips), "summary": s}


def vpi_table(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    points = json.loads(path.read_text(encoding="utf-8"))["points"]
    return {(int(p["length_km"]), int(p["power_dbm"])): p["summary"] for p in points}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--files", nargs="*", default=None, help="prefixos, por exemplo 125km_8dBm")
    ap.add_argument("--steps", type=float, nargs="+", default=[3e-4, 1e-3, 3e-3])
    ap.add_argument("--passes", type=int, default=20)
    ap.add_argument("--dd-steps", type=float, nargs="*", default=[],
                    help="passos da etapa por decisao (LMS) depois do criterio por raio com --mma-step")
    ap.add_argument("--mma-step", type=float, default=1e-3)
    ap.add_argument("--dd-passes", type=int, default=4)
    ap.add_argument("--vpi", type=Path, default=ROOT / "reference" / "vpi" / "replicacao" / "etapa_a_vpi.json")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")

    vpi = vpi_table(args.vpi)
    results = json.loads(args.out.read_text(encoding="utf-8")) if args.out.exists() else []
    done = {(r["file"], r["variant"]) for r in results}
    files = sorted(args.folder.glob("*.npz"))
    if args.files:
        files = [f for f in files if any(f.name.startswith(p + "_") for p in args.files)]

    def save(record):
        results.append(record)
        args.out.write_text(json.dumps(results), encoding="utf-8")
        s = record["summary"]
        ref = vpi.get((int(record["length_km"]), int(record["power_dbm"])))
        extra = (f" | VPI: aditivo {ref['additive']:.5f}, interf. {ref['isi']:.5f}, resto {ref['inv_snr_rest']:.5f}"
                 if ref else "")
        print(f"{record['file'][:-4]:22s} {record['variant']:24s}: SNR {record['snr_db']:6.2f} | "
              f"aditivo {s['additive']:.5f} | interf. {s['isi']:.5f} | resto {s['inv_snr_rest']:.5f} | "
              f"saltos {record['slips']} | {record['seconds']:.0f} s{extra}", flush=True)

    for f in files:
        d = np.load(f)
        x = d["adc"].astype(complex)
        a = np.array([map_bits(row, str(d["modulation"])) for row in d["bits"]])
        guard = int(d["guard"])
        meta = {"file": f.name, "length_km": float(d["length_km"]), "power_dbm": float(d["power_dbm"]),
                "seed": int(d["seed"])}
        if (f.name, "base") not in done:
            t0 = time.perf_counter()
            z, _ = dsp.blind_receiver(x, MODULATION, taps=TAPS, window=WINDOW)
            save({**meta, "variant": "base", **evaluate(z, a, guard), "seconds": time.perf_counter() - t0})
        for step in args.steps:
            names = [f"passo{step:g}_{args.passes}x_{form}" for form in ("congelada", "corrente")]
            if all((f.name, n) in done for n in names):
                continue
            t0 = time.perf_counter()
            outs = gradient_variants(x, step, args.passes)
            seconds = time.perf_counter() - t0
            for name, form in zip(names, ("congelada", "corrente")):
                save({**meta, "variant": name, "step": step, "passes": args.passes, "output": form,
                      **evaluate(outs[form], a, guard), "seconds": seconds})
        for dd_step in args.dd_steps:
            names = [f"raio{args.mma_step:g}_lms{dd_step:g}x{args.dd_passes}_{form}"
                     for form in ("congelada", "corrente")]
            if all((f.name, n) in done for n in names):
                continue
            t0 = time.perf_counter()
            outs = dd_lms_variants(x, args.mma_step, dd_step, args.passes, args.dd_passes)
            seconds = time.perf_counter() - t0
            for name, form in zip(names, ("congelada", "corrente")):
                save({**meta, "variant": name, "step": args.mma_step, "dd_step": dd_step,
                      "passes": args.passes, "dd_passes": args.dd_passes, "output": form,
                      **evaluate(outs[form], a, guard), "seconds": seconds})


if __name__ == "__main__":
    main()
