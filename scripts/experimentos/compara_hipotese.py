"""Placar de uma hipotese contra o VPI, ponto a ponto, nos conjuntos de configuracao conhecida.

Le um arquivo de resultados de experimento (lista de registros com variant, length_km, power_dbm,
ber, ser e rotation.rotated_fraction) e o catalogo do VPI (reference/vpi/catalogo/catalogo_vpi.json).
Para cada variante:
  bordas   quantos pontos em que o VPI perde a fase a variante tambem perde
  meio     quantos pontos em que o VPI converge a variante tambem converge
  BER      a BER da variante sobre a do VPI, onde os dois convergem
  giro     a fracao do bloco girada, da variante e do VPI, onde os dois perdem a fase
"A fase se perdeu" = mais de 2 % do bloco girado, ou BER acima de 0,1 (o mesmo criterio do notebook 09).

Uso, na pasta repo:
  python scripts/experimentos/compara_hipotese.py RESULTADOS.json [--json RESUMO.json]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def lost(ber: float, rotated: float) -> bool:
    return rotated > 0.02 or ber > 0.1


def vpi_points() -> dict:
    points = json.loads((ROOT / "reference/vpi/catalogo/catalogo_vpi.json").read_text(encoding="utf-8"))["points"]
    out = {}
    for p in points:
        if p["set"] == "2026_4ch" or (p["set"] == "400gbps_CH2" and p["length_km"] == 120):
            out[(int(p["length_km"]), int(p["power_dbm"]))] = p
    return out


def score(records: list[dict], vpi: dict) -> dict:
    edge_hit, edge_n, mid_hit, mid_n, ratios, rot = 0, 0, 0, 0, [], []
    misses = []
    for r in records:
        key = (int(r["length_km"]), int(r["power_dbm"]))
        if key not in vpi:
            continue
        v = vpi[key]
        rotated = r["rotation"]["rotated_fraction"] if "rotation" in r else r.get("rotated_fraction", 0.0)
        mine, theirs = lost(r["ber"], rotated), lost(v["ber"], v["rotated_fraction"])
        if theirs:
            edge_n += 1
            edge_hit += mine
            if mine:
                rot.append((key, rotated, v["rotated_fraction"], r["ber"], v["ber"]))
            else:
                misses.append(f"{key[0]}/{key[1]} nao perde")
        else:
            mid_n += 1
            mid_hit += not mine
            if mine:
                misses.append(f"{key[0]}/{key[1]} perde a mais")
            else:
                ratios.append(r["ber"] / v["ber"])
    ratios = np.array(ratios) if ratios else np.array([np.nan])
    return {"edges": f"{edge_hit}/{edge_n}", "middle": f"{mid_hit}/{mid_n}",
            "ber_ratio_median": float(np.median(ratios)), "ber_ratio_min": float(np.min(ratios)),
            "ber_ratio_max": float(np.max(ratios)), "both_lost": rot, "misses": misses,
            "edge_hit": edge_hit, "edge_n": edge_n, "mid_hit": mid_hit, "mid_n": mid_n}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results", type=Path)
    ap.add_argument("--json", type=Path, default=None, help="grava o placar resumido")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    records = json.loads(args.results.read_text(encoding="utf-8"))
    vpi = vpi_points()
    variants = list(dict.fromkeys(r["variant"] for r in records))
    summary = {}
    print(" variante                                  | bordas | meio  | BER/VPI no meio (mediana, faixa) | erros")
    for name in variants:
        s = score([r for r in records if r["variant"] == name], vpi)
        summary[name] = s
        print(f" {name:41s} | {s['edges']:>6} | {s['middle']:>5} | {s['ber_ratio_median']:.2f} "
              f"({s['ber_ratio_min']:.2f} a {s['ber_ratio_max']:.2f})                | {', '.join(s['misses'][:6])}")
    print("\n giro (fracao do bloco) e BER onde os dois perdem a fase: variante contra VPI")
    for name in variants:
        both = summary[name]["both_lost"]
        if both:
            cells = [f"{k[0]}/{k[1]}: {a:.2f} x {b:.2f}, BER {c:.2f} x {d:.2f}" for k, a, b, c, d in both]
            print(f" {name}: " + " | ".join(cells))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")


if __name__ == "__main__":
    main()
