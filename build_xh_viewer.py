"""Collect results_xh/*.json, calibrate xTB -> DFT per bond type, and build the
HTML viewer (results_xh/xh_overtone_viewer.html).

    python build_xh_viewer.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from sugar_nir.simulate import water_real
from xh_molecules import MOLECULES

OUT = Path("results_xh")
TEMPLATE = Path("tools/xh_viewer_template.html")


def load_pairs():
    pairs = []
    files = [Path("results_localmode/ch_overtones_dft.csv"),
             Path("results_localmode/ch_overtones_dft_4-OMe-b-glucose.csv")]
    files += sorted(OUT.glob("dft_pairs_*.csv"))
    for f in files:
        if not f.exists():
            continue
        for r in csv.DictReader(open(f, encoding="utf-8")):
            elem = r.get("elem") or "C"
            pairs.append({"mol": r["molecule"], "form": r.get("form", ""), "label": r["label"],
                          "elem": elem, "cls": r.get("cls", ""),
                          **{f"{m}{v}": float(r[f"{m}_nu{v}"]) for m in ("xtb", "dft") for v in (1, 2, 3)}})
    return pairs


def calibrate(pairs):
    cal = {}
    for elem in sorted({p["elem"] for p in pairs}):
        P = [p for p in pairs if p["elem"] == elem]
        if len(P) < 3:
            continue
        cal[elem] = {}
        for v in (1, 2, 3):
            x = np.array([p[f"xtb{v}"] for p in P])
            y = np.array([p[f"dft{v}"] for p in P])
            b, a = np.polyfit(x, y, 1)
            rmse = float(np.sqrt(np.mean((y - (a + b * x)) ** 2)))
            cal[elem][str(v)] = {"a": float(a), "b": float(b), "rmse": rmse, "n": len(P)}
    return cal


def main():
    pairs = load_pairs()
    cal = calibrate(pairs)
    mols = {}
    for key in MOLECULES:
        f = OUT / f"{key}.json"
        if not f.exists():
            continue
        m = json.load(open(f, encoding="utf-8"))
        for form in m["forms"].values():
            for b in form["bonds"]:
                c = cal.get(b["elem"])
                for v in (1, 2, 3):
                    b[f"cal{v}"] = (round(c[str(v)]["a"] + c[str(v)]["b"] * b[f"nu{v}"], 1)
                                    if c else None)
                b["rel"] = ("보통" if b["elem"] == "C" and c else
                            "낮음 · 물과의 수소결합 미반영" if b["elem"] == "O" and c else
                            "낮음 · DFT 보정 없음")
        m["mw"] = next(iter(m["forms"].values()))["mw"]
        mols[key] = m
    nm = np.arange(900.0, 2502.0, 2.0)
    alpha = water_real(nm, 25.0, path_mm=10.0) * np.log(10)      # Napierian cm^-1
    data = {"molecules": mols, "calib": cal, "pairs": pairs,
            "water": {"nm": nm.tolist(), "alpha": [round(float(a), 5) for a in alpha]}}
    (OUT / "xh_library.json").write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    html = TEMPLATE.read_text(encoding="utf-8").replace("/*DATA*/null", json.dumps(data, ensure_ascii=False))
    (OUT / "xh_overtone_viewer.html").write_text(html, encoding="utf-8")
    print("molecules:", list(mols), "| calibration:",
          {e: {v: (round(c[v]["a"]), round(c[v]["b"], 3), round(c[v]["rmse"]), c[v]["n"])
               for v in c} for e, c in cal.items()})


if __name__ == "__main__":
    main()
