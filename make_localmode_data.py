"""Export local-mode results as JSON for the HTML viewer (results_localmode/viewer_data.json)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

import analyze_localmode as A
from sugar_nir.localmode import SMILES, label_atoms

OUT = Path("results_localmode")


def dih(p, a, b, c, d):
    b0, b1, b2 = p[a] - p[b], p[c] - p[b], p[d] - p[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def main():
    geoms = json.load(open(OUT / "geometries.json"))
    raw = {(r["molecule"], r["label"]): r for r in A.load("xtb")}
    cal = {(r["molecule"], r["label"]): r for r in A.load("cal")}
    ori = A.h_orientation()
    dft = {(r["molecule"], r["label"]): r for r in A.dft_pairs()}
    calib, pairs = A.calibration()
    mols, rows = {}, []
    for name in SMILES:
        mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
        labs = label_atoms(mol)
        P = [np.array(x) for x in geoms[name]["pos"]]
        W = np.array(geoms[name]["weights"])
        sym = [a.GetSymbol() for a in mol.GetAtoms()]
        mols[name] = {"xyz": f"{len(sym)}\n{name}\n" + "\n".join(
            f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, (x, y, z) in zip(sym, P[0])),
            "n_conf": len(P)}
        for h, d in labs.items():
            c = d["carbon"]
            lp, dihs = 0.0, []
            for o in [n for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetSymbol() == "O"]:
                X = [n.GetIdx() for n in o.GetNeighbors() if n.GetIdx() != c][0]
                tag = "O5" if o.GetIdx() == d["o5"] else (
                    "O-H" if mol.GetAtomWithIdx(X).GetSymbol() == "H" else "O-C")
                ds = [dih(p, h, c, o.GetIdx(), X) for p in P]
                lp += float(np.average([1.0 if 25 < abs(x) < 95 else 0.0 for x in ds], weights=W))
                dihs.append({"o": tag, "deg": round(ds[0])})
            k = (name, d["label"])
            dd = dft.get(k)
            rows.append({
                "mol": name, "label": d["label"], "cls": d["cls"], "orient": ori[k],
                "h": int(h), "c": int(c),
                "r": round(float(np.average([np.linalg.norm(p[h] - p[c]) for p in P], weights=W)), 4),
                "lp": round(lp, 1), "dih": dihs,
                "xtb2": round(raw[k]["nu2"], 1), "cal2": round(cal[k]["nu2"], 1),
                "cal3": round(cal[k]["nu3"], 1), "I2": raw[k]["I2"],
                "dft2": round(float(dd["dft_nu2"]), 1) if dd else None,
            })
    out = {
        "molecules": mols, "ch": rows,
        "pop": A.POP, "mw": A.MW,
        "calib": {"a": calib[2][0], "b": calib[2][1], "rmse": calib[2][2], "n": calib[2][3]},
        "pairs": [{"mol": p["molecule"], "label": p["label"], "cls": p["cls"],
                   "xtb": round(float(p["xtb_nu2"]), 1), "dft": round(float(p["dft_nu2"]), 1)}
                  for p in pairs],
    }
    json.dump(out, open(OUT / "viewer_data.json", "w"), ensure_ascii=False)
    print(len(rows), "C-H rows;", sum(len(m["xyz"]) for m in mols.values()), "bytes of xyz")


if __name__ == "__main__":
    main()


def build_html():
    data = (OUT / "viewer_data.json").read_text(encoding="utf-8")
    tpl = Path("tools/ch_viewer_template.html").read_text(encoding="utf-8")
    (OUT / "ch_overtone_viewer.html").write_text(tpl.replace("/*DATA*/null", data), encoding="utf-8")
    print("wrote", OUT / "ch_overtone_viewer.html")


if __name__ == "__main__":
    build_html()
