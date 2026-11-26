"""C-H local-mode overtones of glucose (a/b) and maltose (a/b reducing end).

    python run_localmode.py            # GFN2-xTB / ALPB(water)
Outputs: results_localmode/ch_overtones.csv, report.md, figures
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np

from sugar_nir.localmode import Q, SMILES, XTB, conformers, scan_ch, solve_local_mode

OUT = Path("results_localmode")
N_EMBED = {"a-glucose": 80, "b-glucose": 80, "a-maltose": 400, "b-maltose": 400}


def main():
    OUT.mkdir(exist_ok=True)
    rows, geoms = [], {}
    for name in SMILES:
        t0 = time.time()
        mol, labels, numbers, confs, w, st = conformers(name, n_embed=N_EMBED[name])
        print(f"{name}: {len(confs)} conformers, weights {np.round(w, 2)}", flush=True)
        geoms[name] = {"numbers": numbers.tolist(), "pos": [c.tolist() for c in confs],
                       "weights": w.tolist()}
        eng = XTB(numbers)
        for h, d in sorted(labels.items(), key=lambda t: t[1]["label"]):
            res = []
            for pos in confs:
                E, D = scan_ch(eng, pos, h, d["carbon"])
                if np.isfinite(E).sum() < Q.size:
                    print(f"    {name} {d['label']}: scan q in [{Q[np.isfinite(E)][0]:.2f}, "
                          f"{Q[np.isfinite(E)][-1]:.2f}] A", flush=True)
                nu, I, we, wexe = solve_local_mode(E, D)
                rch = np.linalg.norm(pos[h] - pos[d["carbon"]])
                res.append((nu, I, we, wexe, rch))
            nu = np.average([r[0] for r in res], axis=0, weights=w)
            I = np.average([r[1] for r in res], axis=0, weights=w)
            rows.append(dict(molecule=name, label=d["label"], cls=d["cls"], unit=d["unit"],
                             we=float(np.average([r[2] for r in res], weights=w)),
                             wexe=float(np.average([r[3] for r in res], weights=w)),
                             r_CH=float(np.average([r[4] for r in res], weights=w)),
                             nu1=nu[0], nu2=nu[1], nu3=nu[2], I1=I[0], I2=I[1], I3=I[2],
                             nu2_sd=float(np.std([r[0][1] for r in res]))))
        print(f"  done in {time.time() - t0:.0f} s", flush=True)
    with open(OUT / "ch_overtones_xtb.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    json.dump(geoms, open(OUT / "geometries.json", "w"))


if __name__ == "__main__":
    main()
