"""Compute X-H local-mode overtones (GFN2-xTB, ALPB water) for library molecules.

    python run_xh.py                 # every molecule in xh_molecules.py
    python run_xh.py fructose allulose
Writes results_xh/<key>.json (one file per molecule, all forms).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from sugar_nir.xh import build_form, compute_form
from xh_molecules import MOLECULES

OUT = Path("results_xh")


def main(keys):
    OUT.mkdir(exist_ok=True)
    for key in keys:
        m = MOLECULES[key]
        t0 = time.time()
        res = {"key": key, "name": m["name"], "en": m["en"], "kind": m["kind"], "forms": {}}
        tot = sum(f["pop"] for f in m["forms"].values())
        for fname, spec in m["forms"].items():
            mol = build_form(spec)
            n_embed = spec.get("n_embed", m.get("n_embed", 100))
            r = compute_form(mol, n_embed=n_embed, keep=5, name=f"{key} {fname}",
                             log=lambda s: print(s, flush=True))
            r["pop"] = spec["pop"] / tot
            res["forms"][fname] = r
            json.dump(res, open(OUT / f"{key}.json", "w"), ensure_ascii=False)
        print(f"{key} done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(MOLECULES))
