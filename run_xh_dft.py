"""DFT (B3LYP/6-31G**, ddCOSMO water) check scans for selected X-H bonds of
library molecules, paired with xTB on the same geometry and grid. Used to
calibrate xTB per bond type (C-H, O-H).

    python run_xh_dft.py oh      # O-H set
    python run_xh_dft.py ch5     # furanose C-H set
Writes results_xh/dft_pairs_<set>.csv (appended per bond, so partial runs keep results).
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

from run_localmode_dft import DFT, QD
from sugar_nir.localmode import XTB, scan_ch, solve_local_mode
from sugar_nir.xh import reduced_mass

OUT = Path("results_xh")
SETS = {
    "oh": [("glucose", "β-피라노스", "O1-H"), ("glucose", "β-피라노스", "O2-H"),
           ("glucose", "β-피라노스", "O3-H"), ("glucose", "β-피라노스", "O6-H"),
           ("glucose", "α-피라노스", "O1-H"), ("fructose", "β-피라노스", "O2-H")],
    "ch5": [("fructose", "β-퓨라노스", "C3-H"), ("fructose", "β-퓨라노스", "C5-H"),
            ("fructose", "β-피라노스", "C3-H")],
}


def parse_xyz(xyz):
    from rdkit import Chem
    lines = xyz.strip().splitlines()[2:]
    pt = Chem.GetPeriodicTable()
    numbers = np.array([pt.GetAtomicNumber(l.split()[0]) for l in lines])
    pos = np.array([[float(v) for v in l.split()[1:4]] for l in lines])
    return numbers, pos


def main(set_name):
    path = OUT / f"dft_pairs_{set_name}.csv"
    done = set()
    if path.exists():
        done = {(r["molecule"], r["form"], r["label"]) for r in csv.DictReader(open(path))}
    for key, form, label in SETS[set_name]:
        if (key, form, label) in done:
            continue
        f = json.load(open(OUT / f"{key}.json"))["forms"][form]
        numbers, pos = parse_xyz(f["xyz"])
        b = next(r for r in f["bonds"] if r["label"] == label)
        row = dict(molecule=key, form=form, label=label, elem=b["elem"], cls=b["cls"])
        for meth, eng in (("xtb", XTB(numbers)), ("dft", DFT(numbers))):
            t0 = time.time()
            E, D = scan_ch(eng, pos, b["h"], b["x"], Q=QD)
            nu, I, we, wexe = solve_local_mode(E, D, Q=QD, mu=reduced_mass(b["elem"]))
            row.update({f"{meth}_nu1": nu[0], f"{meth}_nu2": nu[1], f"{meth}_nu3": nu[2],
                        f"{meth}_I2": I[1], f"{meth}_npts": int(np.isfinite(E).sum())})
            print(f"{key} {form} {label} {meth}: nu2={nu[1]:.1f} ({time.time() - t0:.0f} s)", flush=True)
        new = not path.exists()
        with open(path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row))
            if new:
                w.writeheader()
            w.writerow(row)


if __name__ == "__main__":
    main(sys.argv[1])
