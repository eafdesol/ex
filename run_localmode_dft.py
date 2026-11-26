"""DFT check of selected C-H local modes (B3LYP/6-31G**, ddCOSMO water).

Uses the lowest xTB conformer of each molecule (results_localmode/geometries.json)
and repeats the C-H stretch scan with DFT for the requested bonds.

    python run_localmode_dft.py b-glucose b-maltose
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

from sugar_nir.localmode import Q, SMILES, label_atoms, scan_ch, solve_local_mode

OUT = Path("results_localmode")
TARGET = {   # which C-H to check (all of glucose, the informative ones of maltose)
    "b-glucose": None,
    "a-glucose": None,
    "b-maltose": ["non-red. C1-H", "red. C4-H", "red. C1-H", "red. C2-H", "non-red. C2-H",
                  "red. C5-H", "non-red. C5-H"],
}


class DFT:
    def __init__(self, numbers, xc="b3lyp", basis="6-31g**"):
        self.numbers, self.xc, self.basis = numbers, xc, basis
        self.dm = None

    def __call__(self, pos_A, grad=False):
        from pyscf import dft, gto
        sym = [Chem.GetPeriodicTable().GetElementSymbol(int(z)) for z in self.numbers]
        mol = gto.M(atom=list(zip(sym, pos_A)), basis=self.basis, unit="Angstrom", verbose=0)
        mf = dft.RKS(mol).density_fit()
        mf.xc = self.xc
        mf.grids.level = 3
        mf = mf.ddCOSMO()
        mf.with_solvent.eps = 78.36
        e = mf.kernel(dm0=self.dm)
        self.dm = mf.make_rdm1()
        return {"energy": e, "dipole": mf.dip_moment(unit="AU", verbose=0)}


def main(names):
    geoms = json.load(open(OUT / "geometries.json"))
    rows = []
    for name in names:
        g = geoms[name]
        numbers = np.array(g["numbers"])
        pos = np.array(g["pos"][0])
        mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
        labels = label_atoms(mol)
        eng = DFT(numbers)
        for h, d in sorted(labels.items(), key=lambda t: t[1]["label"]):
            if TARGET.get(name) and d["label"] not in TARGET[name]:
                continue
            t0 = time.time()
            eng.dm = None
            E, D = scan_ch(eng, pos, h, d["carbon"])
            nu, I, we, wexe = solve_local_mode(E, D)
            rows.append(dict(molecule=name, label=d["label"], cls=d["cls"], we=we, wexe=wexe,
                             nu1=nu[0], nu2=nu[1], nu3=nu[2], I1=I[0], I2=I[1], I3=I[2]))
            print(f"{name} {d['label']:16s} nu1={nu[0]:7.1f} nu2={nu[1]:7.1f} "
                  f"({time.time() - t0:.0f} s)", flush=True)
            with open(OUT / "ch_overtones_dft.csv", "w", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=list(rows[0]))
                wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1:] or ["b-glucose", "b-maltose"])
