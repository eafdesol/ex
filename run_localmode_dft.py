"""DFT check of key C-H local modes (B3LYP/6-31G**, ddCOSMO water).

Maltose itself (45 atoms) is too slow for DFT scans here (> 8 min per point),
so its two glycosidic C-H environments are checked on small model compounds:
  me-a-glucoside  C1-H : H-C(O5)(O-C)  = maltose non-reducing C1'-H
  4-OMe-b-glucose C4-H : H-C-O-C       = maltose reducing-unit C4-H
plus reference C-H of alpha/beta-glucose. xTB is re-run on the same geometry
and grid so the two methods can be compared one-to-one.

    python run_localmode_dft.py
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

from sugar_nir.localmode import (XTB, conformers, label_atoms, scan_ch, smiles_of,
                                 solve_local_mode)

OUT = Path("results_localmode")
QD = np.linspace(-0.25, 0.55, 11)
TARGETS = {
    "b-glucose": ["C1-H", "C2-H", "C4-H", "C5-H", "C6-Ha"],
    "a-glucose": ["C1-H", "C2-H"],
    "me-a-glucoside": ["C1-H"],
    "4-OMe-b-glucose": ["C4-H"],
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
        mf = mf.ddCOSMO()
        mf.with_solvent.eps = 78.36
        e = mf.kernel(dm0=self.dm)
        self.dm = mf.make_rdm1()
        return {"energy": e, "dipole": mf.dip_moment(unit="AU", verbose=0)}


def geometry(name, geoms):
    if name in geoms:
        return np.array(geoms[name]["numbers"]), np.array(geoms[name]["pos"][0])
    mol, labels, numbers, confs, w, st = conformers(name, n_embed=80, keep=1)
    return numbers, confs[0]


def main():
    geoms = json.load(open(OUT / "geometries.json"))
    rows = []
    for name, targets in TARGETS.items():
        numbers, pos = geometry(name, geoms)
        labels = label_atoms(Chem.AddHs(Chem.MolFromSmiles(smiles_of(name))))
        for h, d in labels.items():
            if d["label"] not in targets:
                continue
            row = dict(molecule=name, label=d["label"], cls=d["cls"])
            for meth, eng in (("xtb", XTB(numbers)), ("dft", DFT(numbers))):
                t0 = time.time()
                E, D = scan_ch(eng, pos, h, d["carbon"], Q=QD)
                nu, I, we, wexe = solve_local_mode(E, D, Q=QD)
                row.update({f"{meth}_nu1": nu[0], f"{meth}_nu2": nu[1], f"{meth}_nu3": nu[2],
                            f"{meth}_we": we, f"{meth}_wexe": wexe,
                            f"{meth}_npts": int(np.isfinite(E).sum())})
                print(f"{name:16s} {d['label']:6s} {meth}: nu1={nu[0]:7.1f} nu2={nu[1]:7.1f} "
                      f"we={we:6.0f} wexe={wexe:5.1f} ({time.time() - t0:.0f} s)", flush=True)
            rows.append(row)
            with open(OUT / "ch_overtones_dft.csv", "w", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=list(rows[0]))
                wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
