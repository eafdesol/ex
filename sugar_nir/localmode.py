"""C-H local-mode overtone calculation for glucose and maltose.

For every C-H bond:
  1. conformers: RDKit ETKDG + MMFF, re-optimised with GFN2-xTB (ALPB water)
  2. scan the H atom along the C-H bond (q = -0.30 ... +0.70 A), energy and
     dipole at each point (xTB, or DFT for the check)
  3. solve the 1D Schrodinger equation on that potential with the C-H reduced
     mass -> 0->v transition wavenumbers and intensities (no Morse assumption)

This is the standard local-mode / HCAO picture of X-H overtones: the overtone
of each C-H depends on its own anharmonic potential, which is set by its
chemical environment (O-H vs O-C neighbour, axial/equatorial, anomeric...).
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize

BOHR = 0.529177210903          # A
HARTREE_CM = 219474.6313705
B_CONST = 16.857629            # hbar^2/(2 amu A^2) in cm^-1
M_H, M_C = 1.007825, 12.0

SMILES = {   # PubChem isomeric SMILES
    "a-glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
    "b-glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
    "b-maltose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O",
    "a-maltose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](O[C@@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O",
}


_B = SMILES["b-glucose"]
_A = SMILES["a-glucose"]
MODELS = {   # small model compounds for the DFT check of maltose-specific C-H
    # methyl alpha-glucoside: C1-H = H-C(O5)(O-C), like the non-reducing C1'-H of maltose
    "me-a-glucoside": _A.replace("(O1)O)", "(O1)OC)"),
    # 4-O-methyl beta-glucose: C4-H = H-C-O-C, like the reducing-unit C4-H of maltose
    "4-OMe-b-glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)OC)O",
}


def smiles_of(name):
    return SMILES.get(name) or MODELS[name]


# ---------------------------------------------------------------- xTB engine
class XTB:
    def __init__(self, numbers, solvent="water"):
        self.numbers = np.asarray(numbers)
        self.solvent = solvent

    def _new(self, pos_A):
        from tblite.interface import Calculator
        c = Calculator("GFN2-xTB", self.numbers, np.asarray(pos_A) / BOHR)
        c.set("verbosity", 0)
        c.set("max-iter", 1000)
        if self.solvent:
            c.add("alpb-solvation", self.solvent)
        return c

    def __call__(self, pos_A, grad=False, restart=None):
        c = self._new(pos_A)
        r = c.singlepoint(restart) if restart is not None else c.singlepoint()
        self.last = r
        out = {"energy": r.get("energy"), "dipole": r.get("dipole")}
        if grad:
            out["gradient"] = r.get("gradient") / BOHR       # Eh/A
        return out


def optimise(engine, pos, gtol=3e-4, maxiter=600):
    shape = pos.shape

    def f(x):
        r = engine(x.reshape(shape), grad=True)
        return r["energy"], r["gradient"].ravel()

    res = minimize(f, pos.ravel(), jac=True, method="L-BFGS-B",
                   options={"gtol": gtol, "maxiter": maxiter})
    return res.x.reshape(shape), res.fun


# ------------------------------------------------------------ structure info
def label_atoms(mol):
    """Return {H index: dict(label, carbon, unit, cls)} for every C-H."""
    ri = mol.GetRingInfo()
    labels = {}
    units = []
    for ring in ri.AtomRings():
        if len(ring) != 6:
            continue
        ring = list(ring)
        o5 = [a for a in ring if mol.GetAtomWithIdx(a).GetSymbol() == "O"][0]
        cs = [a for a in ring if a != o5]
        nb_o5 = [n.GetIdx() for n in mol.GetAtomWithIdx(o5).GetNeighbors()]
        c1 = [c for c in nb_o5 if sum(n.GetSymbol() == "O" for n in
                                      mol.GetAtomWithIdx(c).GetNeighbors()) == 2][0]
        c5 = [c for c in nb_o5 if c != c1][0]
        order = [c1]
        while len(order) < 5:
            prev = order[-1]
            nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(prev).GetNeighbors()
                   if n.GetIdx() in cs and n.GetIdx() not in order]
            order.append(nxt[0])
        c6 = [n.GetIdx() for n in mol.GetAtomWithIdx(c5).GetNeighbors()
              if n.GetSymbol() == "C" and n.GetIdx() not in ring][0]
        units.append({"ring": ring, "o5": o5, "C": order + [c6]})
    # reducing unit = unit whose C1 carries an OH
    def has_oh(c, exclude):
        for n in mol.GetAtomWithIdx(c).GetNeighbors():
            if n.GetSymbol() == "O" and n.GetIdx() not in exclude:
                return any(m.GetSymbol() == "H" for m in n.GetNeighbors())
        return False
    for u in units:
        u["reducing"] = has_oh(u["C"][0], [u["o5"]])
    if len(units) == 1:
        units[0]["name"] = ""
    else:
        for u in units:
            u["name"] = "red." if u["reducing"] else "non-red."
    for u in units:
        for pos, c in enumerate(u["C"], start=1):
            atom = mol.GetAtomWithIdx(c)
            exo_o = [n for n in atom.GetNeighbors() if n.GetSymbol() == "O" and n.GetIdx() != u["o5"]]
            o_types = []
            for o in exo_o:
                o_types.append("OH" if any(m.GetSymbol() == "H" for m in o.GetNeighbors()) else "O-C(glyc)")
            ring_o = any(n.GetIdx() == u["o5"] for n in atom.GetNeighbors())
            if pos == 1:
                cls = "anomeric H-C(O5)(" + o_types[0] + ")"
            elif pos == 6:
                cls = "CH2-OH"
            elif ring_o:
                cls = "H-C-O-C (ring O5)"
            else:
                cls = "H-C-" + o_types[0].replace("OH", "O-H")
            hs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "H"]
            for k, h in enumerate(hs):
                lab = f"C{pos}-H" + ("" if len(hs) == 1 else "ab"[k])
                if u["name"]:
                    lab = f"{u['name']} {lab}"
                labels[h] = dict(label=lab, carbon=c, pos=pos, unit=u["name"], cls=cls,
                                 ring=u["ring"], o5=u["o5"], ringC=u["C"][:5])
    return labels


def ring_orientation(pos, ring, atom, sub):
    """'ax' or 'eq' for substituent `sub` on ring atom `atom`."""
    P = pos[ring]
    c = P.mean(0)
    _, _, vt = np.linalg.svd(P - c)
    n = vt[2]
    v = pos[sub] - pos[atom]
    v /= np.linalg.norm(v)
    return "ax" if abs(v @ n) > 0.75 else "eq"


def is_chair(pos, order):
    """order = [O5, C1..C5]; all ring torsions 40-75 deg with alternating sign."""
    def dih(a, b, c, d):
        b0, b1, b2 = pos[a] - pos[b], pos[c] - pos[b], pos[d] - pos[c]
        b1 = b1 / np.linalg.norm(b1)
        v = b0 - (b0 @ b1) * b1
        w = b2 - (b2 @ b1) * b1
        return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))
    t = [dih(*[order[(i + k) % 6] for k in range(4)]) for i in range(6)]
    return all(40 < abs(x) < 75 for x in t) and all(np.sign(t[i]) != np.sign(t[(i + 1) % 6])
                                                      for i in range(6))


def all_4c1(mol, pos, labels):
    units = {d["unit"]: [d["o5"]] + d["ringC"] for d in labels.values()}
    if not all(is_chair(pos, o) for o in units.values()):
        return False
    return not any(o == "ax" and p != 1 for _, p, o in check_stereo(mol, pos, labels))


def check_stereo(mol, pos, labels):
    """Orientation of exocyclic heavy substituents on each ring carbon."""
    out = []
    seen = set()
    for h, d in labels.items():
        if d["pos"] > 5 or (d["carbon"] in seen):
            continue
        seen.add(d["carbon"])
        atom = mol.GetAtomWithIdx(d["carbon"])
        sub = [n.GetIdx() for n in atom.GetNeighbors()
               if n.GetSymbol() != "H" and n.GetIdx() not in d["ring"]]
        out.append((d["unit"], d["pos"], ring_orientation(pos, d["ring"], d["carbon"], sub[0])))
    return sorted(out)


# --------------------------------------------------------------- conformers
def conformers(name, n_embed=60, keep=5, window_kcal=2.0, seed=7, log=print):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles_of(name)))
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, randomSeed=seed))
    mm = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    e_mm = np.array([e for _, e in mm])
    order = np.argsort(e_mm)
    labels = label_atoms(mol)
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    eng = XTB(numbers)
    out = []
    cand = [ci for ci in order
            if all_4c1(mol, mol.GetConformer(int(cids[ci])).GetPositions(), labels)]
    log(f"  {name}: {len(cids)} embedded, {len(cand)} with 4C1 chairs")
    for ci in cand[: max(3 * keep, 12)]:
        pos = mol.GetConformer(int(cids[ci])).GetPositions()
        # require 4C1 chairs: only the anomeric C1 substituent may be axial
        if not all_4c1(mol, pos, labels):
            continue
        p_opt, e = optimise(eng, pos)
        if not all_4c1(mol, p_opt, labels):
            continue
        st = check_stereo(mol, p_opt, labels)
        if any(abs(e - e2) * 627.5 < 0.05 for _, e2, _ in out):
            continue
        out.append((p_opt, e, st))
        log(f"  {name}: conformer E={e:.6f} Eh  stereo={st}")
    out.sort(key=lambda t: t[1])
    e0 = out[0][1]
    out = [t for t in out if (t[1] - e0) * 627.5 <= window_kcal][:keep]
    w = np.exp(-np.array([(t[1] - e0) * 627.5 for t in out]) / 0.593)
    return mol, labels, numbers, [t[0] for t in out], w / w.sum(), out[0][2]


# ------------------------------------------------------------- local modes
Q = np.linspace(-0.30, 0.70, 15)       # A


def scan_ch(engine, pos, h, c, Q=Q):
    """Energies/dipoles along the C-H stretch; failed points are NaN.
    Scans outward from q=0 so each point restarts from its neighbour."""
    u = pos[h] - pos[c]
    u /= np.linalg.norm(u)
    E = np.full(Q.size, np.nan)
    D = np.full((Q.size, 3), np.nan)
    i0 = int(np.argmin(np.abs(Q)))
    for idx in (list(range(i0, Q.size)), list(range(i0, -1, -1))):
        restart = None
        for i in idx:
            p = pos.copy()
            p[h] = pos[h] + Q[i] * u
            try:
                r = engine(p, restart=restart) if isinstance(engine, XTB) else engine(p)
            except Exception:
                break                       # stop going further out on this side
            E[i], D[i] = r["energy"], r["dipole"]
            restart = getattr(engine, "last", None)
    return E, D


def solve_local_mode(E, D, n_levels=5, mu=M_C * M_H / (M_C + M_H), n_grid=1200, Q=Q):
    ok = np.isfinite(E)
    q, E, D = Q[ok], E[ok], D[ok]
    V = (E - E.min()) * HARTREE_CM
    qg = np.linspace(q[0], q[-1], n_grid)
    h = qg[1] - qg[0]
    Vg = CubicSpline(q, V)(qg)
    B = B_CONST / mu
    ev, psi = eigh_tridiagonal(2 * B / h ** 2 + Vg, -B / h ** 2 * np.ones(n_grid - 1),
                               select="i", select_range=(0, n_levels - 1))
    psi /= np.sqrt(h)
    Dg = np.column_stack([CubicSpline(q, D[:, k])(qg) for k in range(3)])
    nu = ev[1:] - ev[0]
    M = np.array([np.linalg.norm(((psi[:, 0] * psi[:, v])[:, None] * Dg).sum(0) * h)
                  for v in range(1, n_levels)])
    inten = nu * M ** 2
    # local-mode Morse parameters from the first two levels
    wexe = (2 * nu[0] - nu[1]) / 2
    we = nu[0] + 2 * wexe
    return nu, inten, we, wexe
