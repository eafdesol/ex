"""General X-H local-mode overtone calculation for any molecule.

Extends the glucose/maltose C-H workflow to every X-H bond (C-H, O-H, N-H,
S-H) of any molecule given as SMILES, including several forms (anomers,
tautomers) with populations.

A form is built either from a verified isomeric SMILES ("smiles") whose atoms
get labels from a mapped template, or from a mapped template plus target CIP
descriptors ("cip"), in which case the matching stereoisomer is picked
automatically. Atom-map numbers on carbons give readable labels:
1..9 -> C1..C9, 11..19 -> C1'..C9', 21..29 -> C1''..C9''.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdCIPLabeler
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

from .localmode import XTB, Q, is_chair, optimise, ring_orientation, scan_ch, solve_local_mode

MASS = {"H": 1.007825, "C": 12.0, "N": 14.003074, "O": 15.994915, "S": 31.972071}


# ------------------------------------------------------------ building forms
def _strip_maps(m):
    m = Chem.Mol(m)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return m


def form_from_cip(template, cip):
    """Stereoisomer of `template` whose mapped atoms have the given CIP codes."""
    base = Chem.MolFromSmiles(template)
    hits = []
    for iso in EnumerateStereoisomers(base):
        rdCIPLabeler.AssignCIPLabels(iso)
        ok = True
        for a in iso.GetAtoms():
            want = cip.get(a.GetAtomMapNum())
            if want and (not a.HasProp("_CIPCode") or a.GetProp("_CIPCode") != want):
                ok = False
                break
        if ok:
            hits.append(iso)
    if len(hits) != 1:
        raise ValueError(f"{len(hits)} stereoisomers match {cip} for {template}")
    return hits[0]


def form_from_smiles(smiles, template):
    """Verified isomeric SMILES; atom maps copied from the template."""
    m = Chem.MolFromSmiles(smiles)
    t = Chem.MolFromSmiles(template)
    match = m.GetSubstructMatch(_strip_maps(t))
    if not match:
        raise ValueError("template does not match smiles")
    for ti, mi in enumerate(match):
        m.GetAtomWithIdx(mi).SetAtomMapNum(t.GetAtomWithIdx(ti).GetAtomMapNum())
    return m


def build_form(spec):
    if "smiles" in spec:
        return form_from_smiles(spec["smiles"], spec["template"])
    return form_from_cip(spec["template"], spec["cip"])


def iso_smiles(m):
    return Chem.MolToSmiles(_strip_maps(m))


# ------------------------------------------------------------ labels
def carbon_name(mapnum):
    if not mapnum:
        return None
    unit, pos = divmod(mapnum, 10)
    return f"C{pos}" + "'" * unit


def xh_bonds(mol):
    """All X-H bonds: list of dict(h, x, elem, label, kind, ...)."""
    out = []
    ri = mol.GetRingInfo()
    for x in mol.GetAtoms():
        if x.GetSymbol() == "H":
            continue
        hs = [n.GetIdx() for n in x.GetNeighbors() if n.GetSymbol() == "H"]
        if not hs:
            continue
        el = x.GetSymbol()
        heavy = [n for n in x.GetNeighbors() if n.GetSymbol() != "H"]
        if el == "C":
            base = carbon_name(x.GetAtomMapNum()) or f"C#{x.GetIdx()}"
            n_o = sum(n.GetSymbol() == "O" for n in heavy)
            kind = {1: "CH", 2: "CH2", 3: "CH3"}.get(len(hs), "CH4")
            ring = ri.NumAtomRings(x.GetIdx()) > 0
            if n_o >= 2:
                cls = "아노머 C-H (O 2개)"
            elif kind == "CH2" and n_o == 1:
                cls = "CH2-O (고리)" if ring else "CH2OH (고리 밖)"
            elif n_o == 1:
                cls = "고리 C-H (O 1개)" if ring else "사슬 C-H (O 1개)"
            else:
                cls = f"{kind} (O 없음)"
            for k, h in enumerate(hs):
                lab = f"{base}-H" + ("" if len(hs) == 1 else "ab"[k] if k < 2 else str(k))
                out.append(dict(h=h, x=x.GetIdx(), elem="C", label=lab, kind=kind, cls=cls))
        elif el == "O":
            c = heavy[0] if heavy else None
            cname = carbon_name(c.GetAtomMapNum()) if c is not None else None
            lab = (f"O{cname[1:]}-H" if cname else f"O#{x.GetIdx()}-H")
            if c is not None and c.GetSymbol() == "C":
                c_o = sum(n.GetSymbol() == "O" for n in c.GetNeighbors())
                c_h = sum(n.GetSymbol() == "H" for n in c.GetNeighbors())
                if c_o >= 2:
                    cls = "아노머 O-H (헤미아세탈·헤미케탈)"
                elif c_h == 2:
                    cls = "1차 O-H (CH2OH)"
                elif c_h == 1:
                    cls = "2차 O-H"
                else:
                    cls = "3차 O-H"
            else:
                cls = "기타 O-H"
            out.append(dict(h=hs[0], x=x.GetIdx(), elem="O", label=lab, kind="OH", cls=cls))
        else:
            for k, h in enumerate(hs):
                out.append(dict(h=h, x=x.GetIdx(), elem=el, label=f"{el}#{x.GetIdx()}-H{k}",
                                kind=f"{el}H", cls=f"{el}-H"))
    return out


# ------------------------------------------------------------ conformers
def conformers(mol, n_embed=80, keep=5, window_kcal=2.0, seed=7, log=print, name=""):
    mol = Chem.AddHs(mol)
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, randomSeed=seed))
    mm = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    order = np.argsort([e for _, e in mm])
    rings6 = [list(r) for r in mol.GetRingInfo().AtomRings() if len(r) == 6]
    chairs = lambda p: all(is_chair(p, r) for r in rings6)
    cand = [i for i in order if chairs(mol.GetConformer(int(cids[i])).GetPositions())]
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    eng = XTB(numbers)
    found = []
    for i in cand[: max(3 * keep, 12)]:
        p, e = optimise(eng, mol.GetConformer(int(cids[i])).GetPositions())
        if not chairs(p) or any(abs(e - e2) * 627.5 < 0.05 for _, e2 in found):
            continue
        found.append((p, e))
    found.sort(key=lambda t: t[1])
    e0 = found[0][1]
    found = [t for t in found if (t[1] - e0) * 627.5 <= window_kcal][:keep]
    w = np.exp(-np.array([(e - e0) * 627.5 for _, e in found]) / 0.593)
    log(f"  {name}: {len(cids)} embedded, {len(cand)} usable, {len(found)} kept")
    return mol, numbers, [p for p, _ in found], w / w.sum()


# ------------------------------------------------------------ descriptors
def _dih(p, a, b, c, d):
    b0, b1, b2 = p[a] - p[b], p[c] - p[b], p[d] - p[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def descriptors(mol, bond, confs, w):
    h, x = bond["h"], bond["x"]
    r = float(np.average([np.linalg.norm(p[h] - p[x]) for p in confs], weights=w))
    out = {"r": round(r, 4)}
    xa = mol.GetAtomWithIdx(x)
    if bond["elem"] == "C":
        lp = 0.0
        for o in [n for n in xa.GetNeighbors() if n.GetSymbol() == "O"]:
            ys = [n.GetIdx() for n in o.GetNeighbors() if n.GetIdx() != x]
            if ys:
                ds = [_dih(p, h, x, o.GetIdx(), ys[0]) for p in confs]
                lp += float(np.average([1.0 if 25 < abs(d) < 95 else 0.0 for d in ds], weights=w))
        out["lp"] = round(lp, 2)
        ring6 = [list(rg) for rg in mol.GetRingInfo().AtomRings() if len(rg) == 6 and x in rg]
        out["orient"] = (ring_orientation(confs[0], ring6[0], x, h) if ring6 and bond["kind"] == "CH"
                         else "")
    else:
        others = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("O", "N")
                  and a.GetIdx() != x]
        dmin = [min(np.linalg.norm(p[h] - p[o]) for o in others) if others else 9.9 for p in confs]
        out["hb"] = round(float(np.average([d < 2.5 for d in dmin], weights=w)), 2)
        out["hb_d"] = round(float(np.average(dmin, weights=w)), 2)
    return out


# ------------------------------------------------------------ local modes
def reduced_mass(elem):
    m = MASS[elem]
    return m * MASS["H"] / (m + MASS["H"])


def compute_form(mol0, n_embed=80, keep=5, log=print, name=""):
    """Full X-H local-mode result for one form (all conformers, weighted)."""
    mol, numbers, confs, w = conformers(mol0, n_embed=n_embed, keep=keep, log=log, name=name)
    bonds = xh_bonds(mol)
    eng = XTB(numbers)
    rows = []
    for b in bonds:
        res = []
        for p in confs:
            E, D = scan_ch(eng, p, b["h"], b["x"])
            nu, I, we, wexe = solve_local_mode(E, D, mu=reduced_mass(b["elem"]))
            res.append((nu, I, we, wexe))
        nu = np.average([r[0] for r in res], axis=0, weights=w)
        I = np.average([r[1] for r in res], axis=0, weights=w)
        row = {k: b[k] for k in ("h", "x", "elem", "label", "kind", "cls")}
        row.update(descriptors(mol, b, confs, w))
        row.update({"nu1": round(float(nu[0]), 1), "nu2": round(float(nu[1]), 1),
                    "nu3": round(float(nu[2]), 1), "I1": float(I[0]), "I2": float(I[1]),
                    "I3": float(I[2]),
                    "we": round(float(np.average([r[2] for r in res], weights=w)), 1),
                    "wexe": round(float(np.average([r[3] for r in res], weights=w)), 2),
                    "nu2_sd": round(float(np.std([r[0][1] for r in res])), 1)})
        rows.append(row)
        log(f"    {name} {b['label']:10s} {row['nu2']:7.1f} cm-1")
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = f"{len(sym)}\n{name}\n" + "\n".join(
        f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, (x, y, z) in zip(sym, confs[0]))
    return {"bonds": rows, "xyz": xyz, "weights": [float(v) for v in w],
            "smiles": iso_smiles(mol0), "mw": Descriptors.MolWt(_strip_maps(mol0))}
