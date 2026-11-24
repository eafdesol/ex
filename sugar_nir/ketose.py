"""Fructose / allulose (D-psicose) model with tautomer equilibria.

Both are ketohexoses C6H12O6 (allulose = C3 epimer of fructose). In water each
exists as an equilibrium of alpha/beta-pyranose, alpha/beta-furanose and a
trace of the open keto form. Because the formula is identical, the NIR
difference between the two sugars comes almost entirely from

  * different tautomer populations (fructose ~68 % beta-pyranose,
    allulose ~39 % alpha-furanose), and
  * different OH orientation / H-bonding (allulose C3-OH axial in pyranose,
    cis-3,4-diol in furanose -> intramolecular H-bond).

Every tautomer is described with its own X-H oscillators; the effective
spectrum of a sugar is the population-weighted sum, and populations drift with
temperature (furanose share increases on heating).

Spectral axis: wavenumber 4000-16000 cm^-1 (625-2500 nm), 8 cm^-1 step.

Populations at 25 C (D2O NMR, literature, rounded):
  fructose: b-pyr 68.4, b-fur 22.4, a-fur 6.2, a-pyr 2.7, keto 0.5 %
  allulose: a-fur 39, a-pyr 22, b-pyr 24, b-fur 15, keto 0.3 %
Conformations / OH orientations are simplified (see TAUTOMERS).
"""
from __future__ import annotations

import numpy as np

from .bands import Group, group_bands, render
from .simulate import noise_sigma, pigment_spectrum, water_spectrum

K = np.arange(4000.0, 16000.0 + 1e-9, 8.0)     # cm^-1, instrument range
WLK = 1e7 / K                                    # nm (descending)
MW = 180.16
SCALE = 0.9
PSV = 0.0062

KGROUPS = {
    # hydroxyls (H-bonded, aqueous)
    "OH_eq":    Group("OH_eq",    3560, 90, 1.0, 120, [(1380, 1.6, 110, "OH str + COH bend"),
                                                       (1050, 0.6, 90, "OH str + CO str")]),
    "OH_ax":    Group("OH_ax",    3595, 88, 1.0, 110, [(1395, 1.6, 105, "OH(ax) str + COH bend"),
                                                       (1060, 0.6, 90, "OH(ax) str + CO str")]),
    "OH_hb":    Group("OH_hb",    3530, 92, 1.1, 130, [(1410, 1.6, 115, "intramol. H-bonded OH str + bend"),
                                                       (1065, 0.6, 95, "H-bonded OH str + CO str")]),
    "OH_prim":  Group("OH_prim",  3575, 90, 1.0, 120, [(1360, 1.6, 110, "CH2OH OH str + bend"),
                                                       (1040, 0.6, 90, "CH2OH OH str + CO str")]),
    "OH_anomK": Group("OH_anomK", 3550, 91, 1.0, 125, [(1400, 1.6, 110, "hemiketal OH str + bend"),
                                                       (1085, 0.7, 90, "hemiketal OH + C2-O str")]),
    # C-H
    "CH_ax":    Group("CH_ax",    2980, 58, 0.55, 50, [(1440, 1.3, 45, "axial CH str + bend"),
                                                       (1110, 0.5, 50, "axial CH + CO str")]),
    "CH_eq":    Group("CH_eq",    3010, 58, 0.55, 50, [(1430, 1.3, 45, "equatorial CH str + bend"),
                                                       (1100, 0.5, 50, "equatorial CH + CO str")]),
    "CH_fur":   Group("CH_fur",   3025, 59, 0.55, 55, [(1420, 1.3, 45, "furanose CH str + bend"),
                                                       (1080, 0.5, 50, "furanose CH + CO str")]),
    "CH2_ring": Group("CH2_ring", 3000, 60, 0.55, 48, [(1470, 1.3, 40, "ring O-CH2 str + scissor"),
                                                       (1350, 0.7, 45, "ring O-CH2 str + wag")]),
    "CH2":      Group("CH2",      2990, 60, 0.55, 50, [(1460, 1.3, 40, "CH2OH str + scissor"),
                                                       (1335, 0.7, 45, "CH2OH str + wag")]),
}

# tautomer -> list of (group, count, d_we).  d_we encodes alpha/beta and
# epimer-specific environment shifts (cm^-1).
TAUTOMERS = {
    "fructose": {
        "b-pyr": [("OH_anomK", 1, +10), ("OH_prim", 1, 0), ("OH_eq", 2, 0), ("OH_ax", 1, 0),
                  ("CH_ax", 2, 0), ("CH_eq", 1, 0), ("CH2_ring", 2, 0), ("CH2", 2, 0)],
        "a-pyr": [("OH_anomK", 1, -10), ("OH_prim", 1, 0), ("OH_eq", 2, 0), ("OH_ax", 1, 0),
                  ("CH_ax", 2, 0), ("CH_eq", 1, 0), ("CH2_ring", 2, 0), ("CH2", 2, +5)],
        "b-fur": [("OH_anomK", 1, +5), ("OH_prim", 2, 0), ("OH_eq", 2, -5),
                  ("CH_fur", 3, 0), ("CH2", 4, 0)],
        "a-fur": [("OH_anomK", 1, -5), ("OH_prim", 2, 0), ("OH_eq", 2, -5),
                  ("CH_fur", 3, -5), ("CH2", 4, 0)],
        "keto":  [("OH_prim", 2, 0), ("OH_eq", 3, 0), ("CH_eq", 3, 0), ("CH2", 4, 0)],
    },
    "allulose": {   # C3 epimer: C3-OH axial in pyranose; cis-3,4-diol in furanose
        "b-pyr": [("OH_anomK", 1, +10), ("OH_prim", 1, 0), ("OH_eq", 1, 0), ("OH_hb", 2, 0),
                  ("CH_ax", 1, 0), ("CH_eq", 2, 0), ("CH2_ring", 2, 0), ("CH2", 2, 0)],
        "a-pyr": [("OH_anomK", 1, -10), ("OH_prim", 1, 0), ("OH_eq", 1, 0), ("OH_hb", 2, 0),
                  ("CH_ax", 1, 0), ("CH_eq", 2, 0), ("CH2_ring", 2, 0), ("CH2", 2, +5)],
        "b-fur": [("OH_anomK", 1, +5), ("OH_prim", 2, 0), ("OH_hb", 2, 0),
                  ("CH_fur", 3, +5), ("CH2", 4, 0)],
        "a-fur": [("OH_anomK", 1, -5), ("OH_prim", 2, 0), ("OH_hb", 2, +5),
                  ("CH_fur", 3, 0), ("CH2", 4, 0)],
        "keto":  [("OH_prim", 2, 0), ("OH_eq", 2, 0), ("OH_ax", 1, 0), ("CH_eq", 3, 0), ("CH2", 4, 0)],
    },
}
SUGARS_K = list(TAUTOMERS)
TAUT_NAMES = ["b-pyr", "a-pyr", "b-fur", "a-fur", "keto"]

# population (fraction) at 25 C and its temperature slope (per C)
POP25 = {
    "fructose": {"b-pyr": .684, "a-pyr": .027, "b-fur": .224, "a-fur": .062, "keto": .005},
    "allulose": {"b-pyr": .240, "a-pyr": .220, "b-fur": .150, "a-fur": .387, "keto": .003},
}
DPOP = {
    "fructose": {"b-pyr": -.0035, "a-pyr": 0.0, "b-fur": +.0025, "a-fur": +.0008, "keto": +.0002},
    "allulose": {"b-pyr": -.0010, "a-pyr": -.0010, "b-fur": +.0005, "a-fur": +.0014, "keto": +.0001},
}


def populations(sugar, T=25.0, pop=None):
    pop = pop or POP25[sugar]
    p = np.array([max(pop[t] + DPOP[sugar][t] * (T - 25.0), 0.0) for t in TAUT_NAMES])
    return dict(zip(TAUT_NAMES, p / p.sum()))


def tautomer_bands(sugar, taut, perturb=None):
    out = []
    for gname, count, d_we in TAUTOMERS[sugar][taut]:
        dw, ws, isc = (0.0, 1.0, 1.0)
        if perturb and gname in perturb:
            dw, ws, isc = perturb[gname]
        for c, I, sig, lab in group_bands(KGROUPS[gname], d_we + dw, ws):
            out.append((c, count * I * isc, sig, f"{sugar} {taut}: {lab}"))
    if taut == "keto":   # C=O 2nd overtone (~3450, out of range) + C=O/CH combination ~4700
        out.append((1735 + 2990 - 20, 0.02, 60, f"{sugar} keto: C=O str + CH str"))
    return out


def tautomer_spectrum(sugar, taut, k=K, perturb=None):
    """Absorbance per 1 g/100 mL of this sugar if it were 100 % this tautomer."""
    return SCALE * render(tautomer_bands(sugar, taut, perturb), 1e7 / k) / MW


def sugar_spectrum(sugar, k=K, T=25.0, perturb=None, pop=None):
    p = populations(sugar, T, pop)
    return sum(p[t] * tautomer_spectrum(sugar, t, k, perturb) for t in TAUT_NAMES)


def band_table_k(kmin=4000, kmax=16000, min_rel=0.02):
    rows = []
    for s in SUGARS_K:
        p = populations(s)
        for t in TAUT_NAMES:
            for c, I, sig, lab in tautomer_bands(s, t):
                rows.append((s, t, c, I * p[t], lab))
    imax = max(r[3] for r in rows)
    rows = [r for r in rows if kmin <= r[2] <= kmax and r[3] / imax >= min_rel]
    return sorted(rows, key=lambda r: r[2])


def simulate_ketose(n, rng, k=K, total=(5.0, 40.0), T_mean=25.0, T_sd=2.0):
    """Aqueous fructose/allulose mixtures.

    Design: total sugar ~ U(5, 40) g/100 mL, allulose fraction ~ U(0, 1)
    (as in fructose -> allulose epimerisation process samples).
    Y columns: fructose, allulose (g/100 mL).  Also returns temperatures.
    """
    tot = rng.uniform(*total, n)
    frac = rng.uniform(0, 1, n)
    Y = np.column_stack([tot * (1 - frac), tot * frac])
    wl = 1e7 / k
    X = np.empty((n, k.size))
    Ts = rng.normal(T_mean, T_sd, n)
    for i in range(n):
        T = Ts[i]
        A = water_spectrum(wl, T) * (1 - PSV * tot[i])
        for j, s in enumerate(SUGARS_K):
            A += Y[i, j] * sugar_spectrum(s, k, T)
        A += rng.uniform(0, 0.4) * pigment_spectrum(wl, rng.uniform(430, 520), rng.uniform(40, 80))
        A *= rng.normal(1.0, 0.02)
        x = k / 10000.0
        A += rng.normal(0, 0.02) + rng.normal(0, 0.01) * x + rng.uniform(0, 0.03) * x ** 2
        A += rng.normal(0, 1, k.size) * noise_sigma(A)
        X[i] = A
    return X, Y, Ts


def theory_selectivity_k(rng, k=K, model_error=True, dominant_only=False):
    """NAS / noise per wavenumber for fructose, allulose, total, allulose fraction.

    Interferents: the other sugar, water (+ temperature derivative), the
    temperature-driven tautomer shift of each sugar, pigment, baselines.
    dominant_only=True emulates a naive model that ignores the tautomer
    equilibrium (each sugar = its single most abundant tautomer).
    """
    from .pls import preprocess
    perturb = None
    if model_error:
        perturb = {g: (rng.normal(0, 15), rng.uniform(0.85, 1.15), rng.uniform(0.8, 1.2))
                   for g in KGROUPS}
    wl = 1e7 / k
    if dominant_only:
        dom = {"fructose": "b-pyr", "allulose": "a-fur"}
        S = {s: tautomer_spectrum(s, dom[s], k, perturb) for s in SUGARS_K}
        dS = {s: np.zeros_like(k) for s in SUGARS_K}
    else:
        S = {s: sugar_spectrum(s, k, 25.0, perturb) for s in SUGARS_K}
        dS = {s: sugar_spectrum(s, k, 30.0, perturb) - S[s] for s in SUGARS_K}
    w0 = water_spectrum(wl, 25.0)
    x = k / 10000.0
    interf = [w0, water_spectrum(wl, 26.0) - w0, pigment_spectrum(wl), pigment_spectrum(wl, 500, 70),
              np.ones_like(k), x, x ** 2] + [v for v in dS.values() if np.any(v)]
    typical = w0 + 10 * S["fructose"] + 10 * S["allulose"]
    P = lambda v: preprocess(np.atleast_2d(v))[0]
    noise = np.abs(P(noise_sigma(typical))) + 0.1 * noise_sigma(typical)
    D = {s: P(v) for s, v in S.items()}
    Di = [P(v) for v in interf]

    def nas(v, others):
        Z = np.array(others).T
        return np.abs(v - Z @ np.linalg.pinv(Z) @ v) / noise

    return {
        "fructose": nas(D["fructose"], [D["allulose"]] + Di),
        "allulose": nas(D["allulose"], [D["fructose"]] + Di),
        "total": nas(0.5 * (D["fructose"] + D["allulose"]), Di),
        "allulose_pct": nas(D["allulose"] - D["fructose"],
                            Di + [0.5 * (D["fructose"] + D["allulose"])]),
    }
