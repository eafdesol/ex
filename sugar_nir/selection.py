"""Wavelength-region selection strategies.

theory_nas : physics prior. Net-analyte-signal of each sugar against all other
             sugars, water (+ temperature derivative), pigment and baseline
             shapes, divided by expected noise. Uses a *perturbed* band model
             to mimic theory-vs-reality error, and needs no measured data.
vip        : VIP > 1 from a full-spectrum PLS model.
ipls       : forward interval PLS on equal-width intervals (data-driven).
hybrid     : forward iPLS restricted to the theory-ranked candidate windows.
"""
from __future__ import annotations

import numpy as np

from .bands import GROUPS, SUGAR_NAMES
from .pls import PLS1, choose_lv, cv_rmse, preprocess
from .simulate import (WL, noise_sigma, pigment_spectrum, pure_sugar_spectrum,
                       water_spectrum)


def mask_to_windows(mask, wl=WL, min_width=10.0, merge_gap=12.0):
    wins, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if (not m or i == len(mask) - 1) and start is not None:
            end = i if m else i - 1
            wins.append([wl[start], wl[end]])
            start = None
    merged = []
    for w in wins:
        if merged and w[0] - merged[-1][1] <= merge_gap:
            merged[-1][1] = w[1]
        else:
            merged.append(w)
    return [tuple(w) for w in merged if w[1] - w[0] >= min_width]


def windows_to_mask(wins, wl=WL):
    m = np.zeros(wl.size, bool)
    for a, b in wins:
        m |= (wl >= a) & (wl <= b)
    return m


def theory_selectivity(rng, wl=WL, model_error=True):
    """Per-wavelength NAS/noise for each sugar (dict name -> array)."""
    perturb = None
    if model_error:  # +-15 cm^-1 position, +-15 % width, +-20 % intensity error
        perturb = {g: (rng.normal(0, 15), rng.uniform(0.85, 1.15),
                       rng.uniform(0.8, 1.2)) for g in GROUPS}
    S = {n: pure_sugar_spectrum(n, wl, perturb) for n in SUGAR_NAMES}
    w0, w1 = water_spectrum(wl, 25.0), water_spectrum(wl, 26.0)
    x = wl / 1000.0
    interf = [w0, w1 - w0, pigment_spectrum(wl), pigment_spectrum(wl, 500, 70),
              np.ones_like(wl), x, 1 / x ** 2]
    typical = w0 + sum(5.0 * s for s in S.values())
    noise = preprocess(noise_sigma(typical)[None] * np.ones((1, wl.size)))
    noise = np.abs(noise[0]) + noise_sigma(typical) * 0.1   # derivative noise
    D = {k: preprocess(v[None])[0] for k, v in S.items()}
    Dint = [preprocess(v[None])[0] for v in interf]
    out = {}
    for name in SUGAR_NAMES:
        Z = np.array([D[o] for o in SUGAR_NAMES if o != name] + Dint).T
        Pz = Z @ np.linalg.pinv(Z)
        nas = D[name] - Pz @ D[name]
        out[name] = np.abs(nas) / noise
    # maltotriose + maltotetraose together (DP3+DP4)
    Z = np.array([D["glucose"], D["maltose"]] + Dint).T
    oligo = 0.5 * (D["maltotriose"] + D["maltotetraose"])
    out["DP3_4_sum"] = np.abs(oligo - Z @ np.linalg.pinv(Z) @ oligo) / noise
    # aggregate targets: total sugar (mean per-gram spectrum) and reducing ends
    # (direction that varies with 1/MW, i.e. moles of reducing end per gram)
    from .bands import SUGARS
    Dm = np.array([D[n] for n in SUGAR_NAMES])
    Zi = np.array(Dint).T
    tot = Dm.mean(0)
    out["total"] = np.abs(tot - Zi @ np.linalg.pinv(Zi) @ tot) / noise
    inv = np.array([1 / SUGARS[n].mw for n in SUGAR_NAMES])
    red = ((inv - inv.mean()) / inv.std()) @ Dm
    Z = np.column_stack([Zi, tot])
    out["reducing_ends"] = np.abs(red - Z @ np.linalg.pinv(Z) @ red) / noise
    out["mean_DP"] = out["reducing_ends"]
    return out


def theory_windows(sel, top_frac=0.15, smooth=7, wl=WL):
    k = np.ones(smooth) / smooth
    s = np.convolve(sel, k, mode="same")
    thr = np.quantile(s, 1 - top_frac)
    return mask_to_windows(s >= thr, wl)


def vip_windows(Xp, y, wl=WL, n_max=15):
    lv = choose_lv(cv_rmse(Xp, y, n_max))
    v = PLS1(n_max).fit(Xp, y).vip(lv)
    v = np.convolve(v, np.ones(7) / 7, mode="same")
    return mask_to_windows(v > 1.0, wl)


def forward_ipls(Xp, y, intervals, n_max=12, max_steps=8):
    """Greedy forward selection of intervals (list of index arrays)."""
    chosen, best = [], np.inf
    for _ in range(max_steps):
        cand = None
        for j, iv in enumerate(intervals):
            if j in chosen:
                continue
            cols = np.concatenate([intervals[c] for c in chosen] + [iv])
            r = cv_rmse(Xp[:, cols], y, n_max).min()
            if cand is None or r < cand[1]:
                cand = (j, r)
        if cand is None or cand[1] >= best * 0.99:
            break
        chosen.append(cand[0]); best = cand[1]
    return chosen, best


def ipls_windows(Xp, y, wl=WL, width_nm=50.0):
    edges = np.arange(wl[0], wl[-1] + width_nm, width_nm)
    ivs = [np.where((wl >= a) & (wl < b))[0] for a, b in zip(edges[:-1], edges[1:])]
    ivs = [i for i in ivs if i.size]
    chosen, _ = forward_ipls(Xp, y, ivs)
    m = np.zeros(wl.size, bool)
    for c in chosen:
        m[ivs[c]] = True
    return mask_to_windows(m, wl, min_width=0)


def backward_ipls(Xp, y, intervals, n_max=12, tol=1.0):
    """Start from all intervals; drop one at a time while RMSECV does not worsen."""
    keep = list(range(len(intervals)))
    cur = cv_rmse(Xp[:, np.concatenate([intervals[k] for k in keep])], y, n_max).min()
    while len(keep) > 1:
        trials = []
        for k in keep:
            rest = [intervals[j] for j in keep if j != k]
            trials.append((cv_rmse(Xp[:, np.concatenate(rest)], y, n_max).min(), k))
        r, k = min(trials)
        if r > cur * tol:
            break
        keep.remove(k); cur = r
    return keep, cur


def hybrid_windows(Xp, y, theory_wins, wl=WL, sub_nm=30.0):
    """Theory windows (split to <=30 nm) as candidates, pruned by backward iPLS."""
    ivs = []
    for a, b in theory_wins:
        e = np.arange(a, b + sub_nm, sub_nm)
        for lo, hi in zip(e[:-1], np.minimum(e[1:], b + 2)):
            i = np.where((wl >= lo) & (wl < hi))[0]
            if i.size >= 3:
                ivs.append(i)
    chosen, _ = backward_ipls(Xp, y, ivs)
    m = np.zeros(wl.size, bool)
    for c in chosen:
        m[ivs[c]] = True
    return mask_to_windows(m, wl, min_width=0)
