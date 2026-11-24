"""Select usable wavenumber regions (4000-16000 cm^-1) for fructose / allulose
PLS models using a tautomer-resolved (alpha/beta, pyranose/furanose) model.

    python run_ketose.py [--repeats 5] [--ncal 150] [--ntest 100]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sugar_nir.ketose import (K, SUGARS_K, TAUT_NAMES, band_table_k, populations,
                              simulate_ketose, sugar_spectrum, tautomer_spectrum,
                              theory_selectivity_k)
from sugar_nir.pls import PLS1, choose_lv, cv_rmse, preprocess
from sugar_nir.selection import (backward_ipls, forward_ipls, mask_to_windows,
                                 windows_to_mask)
from sugar_nir.simulate import water_spectrum

OUT = Path("results_ketose")
TARGETS = ["fructose", "allulose", "total"]
ALL = TARGETS + ["allulose_pct"]   # pct derived from fructose & allulose predictions
KOR = {"fructose": "과당 (g/100 mL)", "allulose": "알룰로스 (g/100 mL)",
       "total": "총 당 (g/100 mL)", "allulose_pct": "알룰로스 비율 % (두 모델 예측값으로 계산)"}
METHODS = ["full_4000_16000", "nir_4000_10000", "theory_tautomer",
           "theory_dominant_only", "VIP", "iPLS", "hybrid"]
MINW, GAP = 24.0, 32.0          # cm^-1
VALID = (K >= K[0] + 56) & (K <= K[-1] - 56)   # drop Savitzky-Golay edge artefacts


def target(Y, s):
    return {"fructose": Y[:, 0], "allulose": Y[:, 1], "total": Y.sum(1),
            "allulose_pct": 100 * Y[:, 1] / Y.sum(1)}[s]


def theory_wins(sel, top_frac=0.15, smooth=5):
    s = np.convolve(sel, np.ones(smooth) / smooth, mode="same")
    return mask_to_windows(s >= np.quantile(s, 1 - top_frac), K, MINW, GAP)


def vip_wins(X, y):
    lv = choose_lv(cv_rmse(X, y))
    v = np.convolve(PLS1(15).fit(X, y).vip(lv), np.ones(5) / 5, mode="same")
    return mask_to_windows(v > 1.0, K, MINW, GAP)


def intervals(edges_list):
    ivs = []
    for a, b in edges_list:
        i = np.where((K >= a) & (K < b))[0]
        if i.size >= 3:
            ivs.append(i)
    return ivs


def ipls_wins(X, y, width=300.0):
    e = np.arange(K[0], K[-1] + width, width)
    ivs = intervals(zip(e[:-1], e[1:]))
    ch, _ = forward_ipls(X, y, ivs)
    m = np.zeros(K.size, bool)
    for c in ch:
        m[ivs[c]] = True
    return mask_to_windows(m, K, 0, GAP)


def hybrid_wins(X, y, tw, sub=100.0):
    pieces = []
    for a, b in tw:
        e = np.arange(a, b + sub, sub)
        pieces += list(zip(e[:-1], np.minimum(e[1:], b + 8)))
    ivs = intervals(pieces)
    ch, _ = backward_ipls(X, y, ivs)
    m = np.zeros(K.size, bool)
    for c in ch:
        m[ivs[c]] = True
    return mask_to_windows(m, K, 0, GAP)


def score(p, yt):
    return dict(rmsep=float(np.sqrt(np.mean((p - yt) ** 2))),
                r2=float(1 - np.sum((p - yt) ** 2) / np.sum((yt - yt.mean()) ** 2)))


def evaluate(Xc, yc, Xt, yt, mask, Xh=None):
    lv = choose_lv(cv_rmse(Xc[:, mask], yc))
    mdl = PLS1(15).fit(Xc[:, mask], yc)
    p = mdl.predict(Xt[:, mask], lv)
    extra = {"pred": p, "pred_h": mdl.predict(Xh[:, mask], lv)} if Xh is not None else {}
    return extra | dict(rmsep=float(np.sqrt(np.mean((p - yt) ** 2))),
                r2=float(1 - np.sum((p - yt) ** 2) / np.sum((yt - yt.mean()) ** 2)),
                lv=lv, n_wl=int(mask.sum()))


def fmt(w):
    return ", ".join(f"{a:.0f}–{b:.0f}" for a, b in w)


def nm(a, b):
    return f"{1e7 / b:.0f}–{1e7 / a:.0f} nm"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--ncal", type=int, default=150)
    ap.add_argument("--ntest", type=int, default=100)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)

    t0 = time.perf_counter()
    simulate_ketose(30, np.random.default_rng(99))
    t_sim = (time.perf_counter() - t0) / 30

    res = {m: {s: [] for s in ALL} for m in METHODS}
    freq = {m: {s: np.zeros(K.size) for s in TARGETS} for m in METHODS[2:]}
    wins0, sel0 = {}, None
    for rep in range(args.repeats):
        rng = np.random.default_rng(rep)
        Xc, Yc, _ = simulate_ketose(args.ncal, rng)
        Xt, Yt, _ = simulate_ketose(args.ntest, rng)
        Xh, Yh, _ = simulate_ketose(args.ntest, rng, T_mean=35.0)   # temperature-shift test
        Xc, Xt, Xh = preprocess(Xc), preprocess(Xt), preprocess(Xh)
        sel = theory_selectivity_k(np.random.default_rng(1000 + rep))
        seld = theory_selectivity_k(np.random.default_rng(1000 + rep), dominant_only=True)
        if rep == 0:
            sel0 = (sel, seld)
        for s in TARGETS:
            yc, yt = target(Yc, s), target(Yt, s)
            tw = theory_wins(sel[s])
            wins = {"theory_tautomer": tw, "theory_dominant_only": theory_wins(seld[s]),
                    "VIP": vip_wins(Xc, yc), "iPLS": ipls_wins(Xc, yc),
                    "hybrid": hybrid_wins(Xc, yc, tw)}
            masks = {"full_4000_16000": np.ones(K.size, bool), "nir_4000_10000": K <= 10000}
            masks.update({k: windows_to_mask(v, K) & VALID for k, v in wins.items()})
            yh = target(Yh, s)
            for m in METHODS:
                r = evaluate(Xc, yc, Xt, yt, masks[m], Xh)
                rh = score(r["pred_h"], yh)
                r.update(rmsep_35C=rh["rmsep"], r2_35C=rh["r2"])
                res[m][s].append(r)
                if m in freq:
                    freq[m][s] += masks[m]
            if rep == 0:
                wins0[s] = wins
            print(f"rep {rep} {s:13s} " + "  ".join(
                f"{m}:{res[m][s][-1]['rmsep']:.3f}/{res[m][s][-1]['rmsep_35C']:.3f}" for m in METHODS), flush=True)
        for m in METHODS:   # conversion ratio from the two concentration models
            f, a = res[m]["fructose"][-1], res[m]["allulose"][-1]
            pct = lambda pf, pa: 100 * pa / np.clip(pf + pa, 1e-6, None)
            r = score(pct(f["pred"], a["pred"]), target(Yt, "allulose_pct"))
            rh = score(pct(f["pred_h"], a["pred_h"]), target(Yh, "allulose_pct"))
            r.update(rmsep_35C=rh["rmsep"], r2_35C=rh["r2"], lv=0, n_wl=max(f["n_wl"], a["n_wl"]))
            res[m]["allulose_pct"].append(r)
        print(f"rep {rep} allulose_pct  " + "  ".join(
            f"{m}:{res[m]['allulose_pct'][-1]['rmsep']:.2f}" for m in METHODS), flush=True)

    summ = {m: {s: {k: float(np.mean([r[k] for r in res[m][s]])) for k in ("rmsep", "r2", "lv", "n_wl", "rmsep_35C", "r2_35C")}
                | {"rmsep_sd": float(np.std([r["rmsep"] for r in res[m][s]]))}
                for s in ALL} for m in METHODS}
    stable = {s: mask_to_windows(freq["hybrid"][s] >= 0.6 * args.repeats, K, MINW, GAP)
              for s in TARGETS}
    for m in METHODS:
        for s in TARGETS:
            for r in res[m][s]:
                r.pop("pred", None); r.pop("pred_h", None)
    json.dump({"summary": summ, "stable_hybrid_windows_cm-1": stable, "windows_rep0": wins0,
               "sim_seconds_per_spectrum": t_sim}, open(OUT / "results.json", "w"),
              indent=2, default=list)
    plot_tautomers()
    plot_sel(sel0, stable)
    plot_perf(summ)
    report(summ, stable, wins0, t_sim, args)
    print(f"simulation {t_sim * 1e3:.1f} ms/spectrum")


def plot_tautomers():
    fig, ax = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for i, s in enumerate(SUGARS_K):
        p = populations(s)
        for j, t in enumerate(TAUT_NAMES):
            ax[i].plot(K, 10 * p[t] * tautomer_spectrum(s, t), color=f"C{j}",
                       label=f"{t} ({100 * p[t]:.1f} %)")
        ax[i].plot(K, 10 * sugar_spectrum(s), "k", lw=1.5, label="equilibrium sum")
        ax[i].set_title(f"{s}: tautomer contributions at 25 °C (10 g/100 mL)")
        ax[i].legend(fontsize=8); ax[i].set_ylabel("A")
    d = preprocess((10 * (sugar_spectrum("allulose") - sugar_spectrum("fructose")))[None])[0]
    ax[2].plot(K, d, "C3", label="d/dk (allulose − fructose)")
    a2 = ax[2].twinx()
    a2.plot(K, water_spectrum(1e7 / K), "k:", lw=1, label="water")
    a2.set_ylabel("water A")
    ax[2].axhline(0, color="gray", lw=.5); ax[2].legend(fontsize=8, loc="upper right")
    ax[2].set_title("Difference spectrum (1st derivative) — where the two sugars differ")
    ax[2].set_xlabel("wavenumber (cm⁻¹)")
    for a in ax:
        a.grid(alpha=.3); a.set_xlim(K[-1], K[0])
    fig.tight_layout(); fig.savefig(OUT / "tautomer_spectra.png", dpi=130); plt.close(fig)


def plot_sel(sel0, stable):
    sel, seld = sel0
    fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [1.3, 1]})
    for i, s in enumerate(TARGETS):
        ax[0].semilogy(K, sel[s] + 1e-4, color=f"C{i}", label=f"{s} (tautomer model)")
        ax[0].semilogy(K, seld[s] + 1e-4, color=f"C{i}", ls=":", lw=1, label=f"{s} (dominant only)")
    ax[0].set_ylabel("NAS / noise"); ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title("Theory selectivity")
    for i, s in enumerate(TARGETS):
        for a, b in stable[s]:
            ax[1].broken_barh([(a, max(b - a, 8))], (i - .35, .7), color=f"C{i}")
    ax[1].set_yticks(range(len(TARGETS))); ax[1].set_yticklabels(TARGETS)
    ax[1].set_title("Stable selected regions (hybrid, ≥60 % of repeats)")
    ax[1].set_xlabel("wavenumber (cm⁻¹)")
    for a in ax:
        a.grid(alpha=.3); a.set_xlim(K[-1], K[0])
    fig.tight_layout(); fig.savefig(OUT / "selectivity_regions.png", dpi=130); plt.close(fig)


def plot_perf(summ):
    fig, ax = plt.subplots(figsize=(12, 4.5))
    w = 0.11
    for i, m in enumerate(METHODS):
        x = np.arange(len(ALL)) + (i - 3) * w
        ax.bar(x, [summ[m][s]["rmsep"] / summ["full_4000_16000"][s]["rmsep"] for s in ALL],
               w, color=f"C{i}", label=f"{m}")
        ax.scatter(x, [summ[m][s]["rmsep_35C"] / summ["full_4000_16000"][s]["rmsep"] for s in ALL],
                   color=f"C{i}", marker="v", edgecolor="k", zorder=3, s=25)
    ax.axhline(1, color="k", lw=.8)
    ax.set_xticks(range(len(ALL))); ax.set_xticklabels(ALL)
    ax.set_ylabel("RMSEP / RMSEP(full, 25 °C)"); ax.legend(fontsize=8, ncol=4)
    ax.grid(axis="y", alpha=.3)
    ax.set_title("Bars: test at 25 °C   ▼: test at 35 °C (tautomer populations shifted)")
    fig.tight_layout(); fig.savefig(OUT / "performance.png", dpi=130); plt.close(fig)


def assign(a, b, table):
    hits = [r for r in table if a - 25 <= r[2] <= b + 25]
    if not hits:
        mid = (a + b) / 2
        r = min(table, key=lambda r: abs(r[2] - mid))
        return f"{r[4]} ({r[2]:.0f} cm⁻¹) 밴드의 경사면"
    hits = sorted(hits, key=lambda r: -r[3])[:4]
    return "; ".join(r[4] for r in hits)


def report(summ, stable, wins0, t_sim, args):
    tab = band_table_k()
    L = ["# 과당 / 알룰로스 영역 선택 결과 (4000–16000 cm⁻¹, 호변이성질체 모델)\n",
         f"- 반복 {args.repeats}회, 검량 {args.ncal}개 / 테스트 {args.ntest}개 (독립)",
         f"- 시뮬레이션 비용: 스펙트럼당 {t_sim * 1e3:.0f} ms",
         "- 시료: 총 당 5–40 g/100 mL, 알룰로스 비율 0–100 %, 온도 25±2 °C (호변이성질체 비율이 온도에 따라 변함)\n",
         "## 호변이성질체 비율 (25 °C, 문헌 NMR 근사)\n",
         "| | β-피라노스 | α-피라노스 | β-퓨라노스 | α-퓨라노스 | 케토 |", "|---|---|---|---|---|---|"]
    for s in SUGARS_K:
        p = populations(s)
        L.append(f"| {s} | " + " | ".join(f"{100 * p[t]:.1f} %" for t in TAUT_NAMES) + " |")
    L += ["\n## 1-b. 온도 이동 테스트 (검량 25 °C → 테스트 35 °C, RMSEP / R²)\n",
          "| 방법 | " + " | ".join(KOR[s] for s in ALL) + " |", "|---|" + "---|" * len(ALL)]
    for m in METHODS:
        L.append(f"| {m} | " + " | ".join(
            f"{summ[m][s]['rmsep_35C']:.3f} / {summ[m][s]['r2_35C']:.3f}" for s in ALL) + " |")
    L += ["\n## 1. 성능 (25 °C 테스트 RMSEP ± SD / R² / 사용 포인트 수)\n",
          "| 방법 | " + " | ".join(KOR[s] for s in ALL) + " |", "|---|" + "---|" * len(ALL)]
    for m in METHODS:
        L.append(f"| {m} | " + " | ".join(
            f"{summ[m][s]['rmsep']:.3f}±{summ[m][s]['rmsep_sd']:.3f} / {summ[m][s]['r2']:.3f} / "
            f"{summ[m][s]['n_wl']:.0f}" for s in ALL) + " |")
    L.append("\n## 2. 추천 영역 (hybrid, 반복 60 % 이상에서 선택)\n")
    um = np.zeros(K.size, bool)
    for s in TARGETS:
        L.append(f"### {KOR[s]}")
        for a, b in stable[s]:
            L.append(f"- **{a:.0f}–{b:.0f} cm⁻¹** ({nm(a, b)}): {assign(a, b, tab)}")
        um |= windows_to_mask(stable[s], K)
        L.append("")
    L += ["### 공통 합집합", fmt(mask_to_windows(um, K, 0, GAP)) +
          f" cm⁻¹ (전체 포인트의 {100 * um.mean():.0f} %)\n", "## 3. 방법별 영역 (반복 0, cm⁻¹)\n"]
    for s in TARGETS:
        L.append(f"**{KOR[s]}**")
        for m, w in wins0[s].items():
            L.append(f"- {m}: {fmt(w)}")
        L.append("")
    (OUT / "report.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
