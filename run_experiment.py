"""Find discriminative 400-2500 nm regions for glucose / maltose / maltotriose /
maltotetraose (glucose-only alpha-1,4 chains) and compare PLS performance: full spectrum vs selected regions.

    python run_experiment.py [--repeats 5] [--ncal 150] [--ntest 100]
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

from sugar_nir.bands import SUGARS, SUGAR_NAMES, band_table
from sugar_nir.pls import PLS1, choose_lv, cv_rmse, preprocess
from sugar_nir.selection import (hybrid_windows, mask_to_windows, ipls_windows, theory_selectivity,
                                 theory_windows, vip_windows, windows_to_mask)
from sugar_nir.simulate import WL, pure_sugar_spectrum, simulate_dataset, water_spectrum

OUT = Path("results")
METHODS = ["full_400_2500", "nir_1000_2500", "theory_NAS", "VIP", "iPLS", "hybrid"]
KOR = {"glucose": "포도당(DP1)", "maltose": "말토오스(DP2)",
       "maltotriose": "말토트리오스(DP3)", "maltotetraose": "말토테트라오스(DP4)",
       "DP3_4_sum": "DP3+DP4 합", "total": "총 당량 (g/100 mL)",
       "reducing_ends": "환원 말단 (mmol/100 mL)", "mean_DP": "평균 중합도 DP"}
TARGETS = SUGAR_NAMES + ["DP3_4_sum", "total", "reducing_ends", "mean_DP"]
MW = np.array([SUGARS[n].mw for n in SUGAR_NAMES])


def target(Y, s):
    if s == "DP3_4_sum":
        return Y[:, 2] + Y[:, 3]
    if s == "total":
        return Y.sum(1)
    if s == "reducing_ends":
        return (Y / MW * 1000).sum(1)
    if s == "mean_DP":   # number-average DP from Mn = 162.14 DP + 18.02
        return (Y.sum(1) / (Y / MW).sum(1) - 18.02) / 162.14
    return Y[:, SUGAR_NAMES.index(s)]


def evaluate(Xc, yc, Xt, yt, mask):
    rm = cv_rmse(Xc[:, mask], yc)
    lv = choose_lv(rm)
    pred = PLS1(15).fit(Xc[:, mask], yc).predict(Xt[:, mask], lv)
    rmsep = float(np.sqrt(np.mean((pred - yt) ** 2)))
    r2 = float(1 - np.sum((pred - yt) ** 2) / np.sum((yt - yt.mean()) ** 2))
    return dict(rmsep=rmsep, r2=r2, lv=lv, n_wl=int(mask.sum()))


def fmt_wins(w):
    return ", ".join(f"{a:.0f}–{b:.0f}" for a, b in w)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--ncal", type=int, default=150)
    ap.add_argument("--ntest", type=int, default=100)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)

    t0 = time.perf_counter()
    simulate_dataset(50, np.random.default_rng(99))
    t_sim = (time.perf_counter() - t0) / 50

    res = {m: {s: [] for s in TARGETS} for m in METHODS}
    freq = {m: {s: np.zeros(WL.size) for s in TARGETS} for m in METHODS[2:]}
    wins0, sel0 = {}, None
    for rep in range(args.repeats):
        rng = np.random.default_rng(rep)
        Xc, Yc = simulate_dataset(args.ncal, rng)
        Xt, Yt = simulate_dataset(args.ntest, rng)
        Xc, Xt = preprocess(Xc), preprocess(Xt)
        sel = theory_selectivity(np.random.default_rng(1000 + rep))
        if rep == 0:
            sel0 = sel
        for s in TARGETS:
            yc, yt = target(Yc, s), target(Yt, s)
            tw = theory_windows(sel[s])
            wins = {
                "theory_NAS": tw,
                "VIP": vip_windows(Xc, yc),
                "iPLS": ipls_windows(Xc, yc),
                "hybrid": hybrid_windows(Xc, yc, tw),
            }
            masks = {"full_400_2500": np.ones(WL.size, bool),
                     "nir_1000_2500": WL >= 1000}
            masks.update({k: windows_to_mask(v) for k, v in wins.items()})
            for m in METHODS:
                r = evaluate(Xc, yc, Xt, yt, masks[m])
                res[m][s].append(r)
                if m in freq:
                    freq[m][s] += masks[m]
            if rep == 0:
                wins0[s] = wins
            print(f"rep {rep} {s:10s} " + "  ".join(
                f"{m}:{res[m][s][-1]['rmsep']:.3f}" for m in METHODS), flush=True)

    # ---------- summary ----------
    summ = {m: {s: {k: float(np.mean([r[k] for r in res[m][s]]))
                    for k in ("rmsep", "r2", "lv", "n_wl")} |
                {"rmsep_sd": float(np.std([r["rmsep"] for r in res[m][s]]))}
                for s in TARGETS} for m in METHODS}
    # stable windows: selected in >= 60 % of repeats by the hybrid method
    stable = {}
    for s in TARGETS:
        stable[s] = mask_to_windows(freq["hybrid"][s] >= 0.6 * args.repeats, WL,
                                    min_width=0)
    json.dump({"summary": summ, "stable_hybrid_windows": stable,
               "windows_rep0": wins0, "sim_seconds_per_spectrum": t_sim},
              open(OUT / "results.json", "w"), indent=2, default=list)

    plot_spectra(sel0, wins0, stable)
    plot_rmsep(summ)
    plot_freq(freq, args.repeats)
    write_report(summ, wins0, stable, t_sim, args)
    print(f"\nsimulation: {t_sim * 1e3:.1f} ms per spectrum")


def plot_spectra(sel, wins0, stable):
    fig, ax = plt.subplots(3, 1, figsize=(11, 12), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1.4]})
    for s in SUGAR_NAMES:
        ax[0].plot(WL, pure_sugar_spectrum(s) * 10, label=f"{s} (10 g/100 mL)")
    ax0b = ax[0].twinx()
    ax0b.plot(WL, water_spectrum(), "k:", lw=1, label="water")
    ax0b.set_ylabel("water A")
    ax[0].set_ylabel("sugar A"); ax[0].legend(fontsize=8, loc="upper left")
    ax[0].set_title("Simulated pure-component spectra (Morse/group model)")
    for s in TARGETS:
        ax[1].semilogy(WL, sel[s] + 1e-3, label=s, ls="--" if s == "DP3_4_sum" else "-")
    ax[1].set_ylabel("NAS / noise (per g/100 mL)"); ax[1].legend(fontsize=8)
    ax[1].set_title("Theory selectivity (net analyte signal vs noise)")
    for k, s in enumerate(TARGETS):
        for a, b in stable[s]:
            ax[2].broken_barh([(a, max(b - a, 2))], (k - 0.35, 0.7), color=f"C{k}")
    ax[2].set_yticks(range(len(TARGETS))); ax[2].set_yticklabels(TARGETS)
    ax[2].set_xlabel("wavelength (nm)")
    ax[2].set_title("Stable selected regions (hybrid, chosen in ≥60 % of repeats)")
    for a in ax:
        a.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / "spectra_selectivity_regions.png", dpi=130)
    plt.close(fig)


def plot_rmsep(summ):
    fig, ax = plt.subplots(figsize=(13, 4.5))
    w = 0.13
    nt = len(TARGETS)
    for i, m in enumerate(METHODS):
        v = [1 - summ[m][s]["r2"] for s in TARGETS]
        e = None
        ax.bar(np.arange(nt) + (i - 2.5) * w, v, w, yerr=e, label=m, capsize=2)
    ax.set_xticks(range(nt)); ax.set_xticklabels(TARGETS, fontsize=8)
    ax.set_ylabel("1 - R² (test, lower is better)"); ax.legend(fontsize=8, ncol=3)
    ax.set_title("PLS test error: full spectrum vs selected regions"); ax.set_yscale("log")
    ax.grid(axis="y", alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / "rmsep_comparison.png", dpi=130)
    plt.close(fig)


def plot_freq(freq, R):
    fig, ax = plt.subplots(len(TARGETS), 1, figsize=(11, 14), sharex=True)
    for k, s in enumerate(TARGETS):
        for m in ("theory_NAS", "VIP", "iPLS", "hybrid"):
            ax[k].plot(WL, freq[m][s] / R, label=m)
        ax[k].set_ylabel(s); ax[k].set_ylim(-0.05, 1.05); ax[k].grid(alpha=.3)
    ax[0].legend(fontsize=8, ncol=4); ax[-1].set_xlabel("wavelength (nm)")
    ax[0].set_title("Selection frequency across repeats")
    fig.tight_layout(); fig.savefig(OUT / "selection_frequency.png", dpi=130)
    plt.close(fig)


def assign(a, b, table):
    labs = sorted({r[3].split(": ", 1)[1] for r in table if a - 15 <= r[1] <= b + 15})
    if labs:
        return "; ".join(labs)
    mid = 0.5 * (a + b)
    near = min(table, key=lambda r: abs(r[1] - mid))
    side = "장파장" if mid > near[1] else "단파장"
    return (f"{near[3].split(': ', 1)[1]} ({near[1]:.0f} nm) 밴드의 {side} 쪽 경사면 "
            "— 물 밴드와 겹치지 않는 미분 신호")


def write_report(summ, wins0, stable, t_sim, args):
    table = band_table(min_rel=0.02)
    L = ["# 당류 판별용 파장 영역 선택 결과\n",
         f"- 반복 {args.repeats}회, 검량 {args.ncal}개 / 테스트 {args.ntest}개 (독립 세트)",
         f"- 시뮬레이션 비용: 스펙트럼 1개당 {t_sim * 1e3:.1f} ms (Morse 슈뢰딩거 수치해 + 밴드 합성)",
         "- 농도: 4종 각 0–10 g/100 mL 독립 랜덤, 수용액 (물 매트릭스, 온도 ±1.5 °C, 색소 간섭, 산란/광로 변동, 이분산 노이즈)\n",
         "## 1. 성능 비교 (RMSEP 평균 ± SD / R² / 파장 수; 단위는 타깃별)\n",
         "| 방법 | " + " | ".join(KOR[s] for s in TARGETS) + " |",
         "|---|" + "---|" * len(TARGETS)]
    for m in METHODS:
        L.append(f"| {m} | " + " | ".join(
            f"{summ[m][s]['rmsep']:.3f}±{summ[m][s]['rmsep_sd']:.3f} / "
            f"{summ[m][s]['r2']:.3f} / {summ[m][s]['n_wl']:.0f}" for s in TARGETS) + " |")
    L += ["\n## 2. 당류별 추천 영역 (hybrid, 반복의 60 % 이상에서 선택된 안정 영역)\n"]
    for s in TARGETS:
        L.append(f"### {KOR[s]} — {s}")
        for a, b in stable[s]:
            L.append(f"- **{a:.0f}–{b:.0f} nm**: {assign(a, b, table)}")
        L.append("")
    um = np.zeros(WL.size, bool)
    for s in TARGETS:
        um |= windows_to_mask(stable[s])
    allw = mask_to_windows(um, WL, min_width=0)
    L += ["### 공통 모델(PLS2 등)용 합집합", fmt_wins(allw) + f"  (총 {um.sum()}개 파장, 전체의 {100 * um.mean():.0f} %)\n",
          "## 3. 방법별 선택 영역 (반복 0)\n"]
    for s in TARGETS:
        L.append(f"**{KOR[s]}**")
        for m, w in wins0[s].items():
            L.append(f"- {m}: {fmt_wins(w)}")
        L.append("")
    L += ["## 그림", "- `spectra_selectivity_regions.png`", "- `rmsep_comparison.png`",
          "- `selection_frequency.png`"]
    (OUT / "report.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
