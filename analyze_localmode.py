"""Analyse C-H local-mode results: per-type shifts and glucose vs maltose
first-overtone (and second-overtone) difference spectra.

    python analyze_localmode.py [xtb|dft]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

from sugar_nir.localmode import SMILES, label_atoms, ring_orientation

OUT = Path("results_localmode")
MW = {"glucose": 180.16, "maltose": 342.30}
# anomer populations in water at ~25 C (approx. literature)
POP = {"glucose": {"a-glucose": 0.36, "b-glucose": 0.64},
       "maltose": {"a-maltose": 0.42, "b-maltose": 0.58}}
FWHM = {2: 60.0, 3: 90.0}       # per-C-H band width in solution (cm^-1), assumed


def dft_pairs():
    rows = []
    for f in sorted(OUT.glob("ch_overtones_dft*.csv")):
        rows += list(csv.DictReader(open(f)))
    return rows


def calibration():
    """Linear map xTB -> DFT for nu2 and nu3 from the paired scans."""
    P = dft_pairs()
    cal = {}
    for v in (2, 3):
        x = np.array([float(r[f"xtb_nu{v}"]) for r in P])
        y = np.array([float(r[f"dft_nu{v}"]) for r in P])
        b, a = np.polyfit(x, y, 1)
        res = y - (a + b * x)
        cal[v] = (a, b, float(np.sqrt(np.mean(res ** 2))), len(x))
    return cal, P


def load(method):
    rows = list(csv.DictReader(open(OUT / "ch_overtones_xtb.csv")))
    for r in rows:
        for k in ("we", "wexe", "nu1", "nu2", "nu3", "I1", "I2", "I3"):
            r[k] = float(r[k])
    if method == "cal":
        cal, _ = calibration()
        for r in rows:
            for v in (2, 3):
                a, b, _, _ = cal[v]
                r[f"nu{v}"] = a + b * r[f"nu{v}"]
    return rows


def h_orientation():
    """Axial/equatorial orientation of every ring C-H (lowest conformer)."""
    geoms = json.load(open(OUT / "geometries.json"))
    out = {}
    for name, g in geoms.items():
        pos = np.array(g["pos"][0])
        mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
        for h, d in label_atoms(mol).items():
            out[(name, d["label"])] = ("CH2" if d["pos"] == 6 else
                                       ring_orientation(pos, d["ring"], d["carbon"], h))
    return out


def spectrum(rows, sugar, v, k):
    A = np.zeros_like(k)
    sig = FWHM[v] / 2.3548
    for name, p in POP[sugar].items():
        for r in rows:
            if r["molecule"] == name:
                A += p * r[f"I{v}"] * np.exp(-0.5 * ((k - r[f"nu{v}"]) / sig) ** 2)
    return A / MW[sugar]


def windows(mask, k, min_w=20):
    out, s = [], None
    for i, m in enumerate(mask):
        if m and s is None:
            s = i
        if (not m or i == len(mask) - 1) and s is not None:
            e = i if m else i - 1
            if k[e] - k[s] >= min_w:
                out.append((k[s], k[e]))
            s = None
    return out


def main(method="xtb"):
    rows = load(method)
    ori = h_orientation()
    for r in rows:
        r["orient"] = ori.get((r["molecule"], r["label"]), "?")
    ref = np.mean([r["nu2"] for r in rows if "glucose" in r["molecule"]])
    title = {"xtb": "GFN2-xTB 원값", "cal": "GFN2-xTB → B3LYP/6-31G** 보정"}[method]
    L = [f"# C-H local-mode 배음 계산 ({title}, 물 용매 모델)\n",
         f"포도당 전체 C-H 1차 배음 평균 (계산 원값): {ref:.0f} cm⁻¹ ({1e7 / ref:.0f} nm)\n",
         "Δ = 해당 C-H의 1차 배음 − 포도당 평균 (cm⁻¹). 절대값보다 Δ가 더 믿을 만함.\n",
         "## 1. C-H별 결과\n",
         "| 분자 | C-H | 유형 | 방향 | ωe | ωexe | 1차 배음 cm⁻¹ | Δ | 2차 배음 cm⁻¹ |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['molecule']} | {r['label']} | {r['cls']} | {r['orient']} | {r['we']:.0f} | "
                 f"{r['wexe']:.1f} | {r['nu2']:.0f} | {r['nu2'] - ref:+.0f} | {r['nu3']:.0f} |")
    # class averages
    L += ["\n## 2. 유형별 평균 (모든 분자, Δ cm⁻¹, ± 표준편차, 개수)\n",
          "| 유형 | 축 C-H | 적도 C-H | 전체 |", "|---|---|---|---|"]
    classes = sorted({r["cls"] for r in rows})
    stats = {}
    for c in classes:
        cells = []
        for o in ("ax", "eq", None):
            v = [r["nu2"] - ref for r in rows if r["cls"] == c and (o is None or r["orient"] == o)]
            cells.append(f"{np.mean(v):+.0f} ± {np.std(v):.0f} (n={len(v)})" if v else "–")
            if o is None:
                stats[c] = np.mean(v)
        L.append(f"| {c} | " + " | ".join(cells) + " |")
    for o in ("ax", "eq", "CH2"):
        v = [r["nu2"] - ref for r in rows if r["orient"] == o]
        L.append(f"| (방향만) {o} | | | {np.mean(v):+.0f} ± {np.std(v):.0f} (n={len(v)}) |")

    # spectra
    k = np.arange(5000.0, 6400.0, 1.0)
    k3 = np.arange(7500.0, 9500.0, 1.0)
    g2, m2 = spectrum(rows, "glucose", 2, k), spectrum(rows, "maltose", 2, k)
    g3, m3 = spectrum(rows, "glucose", 3, k3), spectrum(rows, "maltose", 3, k3)
    n2 = max(g2.max(), m2.max())
    d2 = (m2 - g2) / n2
    dd2 = np.gradient(m2, k) - np.gradient(g2, k)
    # maltose-specific C-H (glycosidic) contribution
    spec_lab = {"non-red. C1-H", "red. C4-H"}
    sig = FWHM[2] / 2.3548
    glyc = np.zeros_like(k)
    for name, p in POP["maltose"].items():
        for r in rows:
            if r["molecule"] == name and r["label"] in spec_lab:
                glyc += p * r["I2"] * np.exp(-0.5 * ((k - r["nu2"]) / sig) ** 2)
    glyc /= MW["maltose"]
    thr = 0.5 * np.abs(d2).max()
    wins_pos = windows(d2 > thr, k)
    wins_neg = windows(d2 < -thr, k)
    L += ["\n## 3. 포도당 vs 말토오스 (그램당, 아노머 비율 반영, C-H당 FWHM 60 cm⁻¹ 가정)\n",
          f"- 말토오스가 더 큰 구간 (차이 > 최대의 50 %): " +
          ", ".join(f"{a:.0f}–{b:.0f} cm⁻¹ (Δ {a - ref:+.0f}~{b - ref:+.0f})" for a, b in wins_pos),
          f"- 포도당이 더 큰 구간: " +
          ", ".join(f"{a:.0f}–{b:.0f} cm⁻¹ (Δ {a - ref:+.0f}~{b - ref:+.0f})" for a, b in wins_neg),
          f"- 말토오스 글리코시드 C-H(비환원 C1-H, 환원 C4-H) 밴드 중심: " +
          ", ".join(f"{r['molecule']} {r['label']} {r['nu2']:.0f} (Δ {r['nu2'] - ref:+.0f})"
                    for r in rows if r["label"] in spec_lab),
          f"- 최대 상대 차이 |Δ스펙트럼|/최대흡광: {np.abs(d2).max():.3f}"]
    if method == "cal":
        cal, P = calibration()
        a, b, rmse, n = cal[2]
        L += ["\n## 4. xTB → DFT 보정 (1차 배음)\n",
              f"DFT = {a:.0f} + {b:.3f} × xTB  (쌍 {n}개, 잔차 RMSE {rmse:.0f} cm⁻¹)\n",
              "| 분자 | C-H | 유형 | xTB | DFT | 보정 xTB |", "|---|---|---|---|---|---|"]
        for r in P:
            x, y = float(r["xtb_nu2"]), float(r["dft_nu2"])
            L.append(f"| {r['molecule']} | {r['label']} | {r['cls']} | {x:.0f} | {y:.0f} | {a + b * x:.0f} |")
        L += ["\n## 5. 파장으로 환산한 결과 (보정값)\n"]
        for lab, ws in (("말토오스가 더 큰 구간", wins_pos), ("포도당이 더 큰 구간", wins_neg)):
            L.append(f"- {lab}: " + ", ".join(f"{a_:.0f}–{b_:.0f} cm⁻¹ = {1e7 / b_:.0f}–{1e7 / a_:.0f} nm"
                                              for a_, b_ in ws))
    (OUT / f"report_{method}.md").write_text("\n".join(L), encoding="utf-8")

    fig, ax = plt.subplots(3, 1, figsize=(11, 11))
    ax[0].plot(k, g2, label="glucose (α36/β64)")
    ax[0].plot(k, m2, label="maltose (α42/β58)")
    ax[0].plot(k, glyc, "C3--", label="maltose glycosidic C-H only")
    for name in ("b-glucose", "b-maltose"):
        for r in rows:
            if r["molecule"] == name:
                ax[0].axvline(r["nu2"], ymax=0.12 if "glucose" in name else 0.08,
                              color="C0" if "glucose" in name else "C1", lw=1)
    ax[0].set_title("C-H 1st overtone (per gram); ticks: individual C-H (β anomers)")
    ax[0].legend(fontsize=8); ax[0].set_xlim(k[-1], k[0])
    ax[1].plot(k, d2, "k", label="(maltose − glucose) / max")
    ax[1].plot(k, dd2 / np.abs(dd2).max() * np.abs(d2).max(), "C2", lw=1,
               label="1st-derivative difference (scaled)")
    for a, b in wins_pos:
        ax[1].axvspan(a, b, color="C1", alpha=.2)
    for a, b in wins_neg:
        ax[1].axvspan(a, b, color="C0", alpha=.2)
    ax[1].axhline(0, color="gray", lw=.5); ax[1].legend(fontsize=8); ax[1].set_xlim(k[-1], k[0])
    ax[1].set_title("Difference spectrum: orange = maltose-rich, blue = glucose-rich")
    ax[2].plot(k3, g3, label="glucose"); ax[2].plot(k3, m3, label="maltose")
    ax[2].set_title("C-H 2nd overtone (per gram)"); ax[2].legend(fontsize=8)
    ax[2].set_xlim(k3[-1], k3[0])
    ax[2].set_xlabel("wavenumber (cm⁻¹, " + ("xTB raw" if method == "xtb" else "calibrated to DFT") + ")")
    for a in ax:
        a.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / f"overtone_{method}.png", dpi=130); plt.close(fig)
    print("\n".join(L[-8:]))


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["xtb"]))
