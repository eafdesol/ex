"""Group-contribution band model for glucose and malto-oligosaccharides (DP1-DP4).

Each sugar is described by its X-H oscillators (counted per molecule).
For every oscillator type the Morse solver gives overtone positions and
intensities; combination bands (stretch + bend / stretch + C-O stretch) are
added with empirical positions, as in standard NIR band assignment tables.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .morse import morse_transitions


@dataclass
class Group:
    name: str
    we: float          # harmonic wavenumber (cm^-1)
    wexe: float        # anharmonicity (cm^-1)
    strength: float    # dipole-derivative scale
    sigma1: float      # Gaussian sigma (cm^-1) of the first overtone
    # (partner frequency cm^-1, intensity relative to 1st overtone, sigma, label)
    combos: list = field(default_factory=list)


# H-bonded oscillators in aqueous solution
GROUPS = {
    "OH_sec":  Group("OH_sec",  3560, 90, 1.00, 120,
                     [(1380, 1.6, 110, "OH str + COH bend"),
                      (1050, 0.6, 90, "OH str + CO str")]),
    "OH_prim": Group("OH_prim", 3575, 90, 1.00, 120,
                     [(1360, 1.6, 110, "OH str + COH bend"),
                      (1040, 0.6, 90, "OH str + CO str")]),
    "OH_anom": Group("OH_anom", 3545, 91, 1.00, 125,
                     [(1405, 1.6, 110, "anomeric OH str + bend"),
                      (1075, 0.7, 90, "anomeric OH + C1-O str")]),
    "CH_ring": Group("CH_ring", 3000, 58, 0.55, 55,
                     [(1440, 1.3, 45, "CH str + CH bend"),
                      (1110, 0.5, 50, "CH str + CO str"),
                      (4207, 0.12, 70, "2CH str + CH bend")]),
    "CH2":     Group("CH2",     2995, 60, 0.55, 50,
                     [(1465, 1.3, 40, "CH2 str + scissor"),
                      (1340, 0.7, 45, "CH2 str + wag")]),
    "CH_anom": Group("CH_anom", 3045, 60, 0.55, 55,
                     [(1400, 1.2, 45, "anomeric CH str + bend"),
                      (1090, 0.6, 50, "anomeric CH + C1-O str")]),
    "CH_link": Group("CH_link", 3035, 60, 0.55, 55,
                     [(1395, 1.2, 45, "C4-H(linked) str + bend"),
                      (1140, 0.8, 55, "C4-H + C-O-C str")]),
    "CH_glyc": Group("CH_glyc", 3055, 61, 0.55, 55,
                     [(1385, 1.2, 45, "glycosidic CH str + bend"),
                      (1150, 0.9, 55, "CH + C-O-C asym str")]),
}


@dataclass
class Sugar:
    name: str
    kind: str
    mw: float
    groups: dict
    oh_shift: float = 0.0   # extra intramolecular H-bond shift of OH (cm^-1)


# Glucose-only alpha-(1->4) chain (malto-oligosaccharides), all reducing sugars.
# Per n glucose units: 1 reducing end (anomeric OH + anomeric CH),
# n-1 glycosidic C1-H, n-1 linked C4-H (next to the glycosidic O).
SUGARS = {
    "glucose":       Sugar("glucose", "DP1 monosaccharide", 180.16,
                           {"OH_anom": 1, "OH_sec": 3, "OH_prim": 1,
                            "CH_anom": 1, "CH_ring": 4, "CH2": 2}, 0.0),
    "maltose":       Sugar("maltose", "DP2 disaccharide", 342.30,
                           {"OH_anom": 1, "OH_sec": 5, "OH_prim": 2,
                            "CH_anom": 1, "CH_glyc": 1, "CH_link": 1,
                            "CH_ring": 7, "CH2": 4}, -6.0),
    "maltotriose":   Sugar("maltotriose", "DP3 trisaccharide", 504.44,
                           {"OH_anom": 1, "OH_sec": 7, "OH_prim": 3,
                            "CH_anom": 1, "CH_glyc": 2, "CH_link": 2,
                            "CH_ring": 10, "CH2": 6}, -9.0),
    "maltotetraose": Sugar("maltotetraose", "DP4 tetrasaccharide", 666.58,
                           {"OH_anom": 1, "OH_sec": 9, "OH_prim": 4,
                            "CH_anom": 1, "CH_glyc": 3, "CH_link": 3,
                            "CH_ring": 13, "CH2": 8}, -11.0),
}
SUGAR_NAMES = list(SUGARS)


def group_bands(g: Group, d_we: float = 0.0, width_scale: float = 1.0,
                n_levels: int = 7):
    """List of (center cm^-1, intensity, sigma cm^-1, label) for one oscillator."""
    nu, rel = morse_transitions(round(g.we + d_we, 1), g.wexe, n_levels)
    bands = []
    i_ot1 = g.strength * rel[1]
    for v in range(2, n_levels):            # overtones (fundamental is < 4000)
        sigma = g.sigma1 * (v / 2) ** 0.8 * width_scale
        bands.append((nu[v - 1], g.strength * rel[v - 1], sigma,
                      f"{g.name} {v - 1}{'st' if v == 2 else 'nd' if v == 3 else 'rd' if v == 4 else 'th'} overtone"))
    for partner, r, sig, label in g.combos:
        bands.append((nu[0] + partner - 15.0, i_ot1 * r, sig * width_scale, label))
    return bands


def sugar_bands(s: Sugar, perturb=None, oh_conc_shift: float = 0.0):
    """All bands of one sugar molecule.

    perturb: optional dict {group_name: (d_we, width_scale, intensity_scale)}
             used to emulate theory-vs-reality model error.
    """
    out = []
    for gname, count in s.groups.items():
        g = GROUPS[gname]
        d_we, ws, isc = (0.0, 1.0, 1.0)
        if perturb is not None and gname in perturb:
            d_we, ws, isc = perturb[gname]
        if gname.startswith("OH"):
            d_we += s.oh_shift + oh_conc_shift
        for c, I, sig, lab in group_bands(g, d_we, ws):
            out.append((c, count * I * isc, sig, f"{s.name}: {lab}"))
    return out


def render(bands, wl_nm: np.ndarray) -> np.ndarray:
    """Sum of Gaussians (in wavenumber space) evaluated on a wavelength grid."""
    k = 1e7 / wl_nm
    A = np.zeros_like(wl_nm, dtype=float)
    for c, I, sig, _ in bands:
        A += I * np.exp(-0.5 * ((k - c) / sig) ** 2)
    return A


def band_table(wl_min=400, wl_max=2500, min_rel=1e-3):
    """Theoretical band assignment table (nm) for all sugars."""
    rows = []
    for s in SUGARS.values():
        b = sugar_bands(s)
        imax = max(x[1] for x in b)
        for c, I, sig, lab in b:
            nm = 1e7 / c
            if wl_min <= nm <= wl_max and I / imax >= min_rel:
                rows.append((s.name, nm, I / imax, lab))
    return sorted(rows, key=lambda r: r[1])
