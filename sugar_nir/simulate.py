"""Synthetic Vis-NIR (400-2500 nm) spectra of aqueous sugar mixtures.

Realism ingredients
  * sugar bands from the Morse/group model (per gram, i.e. divided by MW)
  * water matrix bands + water displacement by dissolved sugar
  * temperature-dependent water band shift
  * concentration-dependent OH H-bond shift (non-linearity)
  * visible-range pigment interferent (browning) unrelated to targets
  * additive / sloped / scatter baselines and path-length multiplicative error
  * heteroscedastic noise: sigma ~ 10^A (fewer photons where absorbance is high)
"""
from __future__ import annotations

import numpy as np

from .bands import SUGARS, SUGAR_NAMES, render, sugar_bands

WL = np.arange(400.0, 2500.0 + 1e-9, 2.0)      # nm
SUGAR_SCALE = 0.9      # absorbance scale for 1 g/100 mL (per MW-normalised count)
PSV = 0.0062           # volume fraction displaced per g/100 mL


def water_bands(T: float = 25.0):
    dT = T - 25.0
    return [
        (10300 + 1.0 * dT, 0.035, 300, "water 2nd OT (970 nm)"),
        (8500 + 1.0 * dT, 0.055, 280, "water comb (1190 nm)"),
        (6900 + 1.5 * dT, 1.00 * (1 - 0.003 * dT), 220, "water nu1+nu3 (1450 nm)"),
        (5600, 0.06, 200, "water (1790 nm)"),
        (5180 + 1.0 * dT, 1.90 * (1 - 0.002 * dT), 140, "water nu2+nu3 (1940 nm)"),
        (4700, 0.08, 250, "water nu2+lib (2130 nm)"),
        (3400 + 1.0 * dT, 60.0, 210, "water nu1/nu3 fundamental tail"),
    ]


def water_spectrum(wl=WL, T=25.0):
    return render(water_bands(T), wl)


def pure_sugar_spectrum(name, wl=WL, perturb=None, oh_conc_shift=0.0):
    s = SUGARS[name]
    return SUGAR_SCALE * render(sugar_bands(s, perturb, oh_conc_shift), wl) / s.mw


def pigment_spectrum(wl=WL, center=470.0, width=60.0):
    return np.exp(-0.5 * ((wl - center) / width) ** 2) + 0.3 * (450.0 / wl) ** 4


def noise_sigma(A, base=2.5e-4):
    return base * 10 ** np.clip(A, 0, 3.2)


def simulate_dataset(n, rng, cmax=10.0, wl=WL):
    """Return X (n x p absorbance), Y (n x 4, g/100 mL), meta."""
    Y = rng.uniform(0, cmax, size=(n, len(SUGAR_NAMES)))
    X = np.empty((n, wl.size))
    for i in range(n):
        T = rng.normal(25.0, 1.5)
        ctot = Y[i].sum()
        oh_shift = -0.5 * ctot
        A = water_spectrum(wl, T) * (1 - PSV * ctot)
        for j, name in enumerate(SUGAR_NAMES):
            A += Y[i, j] * pure_sugar_spectrum(name, wl, oh_conc_shift=oh_shift)
        A += rng.uniform(0, 0.4) * pigment_spectrum(wl, rng.uniform(430, 520),
                                                    rng.uniform(40, 80))
        A *= rng.normal(1.0, 0.02)                          # path length
        x = wl / 1000.0
        A += rng.normal(0, 0.02) + rng.normal(0, 0.01) * x \
            + rng.uniform(0, 0.03) * (1 / x) ** 2           # offset/slope/scatter
        A += rng.normal(0, 1, wl.size) * noise_sigma(A)
        X[i] = A
    return X, Y
