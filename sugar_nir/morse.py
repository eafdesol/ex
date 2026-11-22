"""1D Morse oscillator solved numerically from the Schrodinger equation.

Hamiltonian in the dimensionless coordinate y = a (r - r_e):

    H = -B d^2/dy^2 + D_e (1 - exp(-y))^2

With B = omega_e x_e and D_e = omega_e^2 / (4 omega_e x_e) (all in cm^-1) the
eigenvalues reproduce the textbook Morse term values

    G(v) = omega_e (v + 1/2) - omega_e x_e (v + 1/2)^2

The eigenfunctions are used to compute transition dipole moments
<0|mu|v> with a Mecke-type dipole function mu(y) = y exp(-y / y*), which gives
the realistic fall-off of overtone intensities (roughly 1/10 - 1/30 per order).

A finite-difference tridiagonal eigensolver makes this a few milliseconds per
oscillator, i.e. no expensive ab initio calculation is needed.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.linalg import eigh_tridiagonal


@lru_cache(maxsize=256)
def morse_transitions(we: float, wexe: float, n_levels: int = 7,
                      n_grid: int = 2000, dipole_decay: float = 2.0):
    """Return (wavenumbers, relative intensities) for 0 -> v, v = 1..n_levels-1.

    Intensities are proportional to nu * |<0|mu|v>|^2 and normalised so that the
    fundamental (0 -> 1) has intensity 1.
    """
    B = wexe
    De = we ** 2 / (4.0 * wexe)
    y = np.linspace(-1.5, 8.0, n_grid)
    h = y[1] - y[0]
    V = De * (1.0 - np.exp(-y)) ** 2
    diag = 2.0 * B / h ** 2 + V
    off = -B / h ** 2 * np.ones(n_grid - 1)
    E, psi = eigh_tridiagonal(diag, off, select="i",
                              select_range=(0, n_levels - 1))
    psi = psi / np.sqrt(h)
    mu = y * np.exp(-y / dipole_decay)
    M = np.array([np.sum(psi[:, 0] * mu * psi[:, v]) * h
                  for v in range(1, n_levels)])
    nu = E[1:] - E[0]
    inten = nu * M ** 2
    return nu, inten / inten[0]


def morse_analytic(we: float, wexe: float, v: int) -> float:
    """Analytic 0 -> v transition wavenumber for a Morse oscillator."""
    return v * we - v * (v + 1) * wexe
