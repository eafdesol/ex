"""Minimal PLS1 (NIPALS) giving coefficients for every number of LVs at once."""
from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter


def preprocess(X, window=15, poly=2, deriv=1):
    """Savitzky-Golay derivative (removes offset; reduces slope/scatter)."""
    return savgol_filter(X, window, poly, deriv=deriv, axis=1)


class PLS1:
    def __init__(self, n_max=15):
        self.n_max = n_max

    def fit(self, X, y):
        self.xm, self.ym = X.mean(0), y.mean()
        E, f = X - self.xm, y - self.ym
        A = min(self.n_max, X.shape[0] - 1, X.shape[1])
        W, P, q, T = [], [], [], []
        for _ in range(A):
            w = E.T @ f
            nw = np.linalg.norm(w)
            if nw < 1e-12:
                break
            w /= nw
            t = E @ w
            tt = t @ t
            p = E.T @ t / tt
            qa = f @ t / tt
            E = E - np.outer(t, p)
            f = f - qa * t
            W.append(w); P.append(p); q.append(qa); T.append(t)
        self.W, self.P = np.array(W).T, np.array(P).T
        self.q, self.T = np.array(q), np.array(T).T
        R = self.W @ np.linalg.inv(self.P.T @ self.W)
        self.B = np.cumsum(R * self.q, axis=1)     # p x A, column a = a+1 LVs
        return self

    def predict_all(self, X):
        return (X - self.xm) @ self.B + self.ym     # n x A

    def predict(self, X, n_lv):
        return self.predict_all(X)[:, n_lv - 1]

    def vip(self, n_lv):
        W = self.W[:, :n_lv]
        ss = (self.q[:n_lv] ** 2) * np.sum(self.T[:, :n_lv] ** 2, axis=0)
        return np.sqrt(W.shape[0] * (W ** 2 @ ss) / ss.sum())


def cv_rmse(X, y, n_max=15, folds=5, seed=0):
    """RMSECV for 1..A latent variables (A may be < n_max)."""
    n = len(y)
    idx = np.random.default_rng(seed).permutation(n)
    pred = None
    for k in range(folds):
        te = idx[k::folds]
        tr = np.setdiff1d(idx, te)
        P = PLS1(n_max).fit(X[tr], y[tr]).predict_all(X[te])
        if pred is None:
            pred = np.full((n, P.shape[1]), np.nan)
        A = min(P.shape[1], pred.shape[1])
        pred = pred[:, :A]
        pred[te] = P[:, :A]
    return np.sqrt(np.nanmean((pred - y[:, None]) ** 2, axis=0))


def choose_lv(rmse, tol=1.02):
    """Smallest LV count within 2 % of the minimum RMSECV (parsimony)."""
    return int(np.argmax(rmse <= tol * rmse.min()) + 1)
