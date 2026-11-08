"""Channel-size distribution optimization.

Variables: log s at a few axial knots (PCHIP through the knots gives s(x)). Objective (Nelder-Mead):
    J = max Twg  + 300 max(0, dP - dP_max) + 300 max(0, P_min - P_exit)          [deg C, bar]
robust variant (d > 0): J = max(J_nominal, max Twg at +d) + 300 max(0, P_min - P_exit at -d)
    (+d: features oversize -> thinner wall, lower hc; -d: undersize -> higher dP)
plus a penalty when any knot falls below s_min + d. Hot-wall thickness t(x), channel count and helix stay fixed.
"""
import math

import numpy as np
from scipy.optimize import minimize

from .channel import knots_fun
from .solver1d import ChannelClosed


def default_knots(contour):
    x0, xt, xe = contour.x0, contour.x_t, contour.x_end
    return [x0, x0 + 0.4 * (xt - x0), x0 + 0.8 * (xt - x0), xt, xt + 0.35 * (xe - xt), xe]


class Optimizer:
    def __init__(self, case, knots_x=None, dP_max=None, P_min=None, robust_d=0.0, s_min=0.3, C_hg=None,
                 log=print):
        self.case, self.log = case, log
        self.kx = np.array(knots_x or default_knots(case.contour), float)
        self.dP_max = dP_max if dP_max is not None else (case.dP_max or 10.0)
        self.P_min = P_min if P_min is not None else (case.P_min or 0.0)
        self.d, self.s_min, self.C_hg = robust_d, s_min, C_hg
        self.nfev = 0
        self.best = (math.inf, None)

    def s_fun(self, logs):
        return knots_fun(np.c_[self.kx, np.exp(logs)])

    def evaluate(self, logs, d=0.0):
        D = self.case.design(d, s_fun=self.s_fun(logs))
        return self.case.metrics(self.case.solve(D, C_hg=self.C_hg))

    def objective(self, logs):
        self.nfev += 1
        try:
            m0 = self.evaluate(logs)
            J = m0["Twg"] + 300 * max(0, m0["dP"] - self.dP_max) + 300 * max(0, self.P_min - m0["Pmin"])
            if self.d > 0:
                mp = self.evaluate(logs, +self.d)
                mm = self.evaluate(logs, -self.d)
                J = max(J, mp["Twg"]) + 300 * max(0, self.P_min - mm["Pmin"])
            smin = float(np.exp(logs).min())
            if smin < self.s_min + self.d:
                J += 3000 * (self.s_min + self.d - smin)
        except (ChannelClosed, ValueError, ZeroDivisionError, FloatingPointError):
            J = 1e4
        if J < self.best[0]:
            self.best = (J, np.array(logs))
            if self.nfev > 1:
                self.log("  eval %4d  J %.1f  s = %s" % (self.nfev, J, np.round(np.exp(logs), 3).tolist()))
        return J

    def start(self, target=0.85):
        """Current s(x) at the knots, scaled uniformly so that dP ~ target x dP_max."""
        s0 = np.asarray(self.case.design().s(self.kx), float) + self.d
        lo, hi = 0.3, 4.0
        for _ in range(14):
            k = math.sqrt(lo * hi)
            try:
                dp = self.evaluate(np.log(s0 * k))["dP"]
            except (ChannelClosed, ValueError):
                dp = 1e9
            if dp > target * self.dP_max:
                lo = k
            else:
                hi = k
        return np.log(s0 * hi)

    def run(self, maxfev=200, x0=None):
        x0 = self.start() if x0 is None else np.asarray(x0, float)
        self.log("  start s = %s" % np.round(np.exp(x0), 3).tolist())
        r = minimize(self.objective, x0, method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=0.01, fatol=0.5))
        logs = self.best[1] if self.best[1] is not None else r.x
        s = np.exp(logs)
        out = dict(knots=np.c_[self.kx, s].round(4).tolist(), J=float(self.best[0]), nfev=self.nfev,
                   nominal=self.evaluate(logs))
        if self.d > 0:
            out["plus_d"] = self.evaluate(logs, +self.d)
            out["minus_d"] = self.evaluate(logs, -self.d)
        return out


def toml_snippet(knots):
    return "s_mm = [%s]\n" % ", ".join("[%.2f, %.4f]" % (a, b) for a, b in knots)
