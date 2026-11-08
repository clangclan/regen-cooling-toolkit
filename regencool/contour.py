"""Hot-gas wall contour r(x) in mm. x = 0 at the injector face, increasing toward the nozzle exit."""
import csv
import math

import numpy as np
from scipy.interpolate import PchipInterpolator


class Contour:
    """Shape-preserving (PCHIP) interpolation of the hot-gas wall radius.

    rc_throat: throat wall radius of curvature [mm] used in the Bartz correlation. If it is not given it is
    estimated from a quadratic fit r = r_t + (x - x_t)^2 / (2 rc) over x_t +- r_t/2.
    """

    def __init__(self, x, r, rc_throat=None):
        x = np.asarray(x, float)
        r = np.asarray(r, float)
        o = np.argsort(x, kind="stable")
        x, r = x[o], r[o]
        x, iu = np.unique(x, return_index=True)
        r = r[iu]
        if len(x) < 4:
            raise ValueError("contour needs at least 4 points")
        self.x, self.r = x, r
        self.R = PchipInterpolator(x, r)
        self.dR = self.R.derivative()
        self.x0, self.x_end = float(x[0]), float(x[-1])
        i = int(np.argmin(r))
        xx = np.linspace(x[max(i - 1, 0)], x[min(i + 1, len(x) - 1)], 4001)
        k = int(np.argmin(self.R(xx)))
        self.x_t, self.r_t = float(xx[k]), float(self.R(xx[k]))
        if not (self.x0 < self.x_t < self.x_end):
            raise ValueError("the throat (minimum radius) must lie inside the contour")
        self.rc_throat = float(rc_throat) if rc_throat else self._throat_curvature()

    def _throat_curvature(self):
        h = 0.5 * self.r_t
        xx = np.linspace(max(self.x_t - h, self.x0), min(self.x_t + h, self.x_end), 101)
        a = np.polyfit(xx - self.x_t, self.R(xx), 2)[0]
        return 1.0 / (2.0 * a) if a > 1e-9 else 1.5 * self.r_t

    def area_ratio(self, x):
        """A/A_t (>= 1)."""
        return np.maximum((self.R(x) / self.r_t) ** 2, 1.0)

    @property
    def eps_exit(self):
        return float((self.R(self.x_end) / self.r_t) ** 2)

    @property
    def eps_max(self):
        return float(np.max((self.r / self.r_t) ** 2))

    @property
    def throat_area_m2(self):
        return math.pi * (self.r_t * 1e-3) ** 2

    def describe(self):
        return ("contour: x %.1f..%.1f mm, throat x=%.2f r=%.3f mm (rc %.2f mm), r_chamber %.2f mm, "
                "contraction %.2f, exit eps %.2f" % (self.x0, self.x_end, self.x_t, self.r_t, self.rc_throat,
                                                     float(self.R(self.x0)), (float(self.R(self.x0)) / self.r_t) ** 2,
                                                     self.eps_exit))

    # ------------------------------------------------------------------ constructors
    @classmethod
    def from_csv(cls, path, rc_throat=None, x_col="x_mm", r_col="r_mm"):
        """CSV with a header row; default columns x_mm, r_mm (hot-gas wall radius)."""
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        return cls([float(o[x_col]) for o in rows], [float(o[r_col]) for o in rows], rc_throat)

    @classmethod
    def parametric(cls, r_chamber, r_throat, l_cyl, conv_half_angle=30.0, r1_ratio=1.5, r2_ratio=1.0,
                   rn_ratio=0.382, eps_exit=4.0, divergent="bell", theta_n=25.0, theta_e=10.0,
                   bell_fraction=0.8, ds=0.1, rc_throat=None):
        """Cylinder + fillet R2 + cone + throat arc R1 | throat arc Rn + bell (quadratic Bezier) or cone.

        Radii R1, R2, Rn are given as ratios to r_throat. Bell length = bell_fraction x 15-deg cone length.
        """
        rt, rc = float(r_throat), float(r_chamber)
        tc = math.radians(conv_half_angle)
        R1, R2, Rn = r1_ratio * rt, r2_ratio * rt, rn_ratio * rt
        re = rt * math.sqrt(eps_exit)
        x2 = float(l_cyl)
        p2 = (x2 + R2 * math.sin(tc), rc - R2 * (1 - math.cos(tc)))       # end of chamber fillet
        r1y = rt + R1 * (1 - math.cos(tc))                                 # start of upstream throat arc
        if p2[1] <= r1y:
            raise ValueError("convergent section has no straight part: reduce r1_ratio/r2_ratio")
        x1 = p2[0] + (p2[1] - r1y) / math.tan(tc)
        xt = x1 + R1 * math.sin(tc)
        pts = []

        def n_of(length):
            return max(int(math.ceil(length / ds)), 2)

        def add(xs, rs):
            pts.extend(zip(np.atleast_1d(xs), np.atleast_1d(rs)))

        m = n_of(x2)
        add(np.linspace(0.0, x2, m), np.full(m, rc))
        a = np.linspace(0.0, tc, n_of(R2 * tc))
        add(x2 + R2 * np.sin(a), rc - R2 * (1 - np.cos(a)))
        s = np.linspace(0.0, 1.0, n_of(math.hypot(x1 - p2[0], p2[1] - r1y)))
        add(p2[0] + s * (x1 - p2[0]), p2[1] + s * (r1y - p2[1]))
        a = np.linspace(tc, 0.0, n_of(R1 * tc))
        add(xt - R1 * np.sin(a), rt + R1 * (1 - np.cos(a)))
        tn = math.radians(theta_n)
        a = np.linspace(0.0, tn, n_of(Rn * tn))
        add(xt + Rn * np.sin(a), rt + Rn * (1 - np.cos(a)))
        xn, rn_ = xt + Rn * math.sin(tn), rt + Rn * (1 - math.cos(tn))
        if re <= rn_:
            raise ValueError("exit radius is smaller than the end of the throat arc: increase eps_exit")
        if divergent == "cone":
            xe = xn + (re - rn_) / math.tan(tn)
            s = np.linspace(0.0, 1.0, n_of(xe - xn))
            add(xn + s * (xe - xn), rn_ + s * (re - rn_))
        elif divergent == "bell":
            te = math.radians(theta_e)
            t15 = math.radians(15.0)
            ln = bell_fraction * ((math.sqrt(eps_exit) - 1) * rt + Rn * (1 / math.cos(t15) - 1)) / math.tan(t15)
            xe = xt + ln
            m1, m2 = math.tan(tn), math.tan(te)
            qx = (re - rn_ + m1 * xn - m2 * xe) / (m1 - m2)
            qr = rn_ + m1 * (qx - xn)
            if not (xn < qx < xe):
                raise ValueError("bell: theta_n / theta_e / bell_fraction are inconsistent")
            s = np.linspace(0.0, 1.0, n_of(xe - xn))
            add((1 - s) ** 2 * xn + 2 * s * (1 - s) * qx + s * s * xe,
                (1 - s) ** 2 * rn_ + 2 * s * (1 - s) * qr + s * s * re)
        else:
            raise ValueError("divergent must be 'bell' or 'cone'")
        P = np.array(pts)
        return cls(P[:, 0], P[:, 1], rc_throat=rc_throat or 0.5 * (R1 + Rn))

    def to_csv(self, path, dx=0.25):
        xg = np.arange(self.x0, self.x_end + 1e-9, dx)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["x_mm", "r_mm"])
            for a, b in zip(xg, self.R(xg)):
                w.writerow(["%.4f" % a, "%.5f" % b])
