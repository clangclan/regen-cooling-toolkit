"""Hot-gas side: Bartz coefficient, Mach number, adiabatic-wall temperature and static pressure along x.

Two models share one interface (on_grid(x) -> dict of arrays h0, M, Taw, g1, p, T):
  IdealGas   -- constant gamma / MW (quick estimates, hand checks)
  CanteraGas -- chamber HP equilibrium + isentropic expansion tables from cantera_gas.compute()
"""
import math

import numpy as np

R_U = 8314.46  # J/kmol/K


def bartz_h0(contour, x, mu, cp, Pr, pc, cstar):
    """Bartz (1957) coefficient without the boundary-layer property factor sigma [W/m^2/K].

    h0 = 0.026/Dt^0.2 * mu^0.2 cp / Pr^0.6 * (pc/c*)^0.8 * (Dt/rc)^0.1 * (At/A)^0.9
    """
    Dt = 2 * contour.r_t * 1e-3
    rc = contour.rc_throat * 1e-3
    ar = contour.area_ratio(x)
    return 0.026 / Dt ** 0.2 * (mu ** 0.2 * cp / Pr ** 0.6) * (pc / cstar) ** 0.8 * (Dt / rc) ** 0.1 * (1 / ar) ** 0.9


def bartz_sigma(Twg, T0, g1, omega=0.6):
    """Bartz property-variation factor; g1 = T0/T_static."""
    return 1.0 / ((0.5 * Twg / T0 * g1 + 0.5) ** (0.8 - 0.2 * omega) * g1 ** (0.2 * omega))


def mach_from_area(ar, gamma, supersonic):
    """Isentropic Mach number from A/A* (vectorized bisection)."""
    ar = np.atleast_1d(np.asarray(ar, float))
    sup = np.broadcast_to(np.asarray(supersonic, bool), ar.shape)
    g = gamma

    def f(M):
        return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1))) - ar

    lo = np.where(sup, 1.0, 1e-4)
    hi = np.where(sup, 8.0, 1.0)
    for _ in range(80):
        m = 0.5 * (lo + hi)
        go_up = (f(m) > 0) == (~sup)
        lo = np.where(go_up, m, lo)
        hi = np.where(go_up, hi, m)
    return 0.5 * (lo + hi)


def cstar_ideal(gamma, MW, T0):
    R = R_U / MW
    return math.sqrt(gamma * R * T0) / (gamma * math.sqrt((2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))))


class HotGas:
    T0 = None
    name = ""

    def sigma(self, Twg, g1):
        return bartz_sigma(Twg, self.T0, g1)

    def on_grid(self, x):
        raise NotImplementedError

    def describe(self):
        return "%s: T0 %.0f K, c*(Bartz) %.0f m/s, cp %.0f J/kg/K, mu %.3e Pa s, Pr %.3f" % (
            self.name, self.T0, self.cstar, self.cp, self.mu, self.Pr)


class IdealGas(HotGas):
    """Constant-gamma gas. mu defaults to Bartz' estimate 1.184e-7 MW^0.5 T0^0.6 [Pa s],
    Pr to the Eucken estimate 4 gamma/(9 gamma - 5). The Bartz mass flux uses pc/c* with
    c* = cstar if given, else eta_cstar x ideal c*."""

    def __init__(self, contour, pc, T0, gamma, MW, Pr=None, mu=None, cstar=None, eta_cstar=1.0):
        self.contour, self.pc, self.T0, self.gamma, self.MW = contour, pc, T0, gamma, MW
        self.R = R_U / MW
        self.cp = gamma * self.R / (gamma - 1)
        self.mu = mu if mu else 1.184e-7 * MW ** 0.5 * T0 ** 0.6
        self.Pr = Pr if Pr else 4 * gamma / (9 * gamma - 5)
        self.cstar_ideal = cstar_ideal(gamma, MW, T0)
        self.cstar = cstar if cstar else eta_cstar * self.cstar_ideal
        self.name = "ideal gas (gamma %.3f, MW %.2f)" % (gamma, MW)

    def on_grid(self, x):
        c, g = self.contour, self.gamma
        x = np.atleast_1d(np.asarray(x, float))
        ar = np.maximum(c.area_ratio(x), 1.000001)
        M = mach_from_area(ar, g, x > c.x_t)
        M = np.where(np.abs(x - c.x_t) < 1e-9, 1.0, M)
        g1 = 1 + (g - 1) / 2 * M * M
        rec = self.Pr ** (1 / 3)
        Taw = self.T0 * (1 + rec * (g - 1) / 2 * M * M) / g1
        p = self.pc * g1 ** (-g / (g - 1))
        return dict(h0=bartz_h0(c, x, self.mu, self.cp, self.Pr, self.pc, self.cstar), M=M, Taw=Taw, g1=g1, p=p,
                    T=self.T0 / g1, gam=np.full_like(x, g))


class CanteraGas(HotGas):
    """Tabulated equilibrium ('eq') or frozen ('fr') expansion. Bartz uses the chamber frozen cp, mu, Pr and
    pc/c* with c* = cstar if given, else eta_cstar x c*_mode. T0 is the ideal (100 %) equilibrium value,
    which is conservative for the wall heat flux."""

    KEYS = ("p", "T", "M", "Taw", "gam_s", "gam_f", "cp_f", "u")

    def __init__(self, contour, table, pc, mode="eq", cstar=None, eta_cstar=1.0):
        if mode not in ("eq", "fr"):
            raise ValueError("mode must be 'eq' or 'fr'")
        ch = table["chamber"]
        self.contour, self.table, self.pc, self.mode = contour, table, pc, mode
        self.T0, self.mu, self.cp, self.Pr = ch["T"], ch["mu"], ch["cp_f"], ch["Pr_f"]
        self.cstar_ideal = table["cstar_" + mode]
        self.cstar = cstar if cstar else eta_cstar * self.cstar_ideal
        self.name = "Cantera %s expansion" % ("equilibrium" if mode == "eq" else "frozen")
        rows = table["exp_" + mode]
        self.branch = {}
        for sup in (False, True):
            rr = [o for o in rows if o["sup"] == sup or abs(o["eps"] - 1) < 1e-9]
            e = np.array([o["eps"] for o in rr])
            o_ = np.argsort(e)
            self.branch[sup] = {k: np.array([o[k] for o in rr])[o_] for k in ("eps",) + self.KEYS}
        if self.branch[True]["eps"].max() < contour.eps_max * 0.999:
            raise ValueError("expansion table ends at eps %.2f < contour eps %.2f: recompute the gas table"
                             % (self.branch[True]["eps"].max(), contour.eps_max))

    def local(self, x):
        ar = float(self.contour.area_ratio(x))
        b = self.branch[bool(x > self.contour.x_t)]
        return {k: float(np.interp(ar, b["eps"], b[k])) for k in self.KEYS}

    def on_grid(self, x):
        x = np.atleast_1d(np.asarray(x, float))
        L = [self.local(v) for v in x]
        arr = {k: np.array([o[k] for o in L]) for k in self.KEYS}
        return dict(h0=bartz_h0(self.contour, x, self.mu, self.cp, self.Pr, self.pc, self.cstar), M=arr["M"],
                    Taw=arr["Taw"], g1=self.T0 / arr["T"], p=arr["p"], T=arr["T"], gam=arr["gam_s"],
                    gam_f=arr["gam_f"])
