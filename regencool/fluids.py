"""Coolant property tables (single phase, one design pressure).

A coolant is a TabulatedFluid on a temperature grid: specific enthalpy H(T), density, viscosity, conductivity
and cp. The 1D solver marches in enthalpy, so compressed liquid -> pseudo-critical -> supercritical states are
handled continuously. The table is evaluated at a single pressure (the coolant inlet pressure by default), so the
coolant must stay above its critical pressure everywhere (there is no boiling model): compare the solver's P_min
with fluid.Pc.

Built-in: isopropanol (IPA). Any other coolant can be supplied as a CSV table (e.g. exported from the NIST
WebBook isobaric data) through TabulatedFluid.from_csv.
"""
import csv
import math

import numpy as np

R_U = 8.314462  # J/mol/K


class TabulatedFluid:
    def __init__(self, name, T, H, rho, mu, k, cp=None, Tc=None, Pc=None, omega=None, p_ref=None):
        T = np.asarray(T, float)
        H = np.asarray(H, float)
        if np.any(np.diff(T) <= 0) or np.any(np.diff(H) <= 0):
            raise ValueError("T and H(T) must increase monotonically")
        self.name, self.Tc, self.Pc, self.omega, self.p_ref = name, Tc, Pc, omega, p_ref
        self.tab = dict(T=T, H=H, rho=np.asarray(rho, float), mu=np.asarray(mu, float), k=np.asarray(k, float),
                        cp=np.gradient(H, T) if cp is None else np.asarray(cp, float))
        self.T, self.H = T, H
        self.T_pc = float(T[np.argmax(self.tab["cp"])])   # pseudo-critical temperature (cp peak)

    def at(self, T, key):
        return float(np.interp(T, self.T, self.tab[key]))

    def H_of_T(self, T):
        return float(np.interp(T, self.T, self.H))

    def T_of_H(self, h):
        return float(np.interp(h, self.H, self.T))

    def psat_bar(self, T):
        """Lee-Kesler vapour pressure [bar], clipped at Tc."""
        if None in (self.Tc, self.Pc, self.omega):
            raise ValueError("psat needs Tc, Pc and omega")
        Tr = np.minimum(np.asarray(T, float) / self.Tc, 0.99999)
        f0 = 5.92714 - 6.09648 / Tr - 1.28862 * np.log(Tr) + 0.169347 * Tr ** 6
        f1 = 15.2518 - 15.6875 / Tr - 13.4721 * np.log(Tr) + 0.43577 * Tr ** 6
        return self.Pc / 1e5 * np.exp(f0 + self.omega * f1)

    @classmethod
    def from_csv(cls, path, name="csv", **kw):
        """Columns: T_K, H_J_kg, rho_kg_m3, mu_Pa_s, k_W_mK [, cp_J_kgK]. Extra keywords: Tc, Pc [Pa], omega."""
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        col = lambda c: np.array([float(o[c]) for o in rows])  # noqa: E731
        cp = col("cp_J_kgK") if rows and "cp_J_kgK" in rows[0] else None
        return cls(name, col("T_K"), col("H_J_kg"), col("rho_kg_m3"), col("mu_Pa_s"), col("k_W_mK"), cp, **kw)


# ---------------------------------------------------------------------------------------------- isopropanol
IPA_TC, IPA_PC, IPA_OM, IPA_MW, IPA_VC = 508.3, 47.62e5, 0.665, 60.096, 220e-6


def _ipa_rho_liq(T):
    T = min(T, 507.0)
    return IPA_MW * 1.1800 / 0.26475 ** (1 + (1 - T / IPA_TC) ** 0.24653)


def _ipa_mu_liq(T):
    return math.exp(-8.23 + 2282.2 / T - 0.98495 * math.log(T))


_CPL_T = np.array([280, 290, 320, 350, 380, 410, 440, 460.0])
_CPL_V = np.array([2.47, 2.56, 2.83, 3.10, 3.36, 3.62, 3.92, 4.15]) * 1e3


def _ipa_cp_liq(T):
    return float(np.interp(T, _CPL_T, _CPL_V))


def _ipa_k_liq(T):
    return 0.2029 - 2.278e-4 * T


_IG_T = np.array([250, 298, 400, 500, 600, 700, 800.0])
_IG_V = np.array([78.0, 89.3, 112.6, 133.5, 151.4, 166.6, 179.5])   # ideal-gas cp, J/mol/K


def _ipa_cp_ig(T):
    return float(np.interp(T, _IG_T, _IG_V))


def _pr_state(T, P):
    """Peng-Robinson molar volume (dense root) and enthalpy departure for IPA."""
    kap = 0.379642 + 1.48503 * IPA_OM - 0.164423 * IPA_OM ** 2 + 0.016666 * IPA_OM ** 3
    ac = 0.45724 * R_U * R_U * IPA_TC * IPA_TC / IPA_PC
    bb = 0.07780 * R_U * IPA_TC / IPA_PC
    al = (1 + kap * (1 - math.sqrt(T / IPA_TC))) ** 2
    a = ac * al
    dadT = -ac * kap * math.sqrt(al) / math.sqrt(T * IPA_TC)
    A = a * P / (R_U * T) ** 2
    B = bb * P / (R_U * T)
    co = [1, -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)]
    rr = sorted(z.real for z in np.roots(co) if abs(z.imag) < 1e-9 and z.real > B)
    Z = rr[0]
    v = Z * R_U * T / P
    s2 = math.sqrt(2)
    Hdep = R_U * T * (Z - 1) + (T * dadT - a) / (2 * s2 * bb) * math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B))
    return v, Hdep


def isopropanol(p=52e5, T_max=800.0, dT=0.5):
    """IPA table at pressure p [Pa], 280 K .. T_max.

    Liquid (T <= 450 K): DIPPR-type correlations. Above: volume-translated Peng-Robinson enthalpy/density
    with Jossi-Stiel-Thodos (viscosity) and Stiel-Thodos (conductivity) dense-fluid transport, blended for
    continuity at 450 K and faded to the unscaled model by 560 K.
    Expected accuracy: liquid ~2-5 %, near/above the pseudo-critical point ~15-30 % (no critical enhancement of k).
    """
    T = np.arange(280.0, T_max + 1e-9, dT)
    mwk = IPA_MW * 1e-3
    v293, _ = _pr_state(293.15, p)
    csh = v293 - mwk / _ipa_rho_liq(293.15) * 1.002                   # volume translation at 293 K
    Hig = np.zeros_like(T)
    for i in range(1, len(T)):
        Hig[i] = Hig[i - 1] + 0.5 * (_ipa_cp_ig(T[i]) + _ipa_cp_ig(T[i - 1])) * (T[i] - T[i - 1])
    Hpr = np.zeros_like(T)
    Vpr = np.zeros_like(T)
    for i, t in enumerate(T):
        v, hd = _pr_state(t, p)
        Hpr[i] = (Hig[i] + hd) / mwk
        Vpr[i] = v - csh
    rho_pr = mwk / Vpr
    Ts = 450.0
    i_s = int(round((Ts - T[0]) / dT))
    H = np.zeros_like(T)
    for i in range(1, i_s + 1):
        H[i] = H[i - 1] + 0.5 * (_ipa_cp_liq(T[i]) + _ipa_cp_liq(T[i - 1])) * (T[i] - T[i - 1])
    # PR underestimates the liquid cp of alcohols: scale the PR increments, fading out above 560 K
    scl = _ipa_cp_liq(Ts) / ((Hpr[i_s + 2] - Hpr[i_s - 2]) / (4 * dT))
    for i in range(i_s + 1, len(T)):
        f = scl + (1 - scl) * min(max((T[i] - Ts) / 110.0, 0), 1)
        H[i] = H[i - 1] + f * (Hpr[i] - Hpr[i - 1])
    rho = np.array([_ipa_rho_liq(t) if t <= Ts else 0.0 for t in T])
    rs = _ipa_rho_liq(Ts) / rho_pr[i_s]
    for i in range(i_s + 1, len(T)):
        f = rs + (1 - rs) * min(max((T[i] - Ts) / 110.0, 0), 1)
        rho[i] = rho_pr[i] * f
    rhoc = mwk / IPA_VC

    def mu0_chung(t):
        ts = 1.2593 * t / IPA_TC
        om = 1.16145 * ts ** -0.14874 + 0.52487 * math.exp(-0.7732 * ts) + 2.16178 * math.exp(-2.43787 * ts)
        fc = 1 - 0.2756 * IPA_OM + 0.059035 * 0.652 ** 4 + 0.143
        return 40.785 * fc * math.sqrt(IPA_MW * t) / ((IPA_VC * 1e6) ** (2 / 3) * om) * 1e-7

    xi = IPA_TC ** (1 / 6) / (IPA_MW ** 0.5 * (IPA_PC / 101325) ** (2 / 3))

    def mu_jst(t, r):
        rr = r / rhoc
        q = 0.1023 + 0.023364 * rr + 0.058533 * rr ** 2 - 0.040758 * rr ** 3 + 0.0093324 * rr ** 4
        return mu0_chung(t) + max(q ** 4 - 1e-4, 0) / xi * 1e-3

    zc = IPA_PC * IPA_VC / (R_U * IPA_TC)
    gam = 210 * (IPA_TC * IPA_MW ** 3 / (IPA_PC / 1e5) ** 4) ** (1 / 6)

    def k_st(t, r):
        rr = r / rhoc
        if rr < 0.5:
            d = 1.22e-2 * (math.exp(0.535 * rr) - 1)
        elif rr < 2.0:
            d = 1.14e-2 * (math.exp(0.67 * rr) - 1.069)
        else:
            d = 2.60e-3 * (math.exp(1.155 * rr) + 2.016)
        return 0.0175 * (t / 400) ** 1.6 + d / (gam * zc ** 5)

    mu = np.zeros_like(T)
    kk = np.zeros_like(T)
    sm = _ipa_mu_liq(Ts) / mu_jst(Ts, rho[i_s])
    sk = _ipa_k_liq(Ts) / k_st(Ts, rho[i_s])
    for i, t in enumerate(T):
        if t <= Ts:
            mu[i], kk[i] = _ipa_mu_liq(t), _ipa_k_liq(t)
        else:
            f = min(max((t - Ts) / 110.0, 0), 1)
            mu[i] = mu_jst(t, rho[i]) * (sm + (1 - sm) * f)
            kk[i] = k_st(t, rho[i]) * (sk + (1 - sk) * f)
    return TabulatedFluid("isopropanol", T, H, rho, mu, kk, Tc=IPA_TC, Pc=IPA_PC, omega=IPA_OM, p_ref=p)


def make(cfg, p_default):
    """Coolant from the [coolant] config table."""
    name = cfg.get("fluid", "isopropanol").lower()
    p = cfg.get("property_pressure_bar", p_default / 1e5) * 1e5
    if name in ("ipa", "isopropanol", "2-propanol"):
        return isopropanol(p)
    if name == "csv":
        kw = {k: cfg[k] for k in ("Tc", "omega") if k in cfg}
        if "Pc_bar" in cfg:
            kw["Pc"] = cfg["Pc_bar"] * 1e5
        return TabulatedFluid.from_csv(cfg["table"], name=cfg.get("name", "csv"), p_ref=p, **kw)
    raise ValueError("unknown coolant fluid %r (use 'isopropanol' or 'csv')" % name)
