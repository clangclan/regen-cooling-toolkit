"""Injector orifice sizing (first cut): mass flux G through a sharp-edged orifice, A = mdot / (Cd G).

Models
  spi        single-phase incompressible   G = sqrt(2 rho (P1 - P2))
  ipa_real   isentropic real-fluid expansion of isopropanol (Peng-Robinson, stable root by minimum Gibbs energy):
             G = max over P2' in [P2, P1] of rho(P2', s1) sqrt(2 (h1 - h(P2', s1)))  -> choking of a hot
             (e.g. regeneratively heated, supercritical) fuel is captured. Volume-translated density.
  n2o_dyer   self-pressurising nitrous oxide, Dyer et al. (2007) non-homogeneous non-equilibrium model
             G = k/(1+k) G_SPI + 1/(1+k) G_HEM,  k = sqrt((P1 - P2)/(Pv - P2)),
             HEM by isentropic flashing along the saturation line; N2O saturation from ESDU 91022 correlations.
             Below Pv <= P2 (no flashing) it reduces to SPI.

Uncertainty: Cd is the main unknown (0.6-0.85 for sharp-edged / short-tube orifices, depends on L/D, inlet
rounding, cavitation); calibrate with cold-flow tests. Dyer's model is an engineering fit (+-20-30 % on G).
"""
import math

import numpy as np
from scipy.optimize import brentq

from . import fluids

R_U = fluids.R_U


def spi(rho, P1, P2):
    return math.sqrt(2 * rho * max(P1 - P2, 0.0))


def diameter(mdot, G, Cd, n=1):
    """Orifice diameter [mm] for n equal holes."""
    A = mdot / (n * Cd * G)
    return math.sqrt(4 * A / math.pi) * 1e3


# ---------------------------------------------------------------------------------------------- IPA real fluid
class _IPA:
    """Peng-Robinson IPA with ideal-gas cp from fluids; h, s relative to 250 K / 1 bar ideal gas."""

    def __init__(self):
        om, Tc, Pc = fluids.IPA_OM, fluids.IPA_TC, fluids.IPA_PC
        self.kap = 0.379642 + 1.48503 * om - 0.164423 * om ** 2 + 0.016666 * om ** 3
        self.ac = 0.45724 * R_U * R_U * Tc * Tc / Pc
        self.b = 0.07780 * R_U * Tc / Pc
        self.mw = fluids.IPA_MW * 1e-3
        Tg = np.arange(250.0, 1200.0, 0.5)
        cp = np.array([fluids._ipa_cp_ig(t) for t in Tg])
        self.Tg = Tg
        self.Hig = np.concatenate([[0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(Tg))])
        self.Sig = np.concatenate([[0], np.cumsum(0.5 * (cp[1:] / Tg[1:] + cp[:-1] / Tg[:-1]) * np.diff(Tg))])
        v293 = self._roots(293.15, 50e5)[0][0] * R_U * 293.15 / 50e5
        self.csh = v293 - self.mw / fluids._ipa_rho_liq(293.15)          # volume translation (density only)
        self._ts = {}

    def _roots(self, T, P):
        """All PR roots at (T, P): list of (Z, ln phi, H_dep [J/mol], S_dep [J/mol/K]), ascending Z."""
        Tc, s2 = fluids.IPA_TC, math.sqrt(2)
        al = (1 + self.kap * (1 - math.sqrt(T / Tc))) ** 2
        a = self.ac * al
        dadT = -self.ac * self.kap * math.sqrt(al) / math.sqrt(T * Tc)
        A = a * P / (R_U * T) ** 2
        B = self.b * P / (R_U * T)
        Zs = sorted(z.real for z in np.roots([1, -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)])
                    if abs(z.imag) < 1e-9 and z.real > B)
        out = []
        for Z in Zs:
            L = math.log((Z + (1 + s2) * B) / (Z + (1 - s2) * B))
            lnphi = (Z - 1) - math.log(Z - B) - A / (2 * s2 * B) * L
            Hd = R_U * T * (Z - 1) + (T * dadT - a) / (2 * s2 * self.b) * L
            Sd = R_U * math.log(Z - B) + dadT / (2 * s2 * self.b) * L
            out.append((Z, lnphi, Hd, Sd))
        return out

    def hs(self, T, P, root="stable"):
        """h [J/kg], s [J/kg/K], rho [kg/m3]. root: 'stable' (minimum Gibbs energy), 'liq' or 'vap'."""
        rr = self._roots(T, P)
        o = rr[0] if root == "liq" else rr[-1] if root == "vap" else min(rr, key=lambda e: e[1])
        Z, _, Hd, Sd = o
        v = Z * R_U * T / P
        h = (np.interp(T, self.Tg, self.Hig) + Hd) / self.mw
        s = (np.interp(T, self.Tg, self.Sig) - R_U * math.log(P / 1e5) + Sd) / self.mw
        return h, s, self.mw / (v - self.csh)

    def tsat(self, P):
        """PR saturation temperature (equal fugacities) at P < Pc, or None."""
        if P >= 0.999 * fluids.IPA_PC:
            return None
        if P not in self._ts:
            self._ts[P] = self._tsat(P)
        return self._ts[P]

    def _tsat(self, P):
        def f(T):
            rr = self._roots(T, P)
            return rr[0][1] - rr[-1][1] if len(rr) > 1 else None

        Ts = [t for t in np.arange(260.0, fluids.IPA_TC, 0.5)]
        prev = None
        for t in Ts:
            v = f(t)
            if v is None:
                prev = None
                continue
            if prev is not None and prev[1] * v <= 0:
                return brentq(lambda tt: f(tt) if f(tt) is not None else 0.0, prev[0], t)
            prev = (t, v)
        return None


_ipa = None


def _ipa_isentrope(s1, p, T_hi):
    """State on the isentrope s1 at pressure p: (T, h, rho, quality or None). HEM mixture inside the dome."""
    Ts = _ipa.tsat(p)
    if Ts is None:
        T2 = brentq(lambda t: _ipa.hs(t, p)[1] - s1, 255.0, T_hi)
        h2, _, r2 = _ipa.hs(T2, p)
        return T2, h2, r2, None
    hl, sl, rl = _ipa.hs(Ts, p, "liq")
    hv, sv, rv = _ipa.hs(Ts, p, "vap")
    if s1 <= sl:
        T2 = brentq(lambda t: _ipa.hs(t, p, "liq")[1] - s1, 255.0, Ts)
        h2, _, r2 = _ipa.hs(T2, p, "liq")
        return T2, h2, r2, None
    if s1 >= sv:
        T2 = brentq(lambda t: _ipa.hs(t, p, "vap")[1] - s1, Ts, max(T_hi, Ts + 1.0))
        h2, _, r2 = _ipa.hs(T2, p, "vap")
        return T2, h2, r2, None
    x = (s1 - sl) / (sv - sl)
    return Ts, hl + x * (hv - hl), 1 / (x / rv + (1 - x) / rl), x


def ipa_real(T1, P1, P2, n=120):
    """Isentropic real-fluid orifice mass flux [kg/m2/s] of IPA from (T1, P1) into P2 (homogeneous equilibrium
    mixture where the isentrope enters the two-phase dome). Returns (G, info)."""
    global _ipa
    if _ipa is None:
        _ipa = _IPA()
    h1, s1, r1 = _ipa.hs(T1, P1)
    best = (0.0, None)
    for p in np.linspace(P1, P2, n)[1:]:
        T2, h2, r2, x = _ipa_isentrope(s1, p, T1 + 1.0)
        G = r2 * math.sqrt(max(2 * (h1 - h2), 0.0))
        if G > best[0]:
            best = (G, dict(P_throat_bar=p / 1e5, T_throat=T2, rho_throat=r2, u_throat=G / r2, quality=x))
    info = dict(best[1], rho1=r1, choked=bool(best[1]["P_throat_bar"] > P2 / 1e5 + 0.05),
                G_spi=spi(r1, P1, P2))
    return best[0], info


# ---------------------------------------------------------------------------------------------- N2O
N2O_TC, N2O_PC, N2O_RHOC = 309.57, 72.51e5, 452.0


def n2o_psat(T):
    Tr = T / N2O_TC
    th = 1 - Tr
    return N2O_PC * math.exp((-6.71893 * th + 1.35966 * th ** 1.5 - 1.3779 * th ** 2.5 - 4.051 * th ** 5) / Tr)


def n2o_rhol(T):
    th = 1 - T / N2O_TC
    return N2O_RHOC * math.exp(1.72328 * th ** (1 / 3) - 0.83950 * th ** (2 / 3) + 0.51060 * th - 0.10412 * th ** (4 / 3))


def n2o_rhov(T):
    u = N2O_TC / T - 1
    return N2O_RHOC * math.exp(-1.00900 * u ** (1 / 3) - 6.28792 * u ** (2 / 3) + 7.50332 * u - 7.90463 * u ** (4 / 3)
                               + 0.629427 * u ** (5 / 3))


def n2o_hl(T):
    th = 1 - T / N2O_TC
    return 1e3 * (-200 + 116.043 * th ** (1 / 3) - 917.225 * th ** (2 / 3) + 794.779 * th - 589.587 * th ** (4 / 3))


def n2o_hv(T):
    th = 1 - T / N2O_TC
    return 1e3 * (-200 + 440.055 * th ** (1 / 3) - 459.701 * th ** (2 / 3) + 434.081 * th - 485.338 * th ** (4 / 3))


_TN = np.arange(200.0, 309.0, 0.05)
_SL = np.concatenate([[0], np.cumsum([
    (n2o_hl(_TN[i + 1]) - n2o_hl(_TN[i]) - (n2o_psat(_TN[i + 1]) - n2o_psat(_TN[i]))
     / (0.5 * (n2o_rhol(_TN[i]) + n2o_rhol(_TN[i + 1])))) / (0.5 * (_TN[i] + _TN[i + 1]))
    for i in range(len(_TN) - 1)])])                      # saturated-liquid entropy: T ds = dh - v dP


def _sl(T):
    return float(np.interp(T, _TN, _SL))


def n2o_tsat(P):
    return brentq(lambda t: n2o_psat(t) - P, 200.0, N2O_TC - 0.01)


def n2o_dyer(T1, P1, P2, n=200):
    """Dyer NHNE mass flux [kg/m2/s] of liquid N2O at (T1, P1) into P2. Returns (G, info)."""
    Pv = n2o_psat(T1)
    if P1 < 0.99 * Pv:
        raise ValueError("N2O upstream pressure %.2f bar is below its vapour pressure %.2f bar at %.1f K"
                         % (P1 / 1e5, Pv / 1e5, T1))
    rl = n2o_rhol(T1)
    G_spi = spi(rl, P1, P2)
    info = dict(Pv_bar=Pv / 1e5, rho_l=rl, G_spi=G_spi, subcooling_K=(n2o_tsat(P1) - T1) if P1 < N2O_PC else None)
    if Pv <= P2:
        return G_spi, dict(info, G_hem=None, k=None)
    h1 = n2o_hl(T1) + (P1 - Pv) / rl
    s1 = _sl(T1)
    G_hem, arg = 0.0, None
    for p in np.linspace(Pv, P2, n)[1:]:
        t2 = n2o_tsat(p)
        sfg = (n2o_hv(t2) - n2o_hl(t2)) / t2
        x = (s1 - _sl(t2)) / sfg
        h2 = n2o_hl(t2) + x * (n2o_hv(t2) - n2o_hl(t2))
        r2 = 1 / (x / n2o_rhov(t2) + (1 - x) / n2o_rhol(t2))
        G = r2 * math.sqrt(max(2 * (h1 - h2), 0.0))
        if G > G_hem:
            G_hem, arg = G, dict(P_throat_bar=p / 1e5, quality=x, rho_throat=r2)
    k = math.sqrt((P1 - P2) / (Pv - P2))
    return k / (1 + k) * G_spi + 1 / (1 + k) * G_hem, dict(info, G_hem=G_hem, hem=arg, k=k)


# ---------------------------------------------------------------------------------------------- config driver
def size(entry, P2, fluid=None, regen=None):
    """One [injector.<name>] table -> sizing dict.

    entry keys: model ('spi' | 'ipa_real' | 'n2o_dyer'), mdot_kg_s, P_up_bar, T_K, Cd, n_holes, rho_kg_m3 (spi),
    from_regen (take T and P upstream from the 1D coolant outlet, minus dP_feed_bar)."""
    model = entry.get("model", "spi")
    T1, P1 = entry.get("T_K"), entry.get("P_up_bar")
    if entry.get("from_regen") and regen is not None:
        T1 = regen["Tout"]
        P1 = regen["P_out_bar"] - entry.get("dP_feed_bar", 0.0)
    if T1 is None or P1 is None:
        raise ValueError("injector entry needs T_K and P_up_bar (or from_regen = true)")
    P1 *= 1e5
    if P1 <= P2:
        raise ValueError("upstream pressure %.2f bar <= chamber %.2f bar" % (P1 / 1e5, P2 / 1e5))
    if model == "spi":
        rho = entry.get("rho_kg_m3") or (fluid.at(T1, "rho") if fluid is not None else None)
        if rho is None:
            raise ValueError("spi needs rho_kg_m3")
        G, info = spi(rho, P1, P2), dict(rho1=rho)
    elif model == "ipa_real":
        G, info = ipa_real(T1, P1, P2)
    elif model == "n2o_dyer":
        G, info = n2o_dyer(T1, P1, P2)
    else:
        raise ValueError("unknown injector model %r" % model)
    Cd, nh = entry.get("Cd", 0.75), int(entry.get("n_holes", 1))
    mdot = entry["mdot_kg_s"]
    return dict(model=model, T1_K=T1, P1_bar=P1 / 1e5, P2_bar=P2 / 1e5, dP_bar=(P1 - P2) / 1e5, G=G, Cd=Cd,
                n_holes=nh, mdot=mdot, d_mm=diameter(mdot, G, Cd, nh), info=info)
