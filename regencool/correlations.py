"""Coolant-side friction and heat-transfer correlations."""
import math


def f_rect_lam(Re, a):
    """Darcy friction factor, fully developed laminar flow in a rectangular duct (Shah & London), a = aspect <= 1."""
    fRe = 96 * (1 - 1.3553 * a + 1.9467 * a ** 2 - 1.7012 * a ** 3 + 0.9564 * a ** 4 - 0.2537 * a ** 5)
    return fRe / Re


def nu_rect_lam(a):
    """Nusselt number, laminar, constant wall temperature-type, rectangular duct (Shah & London)."""
    return 8.235 * (1 - 2.0421 * a + 3.0853 * a ** 2 - 2.4765 * a ** 3 + 1.0578 * a ** 4 - 0.1861 * a ** 5)


def f_colebrook(Re, rr):
    """Darcy friction factor, Colebrook-White (Haaland start), rr = ks/Dh."""
    f = (-1.8 * math.log10((rr / 3.7) ** 1.11 + 6.9 / Re)) ** -2
    for _ in range(30):
        f = (-2 * math.log10(rr / 3.7 + 2.51 / (Re * math.sqrt(f)))) ** -2
    return f


def f_smooth(Re):
    """Petukhov smooth-tube friction factor (used inside Gnielinski)."""
    return (0.79 * math.log(Re) - 1.64) ** -2


def curv_factor(Re, Dh, kappa):
    """Dean-type curvature factor on heat transfer: (Re (Dh/2R)^2)^(+-0.05).

    kappa [1/m] is the normal curvature of the channel path; kappa > 0 means the hot wall is the outer wall
    of the bend (enhancement), kappa < 0 the inner wall (reduction).
    """
    if abs(kappa) < 1e-9:
        return 1.0
    R = 1 / abs(kappa)
    D = Re * (Dh / (2 * R)) ** 2
    if D <= 1:
        return 1.0
    return D ** (0.05 if kappa > 0 else -0.05)


def coolant_nu_f(fluid, Re, Pr, Dh, a, ks, L_from_inlet, Tb, Tw, mub, muw):
    """Darcy f, Nu and regime ('lam' / 'turb') for one channel station.

    laminar : max(fully developed rectangular duct, Sieder-Tate thermal entry with L from the coolant inlet)
    turbulent: Gnielinski x roughness enhancement sqrt(f_rough/f_smooth) (capped at 2) x property-ratio
               correction: subcooled liquid -> (mu_b/mu_w)^0.11 (capped 1.35), near/above the pseudo-critical
               temperature -> Jackson (rho_w/rho_b)^0.3 (cp_avg/cp_b)^0.4 (clipped 0.5..1.5)
    transition 2300 < Re < 4000: linear blend.
    """
    a = min(a, 1 / a)
    nul = nu_rect_lam(a)
    gz = Re * Pr * Dh / max(L_from_inlet, 1e-4)
    nud = 1.86 * gz ** (1 / 3) * (mub / muw) ** 0.14
    nu_lam = max(nul, nud)
    f_lam = f_rect_lam(Re, a)
    if Re < 2300:
        return f_lam, nu_lam, "lam"
    fs = f_smooth(Re)
    fr = f_colebrook(Re, ks / Dh)
    nus = (fs / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(fs / 8) * (Pr ** (2 / 3) - 1))
    E = min((fr / fs) ** 0.5, 2.0)
    if Tb < 0.9 * fluid.T_pc:
        muw_l = fluid.at(min(Tw, 0.97 * fluid.T_pc), "mu")
        corr = min((mub / muw_l) ** 0.11, 1.35)
    else:
        rb, rw = fluid.at(Tb, "rho"), fluid.at(Tw, "rho")
        cpb = fluid.at(Tb, "cp")
        cpa = (fluid.H_of_T(Tw) - fluid.H_of_T(Tb)) / max(Tw - Tb, 1e-3)
        corr = min(max((rw / rb) ** 0.3 * (cpa / cpb) ** 0.4, 0.5), 1.5)
    nu = nus * E * corr
    if Re < 4000:
        w = (Re - 2300) / 1700
        nu = (1 - w) * nu_lam + w * nu
        fr = (1 - w) * f_rect_lam(2300, a) + w * fr
    return fr, max(nu, nul), "turb"
