"""2D thermo-elastic finite elements for one periodic channel sector (transverse r-theta plane at fixed x).

Mesh   : bilinear Q4 elements on a polar grid; the channel is removed pixel-wise (element centroid inside the
         channel). The helical channel cut by a transverse plane is stretched circumferentially by 1/cos(beta).
Thermal: steady conduction; gas convection (hg, Taw) on r = r_in; coolant convection (hc_eff, Tb) on the channel
         wall; outer surface and sector sides adiabatic.
Stress : linear elastic, generalized plane strain (uniform axial strain eps_z solved from the axial force balance
         sum(sigma_zz dA) = Fz); sector sides u_theta = 0 (periodic symmetry, penalty); coolant pressure on the
         channel wall, gas pressure on the hot wall. Load cases 'th' (thermal), 'pr' (pressure + end load), 'tot'.
Units  : mm, MPa, N, K.
"""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from .channel import TEAR_H


class Material:
    def __init__(self, E=120e3, nu=0.33, alpha=17e-6, k=320.0, T_ref=293.15, name="CuCrZr (generic)"):
        self.E, self.nu, self.alpha, self.k, self.T_ref, self.name = E, nu, alpha, k, T_ref, name

    @classmethod
    def from_cfg(cls, c):
        return cls(E=c.get("E_GPa", 120.0) * 1e3, nu=c.get("nu", 0.33), alpha=c.get("alpha_per_K", 17e-6),
                   k=c.get("k_W_mK", 320.0), T_ref=c.get("T_ref_K", 293.15), name=c.get("name", "wall material"))


def lame(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return lam, mu, np.array([[lam + 2 * mu, lam, 0], [lam, lam + 2 * mu, 0], [0, 0, mu]])


def radial_nodes(r_in, r_out, r_fine, h, hmax):
    r = list(np.linspace(r_in, r_fine, int(math.ceil((r_fine - r_in) / h)) + 1))
    d = h
    while r[-1] < r_out - 1e-9:
        d = min(d * 1.12, hmax)
        r.append(min(r[-1] + d, r_out))
    if len(r) > 2 and r[-1] - r[-2] < 0.3 * d:
        r.pop(-2)
    return np.array(r)


def mesh(r_in, r_out, dth, t, shape, s, beta_deg=0.0, we=None, h=0.015, hmax=0.08, channel=True):
    """Sector mesh. t: hot-wall thickness, s: channel size (see channel.py), we: rect channel width [mm]."""
    apex = t + s
    if r_in + apex >= r_out:
        raise ValueError("outer radius is inside the channel")
    rn = radial_nodes(r_in, r_out, min(r_in + apex + 0.3, r_out), h, hmax)
    nt = max(8, int(math.ceil(r_in * dth / h)))
    th = np.linspace(-dth / 2, dth / 2, nt + 1)
    nr = len(rn) - 1
    RR, TT = np.meshgrid(rn, th, indexing="ij")
    XN, YN = (RR * np.sin(TT)).ravel(), (RR * np.cos(TT)).ravel()
    rc = 0.5 * (rn[:-1] + rn[1:])
    tc = 0.5 * (th[:-1] + th[1:])
    RC, TC = np.meshgrid(rc, tc, indexing="ij")
    Xl = RC * np.sin(TC) * math.cos(math.radians(beta_deg))    # normal-to-channel coordinate
    Yl = RC * np.cos(TC) - r_in                                # height above the hot wall
    if shape == "tear":
        R = s / TEAR_H
        yc0 = t + R
        fluid = (Xl ** 2 + (Yl - yc0) ** 2 < R ** 2) | ((Yl >= yc0 + R / math.sqrt(2)) & (np.abs(Xl) <= apex - Yl))
    elif shape == "circ":
        R = s / 2
        fluid = Xl ** 2 + (Yl - t - R) ** 2 < R ** 2
    elif shape == "rect":
        if we is None:
            raise ValueError("rect channel needs its width we")
        fluid = (np.abs(Xl) <= we / 2) & (Yl >= t) & (Yl <= t + s)
    else:
        raise ValueError(shape)
    if not channel:
        fluid[:] = False
    I, J = np.nonzero(~fluid)

    def nid(i, j):
        return i * (nt + 1) + j

    conn = np.stack([nid(I, J), nid(I, J + 1), nid(I + 1, J + 1), nid(I + 1, J)], 1)
    return dict(XN=XN, YN=YN, conn=conn, I=I, J=J, fluid=fluid, nr=nr, nt=nt, rn=rn, th=th, apex=apex,
                r_in=r_in, r_out=r_out, dth=dth, shape=shape)


G3 = 1 / math.sqrt(3)
GPS = [(-G3, -G3), (G3, -G3), (G3, G3), (-G3, G3)]


def shape_fn(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return N, dxi, deta


def grads(Xe, Ye, xi, eta):
    N, dxi, deta = shape_fn(xi, eta)
    J11, J12, J21, J22 = Xe @ dxi, Ye @ dxi, Xe @ deta, Ye @ deta
    det = J11 * J22 - J12 * J21
    dx = (J22[:, None] * dxi - J12[:, None] * deta) / det[:, None]
    dy = (-J21[:, None] * dxi + J11[:, None] * deta) / det[:, None]
    return N, dx, dy, det


def Bmat(dx, dy):
    B = np.zeros((dx.shape[0], 3, 8))
    B[:, 0, 0::2] = dx
    B[:, 1, 1::2] = dy
    B[:, 2, 0::2] = dy
    B[:, 2, 1::2] = dx
    return B


def edges(M):
    """Boundary edges (node_p, node_q, outward nx, ny, length) of the hot wall ('gas') and the channel ('cool')."""
    I, J, conn, fluid, nr, nt = M["I"], M["J"], M["conn"], M["fluid"], M["nr"], M["nt"]
    X, Y = M["XN"], M["YN"]
    out = {"gas": [], "cool": []}
    for e, (di, dj, a, b) in enumerate([(-1, 0, 0, 1), (0, 1, 1, 2), (1, 0, 2, 3), (0, -1, 3, 0)]):
        ii, jj = I + di, J + dj
        inr = (ii >= 0) & (ii < nr) & (jj >= 0) & (jj < nt)
        isf = np.zeros(len(I), bool)
        isf[inr] = fluid[ii[inr], jj[inr]]
        p, q = conn[:, a], conn[:, b]
        dxe, dye = X[q] - X[p], Y[q] - Y[p]
        L = np.hypot(dxe, dye)
        rec = np.stack([p, q, dye / L, -dxe / L, L], 1)
        if e == 0:
            out["gas"].append(rec[I == 0])
        out["cool"].append(rec[isf])
    return {k: np.concatenate(v) for k, v in out.items()}


def thermal(M, hg, Taw, hc_eff, Tb, k=320.0):
    """Nodal temperatures [K]. hg, hc_eff [W/m^2/K], k [W/m/K] (lengths in mm are converted)."""
    X, Y, conn = M["XN"], M["YN"], M["conn"]
    n = len(X)
    Xe, Ye = X[conn], Y[conn]
    Ke = np.zeros((len(conn), 4, 4))
    for xi, eta in GPS:
        _, dx, dy, det = grads(Xe, Ye, xi, eta)
        Ke += k * (dx[:, :, None] * dx[:, None, :] + dy[:, :, None] * dy[:, None, :]) * det[:, None, None]
    rows = np.repeat(conn, 4, 1).ravel()
    cols = np.tile(conn, (1, 4)).ravel()
    K = sp.coo_matrix((Ke.ravel(), (rows, cols)), shape=(n, n)).tocsr()
    f = np.zeros(n)
    Ed = edges(M)
    for key, h, Tinf in (("gas", hg, Taw), ("cool", hc_eff, Tb)):
        E = Ed[key]
        p, q = E[:, 0].astype(int), E[:, 1].astype(int)
        hl = h * E[:, 4] * 1e-3
        K = K + sp.coo_matrix((np.concatenate([hl / 3, hl / 3, hl / 6, hl / 6]),
                               (np.concatenate([p, q, p, q]), np.concatenate([p, q, q, p]))), shape=(n, n)).tocsr()
        np.add.at(f, p, hl * Tinf / 2)
        np.add.at(f, q, hl * Tinf / 2)
    used = np.unique(conn)
    T = np.full(n, np.nan)
    T[used] = spl.spsolve(K[used][:, used].tocsc(), f[used])
    return T, Ed


def structural(M, T, p_cool, p_gas, Fz_ext, mat):
    """Element-centroid stresses per load case: sxx, syy, szz, sxy, srr, stt, vm [MPa] and eps_z."""
    lam, mu, DM = lame(mat.E, mat.nu)
    alpha, TREF = mat.alpha, mat.T_ref
    X, Y, conn = M["XN"], M["YN"], M["conn"]
    n, ne = len(X), len(conn)
    Xe, Ye, Te = X[conn], Y[conn], T[conn]
    Ke = np.zeros((ne, 8, 8))
    fth = np.zeros((ne, 8))
    fez = np.zeros((ne, 8))
    GPd = []
    for xi, eta in GPS:
        N, dx, dy, det = grads(Xe, Ye, xi, eta)
        B = Bmat(dx, dy)
        Ke += np.einsum("eki,kl,elj,e->eij", B, DM, B, det)
        dT = Te @ N - TREF
        cth = -(3 * lam + 2 * mu) * alpha * dT
        fth -= (cth * det)[:, None] * (B[:, 0, :] + B[:, 1, :])
        fez -= (lam * det)[:, None] * (B[:, 0, :] + B[:, 1, :])
        GPd.append((B, det, dT))
    dof = np.stack([2 * conn, 2 * conn + 1], 2).reshape(ne, 8)
    rows = np.repeat(dof, 8, 1).ravel()
    cols = np.tile(dof, (1, 8)).ravel()
    K = sp.coo_matrix((Ke.ravel(), (rows, cols)), shape=(2 * n, 2 * n)).tocsr()
    F = np.zeros((2 * n, 3))                     # columns: thermal, pressure, unit eps_z
    np.add.at(F[:, 0], dof.ravel(), fth.ravel())
    np.add.at(F[:, 2], dof.ravel(), fez.ravel())
    Ed = edges(M)
    for key, p in (("cool", p_cool), ("gas", p_gas)):
        Ee = Ed[key]
        a, b = Ee[:, 0].astype(int), Ee[:, 1].astype(int)
        fx, fy = -p * Ee[:, 2] * Ee[:, 4] / 2, -p * Ee[:, 3] * Ee[:, 4] / 2
        for nd in (a, b):
            np.add.at(F[:, 1], 2 * nd, fx)
            np.add.at(F[:, 1], 2 * nd + 1, fy)
    used = np.unique(conn)
    kp = 1e6 * K.diagonal().max()
    jside = np.r_[np.arange(0, M["nr"] + 1) * (M["nt"] + 1), np.arange(0, M["nr"] + 1) * (M["nt"] + 1) + M["nt"]]
    jside = np.intersect1d(jside, used)
    thn = np.arctan2(X[jside], Y[jside])
    ex, ey = np.cos(thn), -np.sin(thn)
    r_ = np.concatenate([2 * jside, 2 * jside, 2 * jside + 1, 2 * jside + 1])
    c_ = np.concatenate([2 * jside, 2 * jside + 1, 2 * jside, 2 * jside + 1])
    K = K + sp.coo_matrix((kp * np.concatenate([ex * ex, ex * ey, ey * ex, ey * ey]), (r_, c_)),
                          shape=(2 * n, 2 * n)).tocsr()
    ud = np.sort(np.concatenate([2 * used, 2 * used + 1]))
    lu = spl.splu(K[ud][:, ud].tocsc())
    U = np.zeros((2 * n, 3))
    U[ud] = lu.solve(F[ud])

    def Nz(u, ez, therm):
        tot = 0.0
        for B, det, dT in GPd:
            eps = np.einsum("eij,ej->ei", B, u[dof])
            szz = lam * (eps[:, 0] + eps[:, 1]) + (lam + 2 * mu) * ez - (3 * lam + 2 * mu) * alpha * dT * therm
            tot += (szz * det).sum()
        return tot

    N_th, N_pr, N_ez = Nz(U[:, 0], 0, 1), Nz(U[:, 1], 0, 0), Nz(U[:, 2], 1, 0)
    ez = {"th": -N_th / N_ez, "pr": (Fz_ext - N_pr) / N_ez}
    ez["tot"] = ez["th"] + ez["pr"]
    _, dx0, dy0, det0 = grads(Xe, Ye, 0.0, 0.0)
    B0 = Bmat(dx0, dy0)
    dTc = Te.mean(1) - TREF
    xc, yc = Xe.mean(1), Ye.mean(1)
    thc = np.arctan2(xc, yc)
    rc = np.hypot(xc, yc)
    c, s_ = np.cos(thc), np.sin(thc)
    res = {}
    for case, w in (("th", (1, 0)), ("pr", (0, 1)), ("tot", (1, 1))):
        u = w[0] * U[:, 0] + w[1] * U[:, 1] + ez[case] * U[:, 2]
        eps = np.einsum("eij,ej->ei", B0, u[dof])
        therm = w[0]
        s0 = lam * ez[case] - (3 * lam + 2 * mu) * alpha * dTc * therm
        sxx = (lam + 2 * mu) * eps[:, 0] + lam * eps[:, 1] + s0
        syy = lam * eps[:, 0] + (lam + 2 * mu) * eps[:, 1] + s0
        sxy = mu * eps[:, 2]
        szz = lam * (eps[:, 0] + eps[:, 1]) + (lam + 2 * mu) * ez[case] - (3 * lam + 2 * mu) * alpha * dTc * therm
        d = dict(sxx=sxx, syy=syy, szz=szz, sxy=sxy, ez=ez[case])
        d["srr"] = s_ * s_ * sxx + c * c * syy + 2 * s_ * c * sxy
        d["stt"] = c * c * sxx + s_ * s_ * syy - 2 * s_ * c * sxy
        d["vm"] = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * sxy ** 2)
        res[case] = d
    return res, dict(xc=xc, yc=yc, rc=rc, thc=thc, area=det0 * 4)
