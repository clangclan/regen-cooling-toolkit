"""Matplotlib figures (default style, one quantity per axis)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _throat(ax, contour):
    ax.axvline(contour.x_t, color="0.6", lw=0.8, ls="--")


def geometry(case, res, path):
    G, c = res["G"], case.contour
    x = G["x"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    a = ax[0, 0]
    a.plot(x, G["r"], label="hot wall")
    a.plot(x, G["r"] + G["t"] * G["ca"], label="channel bottom")
    a.plot(x, G["r"] + (G["t"] + G["s"]) * G["ca"], label="channel top")
    a.set_ylabel("radius [mm]")
    a.legend()
    a = ax[0, 1]
    a.plot(x, G["t"], label="hot wall t")
    a.plot(x, G["s"], label="channel size s")
    a.plot(x, res["web"], label="web / rib")
    a.set_ylabel("[mm]")
    a.legend()
    a = ax[1, 0]
    a.plot(x, G["beta"])
    a.set_ylabel("helix angle [deg]")
    a.set_xlabel("x [mm]")
    a = ax[1, 1]
    a.plot(x, res["Dh"])
    a.set_ylabel("hydraulic diameter [mm]")
    a.set_xlabel("x [mm]")
    for a in ax.ravel():
        _throat(a, c)
        a.grid(alpha=0.3)
    _save(fig, path)


def gas(case, path):
    x = np.arange(case.contour.x0, case.contour.x_end + 1e-9, case.dx)
    g = case.gas.on_grid(x)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    ax[0, 0].plot(x, g["M"])
    ax[0, 0].set_ylabel("Mach number")
    ax[0, 1].plot(x, g["p"] / 1e5)
    ax[0, 1].set_ylabel("static pressure [bar]")
    ax[1, 0].plot(x, g["Taw"], label="adiabatic wall")
    ax[1, 0].plot(x, g["T"], label="static")
    ax[1, 0].set_ylabel("gas temperature [K]")
    ax[1, 0].legend()
    ax[1, 1].plot(x, g["h0"] / 1e3)
    ax[1, 1].set_ylabel("Bartz h0 (sigma = 1) [kW/m2/K]")
    for a in ax.ravel():
        _throat(a, case.contour)
        a.grid(alpha=0.3)
    for a in ax[1]:
        a.set_xlabel("x [mm]")
    _save(fig, path)


def solution(case, res, path, Tlim_C=None):
    x = res["x"]
    fig, ax = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    a = ax[0, 0]
    a.plot(x, res["Twg"] - 273.15, label="Twg (hot wall)")
    a.plot(x, res["Twc"] - 273.15, label="Twc (coolant wall)")
    a.plot(x, res["Tb"] - 273.15, label="Tb (coolant bulk)")
    if Tlim_C:
        a.axhline(Tlim_C, color="k", lw=0.8, ls=":", label="limit")
    a.set_ylabel("temperature [C]")
    a.legend()
    ax[0, 1].plot(x, res["q"] / 1e6)
    ax[0, 1].set_ylabel("hot-wall heat flux [MW/m2]")
    ax[1, 0].plot(x, res["P"] / 1e5)
    if case.fluid.Pc:
        ax[1, 0].axhline(case.fluid.Pc / 1e5, color="k", lw=0.8, ls=":", label="coolant Pc")
        ax[1, 0].legend()
    ax[1, 0].set_ylabel("coolant static pressure [bar]")
    ax[1, 1].semilogy(x, res["hg"], label="hg (gas)")
    ax[1, 1].semilogy(x, res["hc"], label="hc (coolant)")
    ax[1, 1].set_ylabel("heat-transfer coefficient [W/m2/K]")
    ax[1, 1].legend()
    ax[2, 0].plot(x, res["v"])
    ax[2, 0].set_ylabel("coolant velocity [m/s]")
    ax[2, 1].semilogy(x, res["Re"])
    ax[2, 1].axhline(2300, color="k", lw=0.8, ls=":")
    ax[2, 1].set_ylabel("Reynolds number")
    for a in ax.ravel():
        _throat(a, case.contour)
        a.grid(alpha=0.3)
    for a in ax[2]:
        a.set_xlabel("x [mm] (coolant flows from right to left)")
    _save(fig, path)


def station_field(st, path):
    F = st["field"]
    X, Y, conn = F["XN"], F["YN"], F["conn"]
    Tn = F["T"] - 273.15
    verts = np.stack([X[conn], Y[conn]], -1)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    for a, val, lab in ((ax[0], np.nanmean(Tn[conn], 1), "temperature [C]"), (ax[1], F["vm"], "von Mises [MPa]")):
        pc = PolyCollection(verts, array=val, cmap="inferno" if a is ax[0] else "viridis", edgecolors="face",
                            linewidths=0.05)
        a.add_collection(pc)
        a.autoscale_view()
        a.set_aspect("equal")
        a.set_xlabel("circumferential [mm]")
        a.set_ylabel("radial [mm]")
        fig.colorbar(pc, ax=a, label=lab, shrink=0.85)
    _save(fig, path)


def optimization(case, before, after, path):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for res, lab in ((before, "start"), (after, "optimized")):
        ax[0].plot(res["x"], res["G"]["s"], label=lab)
        ax[1].plot(res["x"], res["Twg"] - 273.15, label=lab)
    ax[0].set_ylabel("channel size s [mm]")
    ax[1].set_ylabel("hot-wall temperature [C]")
    for a in ax:
        _throat(a, case.contour)
        a.set_xlabel("x [mm]")
        a.grid(alpha=0.3)
        a.legend()
    _save(fig, path)


def tolerance(rows, path):
    lab = [r["case"] for r in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(1, 3, figsize=(13, 0.45 * len(rows) + 1.8), sharey=True)
    for a, k, xl in ((ax[0], "Twg", "max hot-wall temperature [C]"), (ax[1], "dP", "pressure drop [bar]"),
                     (ax[2], "Pmin", "channel exit pressure [bar]")):
        a.axvline(rows[0][k], color="0.6", lw=0.8, ls="--")       # rows[0] = nominal
        a.plot([r[k] for r in rows], y, "o", color="C0")
        a.set_xlabel(xl)
        a.grid(alpha=0.3)
        a.margins(x=0.15)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(lab)
    ax[0].invert_yaxis()
    _save(fig, path)


def manifold(ang, off, mf_nom, mc, path):
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    o = np.argsort(ang)
    ax[0].plot(ang[o], off[o], "o-")
    ax[0].set_xlabel("channel angle [deg]")
    ax[0].set_ylabel("header static offset [bar]")
    ax[1].plot(ang[o], mf_nom[o], "o-")
    ax[1].set_xlabel("channel angle [deg]")
    ax[1].set_ylabel("flow factor (nominal channels)")
    for m in mc:
        ax[2].hist(m["_Tmax"], bins=40, histtype="step", label="sigma %.3f mm" % m["sigma"])
    ax[2].set_xlabel("max hot-wall temperature over channels [C]")
    ax[2].set_ylabel("samples")
    ax[2].legend()
    for a in ax:
        a.grid(alpha=0.3)
    _save(fig, path)


def stress_summary(stations, path):
    x = [s["x"] for s in stations]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    w = 0.25
    i = np.arange(len(x))
    for k, (case, lab) in enumerate((("th", "thermal"), ("pr", "pressure"), ("tot", "total"))):
        ax[0].bar(i + (k - 1) * w, [s[case]["vm_max"] for s in stations], w, label=lab)
    ax[0].set_ylabel("max von Mises [MPa]")
    ax[1].bar(i - w / 2, [s["Twg_max"] for s in stations], w, label="2D")
    ax[1].bar(i + w / 2, [s["Twg_1d"] for s in stations], w, label="1D")
    ax[1].set_ylabel("max hot-wall temperature [C]")
    for a in ax:
        a.set_xticks(i)
        a.set_xticklabels(["x = %.1f" % v for v in x])
        a.legend()
        a.grid(alpha=0.3, axis="y")
    _save(fig, path)
