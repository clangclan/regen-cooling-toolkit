"""Case = one engine + coolant + channel design, read from a TOML file (see docs/config_reference.md)."""
import hashlib
import json
import math
import os
import tomllib

import numpy as np

from . import fluids, solver1d
from .channel import ChannelDesign, knots_fun
from .contour import Contour
from .fe2d import Material
from .gas import CanteraGas, IdealGas


def profile(v, log=False):
    """Number -> constant, [[x, value], ...] -> PCHIP through the points (constant outside)."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        c = float(v)
        return lambda x: np.full_like(np.asarray(x, float), c)
    return knots_fun(v, log=log)


class Case:
    def __init__(self, cfg, base_dir=".", out_dir=None, log=print):
        self.cfg, self.base, self.log = cfg, base_dir, log
        c = cfg.get("case", {})
        self.name = c.get("name", "case")
        self.dx = float(c.get("dx_mm", 0.25))
        self.C_hg = float(c.get("C_hg", 1.0))
        self.p_amb = float(c.get("p_ambient_bar", 1.01325)) * 1e5
        self.out_dir = out_dir or os.path.join("results", self.name)
        self.contour = self._contour(cfg["contour"])
        k = cfg["coolant"]
        self.mdot = float(k["mdot_kg_s"])
        self.Tin = float(k.get("T_in_K", 293.15))
        self.Pin = float(k["P_in_bar"]) * 1e5
        self.ks = float(k.get("roughness_um", 30.0)) * 1e-6
        self.K_in, self.K_out = float(k.get("K_in", 0.3)), float(k.get("K_out", 1.5))
        if k.get("fluid", "isopropanol").lower() == "csv":
            k = dict(k, table=self._path(k["table"]))
        self.fluid = fluids.make(k, self.Pin)
        self.material = Material.from_cfg(cfg.get("material", {}))
        lim = cfg.get("limits", {})
        self.dP_max = lim.get("dP_max_bar")
        pmin_default = self.fluid.Pc / 1e5 + 1.0 if self.fluid.Pc else None
        self.P_min = lim.get("P_min_bar", pmin_default)
        self.dh_fuel = 0.0
        self._gas = None
        self._gg = None

    # ------------------------------------------------------------------ construction
    @classmethod
    def load(cls, path, out_dir=None, log=print):
        with open(path, "rb") as f:
            cfg = tomllib.load(f)
        return cls(cfg, os.path.dirname(os.path.abspath(path)), out_dir, log)

    def _path(self, p):
        return p if os.path.isabs(p) else os.path.join(self.base, p)

    def _contour(self, c):
        kind = c.get("type", "parametric")
        rc = c.get("rc_throat_mm")
        if kind == "csv":
            return Contour.from_csv(self._path(c["file"]), rc, c.get("x_col", "x_mm"), c.get("r_col", "r_mm"))
        if kind == "parametric":
            return Contour.parametric(
                c["r_chamber_mm"], c["r_throat_mm"], c["l_cyl_mm"], c.get("conv_half_angle_deg", 30.0),
                c.get("r1_ratio", 1.5), c.get("r2_ratio", 1.0), c.get("rn_ratio", 0.382), c.get("eps_exit", 4.0),
                c.get("divergent", "bell"), c.get("theta_n_deg", 25.0), c.get("theta_e_deg", 10.0),
                c.get("bell_fraction", 0.8), rc_throat=rc)
        raise ValueError("contour.type must be 'parametric' or 'csv'")

    # ------------------------------------------------------------------ hot gas
    @property
    def gas(self):
        if self._gas is None:
            self._gas = self._make_gas()
            self._gg = None
        return self._gas

    def set_dh_fuel(self, dh):
        self.dh_fuel = float(dh)
        self._gas = None

    def _make_gas(self):
        g = self.cfg["gas"]
        pc = float(g["pc_bar"]) * 1e5
        model = g.get("model", "ideal")
        if model == "ideal":
            i = g.get("ideal", {})
            return IdealGas(self.contour, pc, i["T0_K"], i["gamma"], i["MW"], Pr=i.get("Pr"), mu=i.get("mu_Pa_s"),
                            cstar=g.get("cstar_m_s"), eta_cstar=g.get("eta_cstar", 1.0))
        if model == "cantera":
            table = self.gas_table()
            return CanteraGas(self.contour, table, pc, mode=g["cantera"].get("mode", "eq"), cstar=g.get("cstar_m_s"),
                              eta_cstar=g.get("eta_cstar", 1.0))
        raise ValueError("gas.model must be 'ideal' or 'cantera'")

    def gas_table(self):
        from . import cantera_gas
        g = self.cfg["gas"]
        c = g["cantera"]
        pc = float(g["pc_bar"]) * 1e5
        eps = self.contour.eps_max
        key = json.dumps([c, pc, round(self.dh_fuel, -2), round(eps, 3)], sort_keys=True)
        fn = os.path.join(self.out_dir, "gas_cache", "ct_%s.json" % hashlib.sha1(key.encode()).hexdigest()[:12])
        if os.path.exists(fn):
            with open(fn) as f:
                return json.load(f)
        self.log("  Cantera: chamber equilibrium + expansion tables (dh_fuel %.0f kJ/kg) ..." % (self.dh_fuel / 1e3))
        tab = cantera_gas.compute(pc, c["of"], c["fuel"], c["oxidizer"], eps, c.get("mechanism", "gri30_highT.yaml"),
                                  self.dh_fuel, c.get("npts", 600), modes=(c.get("mode", "eq"),))
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        with open(fn, "w") as f:
            json.dump(tab, f)
        return tab

    def gas_grid(self):
        if self._gg is None:
            x = np.arange(self.contour.x0, self.contour.x_end + 1e-9, self.dx)
            self._gg = self.gas.on_grid(x)
        return self._gg

    @property
    def regen_coupled(self):
        g = self.cfg["gas"]
        return g.get("model") == "cantera" and g["cantera"].get("regen_coupling", False)

    def couple(self, iterations=3, tol=5e3):
        """Cantera + regen coupling: the fuel enthalpy pickup in the jacket raises the chamber enthalpy.
        Iterates dh_fuel = dh_coolant * mdot_coolant / mdot_fuel. Returns the last 1D result."""
        res = self.solve()
        if not self.regen_coupled:
            return res
        mf = self.cfg["gas"]["cantera"].get("fuel_mdot_kg_s", self.mdot)
        for _ in range(iterations):
            dh = res["dh_coolant"] * self.mdot / mf
            if abs(dh - self.dh_fuel) < tol:
                break
            self.set_dh_fuel(dh)
            res = self.solve()
        self.log("  regen coupling: dh_fuel %.0f kJ/kg, T0 %.0f K" % (self.dh_fuel / 1e3, self.gas.T0))
        return res

    def nozzle_flow_check(self):
        """Exit pressure vs. ambient: Summerfield (pe/pa < ~0.25-0.4) and Schmucker p_sep = pa (1.88 Me - 1)^-0.64."""
        g = self.gas.on_grid([self.contour.x_end])
        pe, Me = float(g["p"][0]), float(g["M"][0])
        psep = self.p_amb * (1.88 * Me - 1) ** -0.64 if Me > 0.6 else float("nan")
        return dict(pe_bar=pe / 1e5, Me=Me, pe_over_pa=pe / self.p_amb, p_sep_schmucker_bar=psep / 1e5,
                    separated_schmucker=bool(pe < psep), summerfield_risk=bool(pe / self.p_amb < 0.4))

    # ------------------------------------------------------------------ channels
    def design(self, d=0.0, s_fun=None, t_fun=None):
        ch = self.cfg["channels"]
        beta = ch.get("helix_deg", 0.0)
        beta_fun = None if (isinstance(beta, (int, float)) and beta == 0) else profile(beta)
        D = ChannelDesign(t_fun or profile(ch["t_mm"]), s_fun or profile(ch["s_mm"]), ch["n"],
                          shape=ch.get("shape", "tear"), rib=ch.get("rib_mm", 1.0), land_deg=ch.get("land_deg", 0.0),
                          n_land=ch.get("n_land", 0), beta_fun=beta_fun, xs=ch.get("split_x_mm"),
                          rib_s=ch.get("split_rib_mm"), w_split=ch.get("w_split_mm"), name=self.name)
        D = D.clamp_web(self.contour, ch.get("web_min_mm", 0.0), self.dx)
        return D.offset(d) if d else D

    # ------------------------------------------------------------------ solve
    def solve(self, D=None, mdot=None, C_hg=None, G=None):
        return solver1d.solve(D or self.design(), self.contour, self.gas, self.fluid,
                              self.mdot if mdot is None else mdot, self.Tin, self.Pin, self.ks, self.K_in, self.K_out,
                              self.C_hg if C_hg is None else C_hg, self.material.k, self.dx, G=G, GG=self.gas_grid())

    def metrics(self, res):
        return solver1d.metrics(res, self.contour, self.fluid)

    def stations(self, res):
        st = self.cfg.get("stress", {})
        xs = []
        for v in st.get("x_mm", ["throat"]):
            if v == "throat":
                xs.append(self.contour.x_t)
            elif v == "max_Twg":
                xs.append(float(res["x"][int(np.argmax(res["Twg"]))]))
            else:
                xs.append(float(v))
        return xs

    def describe(self):
        lines = [self.contour.describe(), self.gas.describe(),
                 "coolant: %s, mdot %.4f kg/s, Tin %.1f K, Pin %.2f bar, ks %.0f um, K_in %.2f, K_out %.2f"
                 % (self.fluid.name, self.mdot, self.Tin, self.Pin / 1e5, self.ks * 1e6, self.K_in, self.K_out)]
        D = self.design()
        lines.append("channels: %d x %s, Neff %.2f, lands %d x %.1f deg, helix %s" % (
            D.nch, D.shape, D.Neff, D.NL, D.land, "yes" if D.beta_fun else "axial"))
        if self.fluid.Pc:
            lines.append("coolant critical pressure %.2f bar -> P_min limit %s bar" % (
                self.fluid.Pc / 1e5, "%.2f" % self.P_min if self.P_min else "-"))
        return "\n".join(lines)


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items() if k not in ("G", "field")}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o
