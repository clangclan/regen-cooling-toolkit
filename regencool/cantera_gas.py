"""Combustion-gas properties with Cantera.

Chamber : adiabatic HP equilibrium at pc from the propellant feed enthalpies (+ optional regen heat pickup of the
          fuel). Frozen cp, viscosity and Prandtl number from Cantera's mixture-averaged transport.
Nozzle  : isentropic expansion s = s0 from pc downwards, either shifting equilibrium ('eq', SP equilibrium at
          every pressure) or frozen composition ('fr'). The throat is the pressure of maximum mass flux
          G = rho u with u = sqrt(2 (h0 - h)); the area ratio of every state is eps = G_t / G.
          Taw = T(h + r u^2/2, p) with the recovery factor r = Pr^(1/3).

Supported elements: C, H, O, N (the default mechanism is GRI-Mech 3.0 with high-temperature thermo).
"""
import math
import re

import numpy as np

ATOMIC = dict(C=12.011, H=1.008, O=15.999, N=14.007)


def _cantera():
    try:
        import cantera as ct
    except ImportError as e:  # pragma: no cover
        raise ImportError("gas model 'cantera' needs the cantera package:  pip install cantera") from e
    return ct


def parse_formula(formula):
    """'C3H8O' -> {'C': 3, 'H': 8, 'O': 1}."""
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        if el not in ATOMIC:
            raise ValueError("only C, H, O, N are supported in a formula (got %r)" % el)
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    if not out:
        raise ValueError("cannot parse formula %r" % formula)
    return out


def propellant(gas, spec):
    """Molar feed enthalpy [J/kmol], molar mass [kg/kmol] and atoms per molecule of one propellant.

    spec (dict) is either
      {species = "N2O", T_K = 293.15, dh_vap_kJ_kg = 170}   a species of the mechanism, gas-phase enthalpy at T_K
                                                           minus the vaporization enthalpy (0 for a gas feed)
      {formula = "C3H8O", hf_kJ_mol = -318.1, cp_J_kgK = 2600, T_K = 293.15}
                                                           liquid with its formation enthalpy at 298.15 K
    Cryogenic feeds: give T_K = 298.15 and dh_vap_kJ_kg = h_gas(298.15 K) - h_liquid(T_feed).
    """
    T = float(spec.get("T_K", 298.15))
    if "species" in spec:
        sp = spec["species"]
        gas.TPX = T, 1e5, sp + ":1"
        MW = gas.mean_molecular_weight
        h = gas.enthalpy_mole - float(spec.get("dh_vap_kJ_kg", 0.0)) * 1e3 * MW
        atoms = {e: gas.n_atoms(sp, e) for e in gas.element_names if gas.n_atoms(sp, e) > 0}
    elif "formula" in spec:
        atoms = parse_formula(spec["formula"])
        MW = sum(ATOMIC[e] * n for e, n in atoms.items())
        h = float(spec["hf_kJ_mol"]) * 1e6 + float(spec.get("cp_J_kgK", 0.0)) * MW * (T - 298.15)
    else:
        raise ValueError("propellant needs 'species' or 'formula'")
    return h, MW, atoms


def start_composition(el):
    """Product-like composition (CO, CO2, H2O, H2, O2, N2) with exactly the given element amounts."""
    extra = set(el) - set("CHON")
    if extra:
        raise ValueError("elements %s are not supported by the start composition" % sorted(extra))
    C, H, O, N = (el.get(k, 0.0) for k in "CHON")
    if O < C - 1e-12:
        raise ValueError("fewer O than C atoms: mixture too fuel-rich for the built-in start composition")
    CO, Orem = C, O - C
    H2O = min(H / 2, Orem)
    Orem -= H2O
    H2 = H / 2 - H2O
    CO2 = min(C, Orem)
    Orem -= CO2
    CO -= CO2
    X = dict(CO=CO, CO2=CO2, H2O=H2O, H2=H2, O2=Orem / 2, N2=N / 2)
    return {k: v for k, v in X.items() if v > 1e-12}


def compute(pc, of, fuel, oxidizer, eps_max, mechanism="gri30_highT.yaml", dh_fuel=0.0, npts=600,
            modes=("eq", "fr")):
    """Chamber state and expansion tables.

    pc [Pa], of = oxidizer/fuel mass ratio, dh_fuel [J/kg] enthalpy added to the fuel before injection
    (regenerative heating), eps_max: largest area ratio needed (the table stops a little beyond it).
    Returns a JSON-serializable dict (see CanteraGas for the keys used).
    """
    ct = _cantera()
    gas = ct.Solution(mechanism)
    hf, MWf, af = propellant(gas, fuel)
    hf += dh_fuel * MWf
    ho, MWo, ao = propellant(gas, oxidizer)
    n = of * MWf / MWo                                   # kmol oxidizer per kmol fuel
    H = (hf + n * ho) / (MWf + n * MWo)                  # J/kg
    el = {e: af.get(e, 0.0) + n * ao.get(e, 0.0) for e in set(af) | set(ao)}
    gas.TPX = 3000.0, pc, start_composition(el)
    gas.HP = H, pc
    gas.equilibrate("HP")
    st, X0 = gas.state.copy(), gas.X.copy()

    def frozen():
        cp, cv = gas.cp_mass, gas.cv_mass
        mu, k = gas.viscosity, gas.thermal_conductivity
        return dict(T=gas.T, MW=gas.mean_molecular_weight, cp_f=cp, gam_f=cp / cv, mu=mu, k_f=k, Pr_f=mu * cp / k,
                    rho=gas.density)

    def state_sp(s, p, X, eq):
        gas.SPX = s, p, X
        if eq:
            gas.equilibrate("SP")
        return gas.T, gas.enthalpy_mass, gas.density, gas.X.copy()

    def sound_eq(s, p, X):
        d = 1e-4
        r = [state_sp(s, p * f, X, True)[2] for f in (1 - d, 1 + d)]
        return math.sqrt(2 * d * p / (r[1] - r[0]))

    gas.state = st
    ch = frozen()
    T0, s0, h0 = gas.T, gas.entropy_mass, gas.enthalpy_mass
    hs = []
    for dT in (-5.0, 5.0):
        gas.TPX = T0 + dT, pc, X0
        gas.equilibrate("TP")
        hs.append(gas.enthalpy_mass)
    ch["cp_eq"] = (hs[1] - hs[0]) / 10.0
    gas.state = st
    ch["gam_s"] = gas.density * sound_eq(s0, pc, X0) ** 2 / pc
    gas.state = st
    top = sorted(zip(gas.species_names, gas.X), key=lambda t: -t[1])[:9]
    ch["X"] = {k: float(v) for k, v in top}
    out = dict(pc=pc, of=of, dh_fuel_kJkg=dh_fuel / 1e3, mechanism=mechanism, eps_max=eps_max, chamber=ch)
    rec = ch["Pr_f"] ** (1 / 3)

    for mode in modes:
        eq = mode == "eq"
        gas.state = st

        def G_of_p(p):
            T, h, rho, _ = state_sp(s0, p, X0, eq)
            return rho * math.sqrt(max(2 * (h0 - h), 0.0))

        a, b = 0.40 * pc, 0.80 * pc
        gr = (math.sqrt(5) - 1) / 2
        c, d = b - gr * (b - a), a + gr * (b - a)
        for _ in range(60):
            if G_of_p(c) > G_of_p(d):
                b = d
            else:
                a = c
            c, d = b - gr * (b - a), a + gr * (b - a)
        pt = 0.5 * (a + b)
        Gt = G_of_p(pt)
        ps = np.unique(np.r_[np.geomspace(pc * 0.9995, pc / 2000.0, npts), pt])[::-1]
        rows, X = [], X0.copy()
        for p in ps:
            T, h, rho, X = state_sp(s0, p, X, eq)
            fr = frozen()
            u = math.sqrt(max(2 * (h0 - h), 0.0))
            G = rho * u
            a_s = sound_eq(s0, p, X) if eq else math.sqrt(fr["gam_f"] * p / rho)
            gas.HPX = h + rec * 0.5 * u * u, p, X
            if eq:
                gas.equilibrate("HP")
            eps = Gt / G if G > 0 else 1e9
            rows.append(dict(p=p, T=T, rho=rho, u=u, eps=eps, M=u / a_s, gam_s=rho * a_s * a_s / p,
                             gam_f=fr["gam_f"], cp_f=fr["cp_f"], MW=fr["MW"], Taw=gas.T, sup=bool(p < pt)))
            if p < pt and eps > 1.05 * eps_max:
                break
        out["cstar_" + mode] = pc / Gt
        out["pt_" + mode] = pt
        out["exp_" + mode] = rows
    return out
