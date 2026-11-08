"""Manufacturing-tolerance sensitivity.

uniform offset d (all channels): channel features +d, rect ribs -d, hot wall -d/2 (ChannelDesign.offset).
one channel at -d          : a single undersized channel sees the same manifold-to-manifold pressure drop as the
                             nominal channels, so its flow factor mf is found by bisection on dP(mf) = dP_nominal.
                             This is the hot-streak case (reduced flow -> highest wall temperature).
"""
from .solver1d import ChannelClosed


def matched_flow(case, D, dP_target, C_hg=None, lo=0.15, hi=1.3, n=22):
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        try:
            bad = case.solve(D, mdot=case.mdot * mid, C_hg=C_hg)["dP_tot"] > dP_target
        except (ChannelClosed, ValueError):
            bad = True
        if bad:
            hi = mid
        else:
            lo = mid
    return lo, case.metrics(case.solve(D, mdot=case.mdot * lo, C_hg=C_hg))


def study(case, offsets=(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2), one_channel=(0.05, 0.1, 0.2), C_hg=None, log=print):
    rows = []
    base = case.metrics(case.solve(C_hg=C_hg))
    rows.append(dict(case="nominal", d=0.0, flow=1.0, **base))
    log("  %-22s dP %5.2f  Pmin %6.2f  Twg %6.0f C @ %.1f" % ("nominal", base["dP"], base["Pmin"], base["Twg"],
                                                              base["x_Twg"]))
    for d in offsets:
        try:
            m = case.metrics(case.solve(case.design(d), C_hg=C_hg))
        except (ChannelClosed, ValueError) as e:
            log("  all channels %+.3f mm: failed (%s)" % (d, e))
            continue
        rows.append(dict(case="all %+.3f mm" % d, d=d, flow=1.0, **m))
        log("  %-22s dP %5.2f  Pmin %6.2f  Twg %6.0f C @ %.1f" % (rows[-1]["case"], m["dP"], m["Pmin"], m["Twg"],
                                                                  m["x_Twg"]))
    for d in one_channel:
        try:
            f, m = matched_flow(case, case.design(-abs(d)), base["dP"], C_hg)
        except (ChannelClosed, ValueError) as e:
            log("  one channel -%.3f mm: failed (%s)" % (d, e))
            continue
        rows.append(dict(case="one channel -%.3f mm" % abs(d), d=-abs(d), flow=f, **m))
        log("  %-22s flow x%.3f  Twg %6.0f C @ %.1f  Tout %.0f C" % (rows[-1]["case"], f, m["Twg"], m["x_Twg"],
                                                                     m["Tout"]))
    return rows
