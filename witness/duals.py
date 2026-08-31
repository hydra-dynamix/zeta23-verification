"""Correct dual extraction for the CG master — validated, not guessed.

The old inline heuristic tried four sign conventions and accepted the first whose reduced
cost was small ON A SINGLE COLUMN (the heaviest basic one). With a 256-column support that
test is far too weak: a wrong convention can pass it by luck, after which every reduced cost
the pricer computes is meaningless. Diagnosed 2026-08-24: the pricer had been running on a
convention giving basic-column |rc| up to 2.0 instead of ~1e-14.

This version scores all four conventions against EVERY basic column and refuses to return
one that does not satisfy complementary slackness.
"""
import numpy as np

def extract(r, A, cost, tol=1e-6, f256=None):
    """Return (q, mu, report) with q the net row dual and mu the convexity dual.

    The master now carries a 511th inequality row enforcing the D(1) bound
    (sum_c w_c F_c(256) <= F256_CAP).  Its dual lam must enter the reduced cost:

        rc(c) = sigma_c + q.F_c + lam*F_c(256) + mu

    or the pricer optimises the wrong objective.  When f256 is supplied the lam term is
    scored alongside the band rows and returned in report["lam"].

    Raises RuntimeError if no sign convention satisfies rc ~ 0 on all basic columns,
    rather than silently returning a wrong one.
    """
    mu0 = float(r.eqlin.marginals[0])
    marg = np.asarray(r.ineqlin.marginals)
    n = A.shape[0]
    lam0 = 0.0
    if marg.size == 2*n + 1:
        lam0 = float(marg[2*n])
        marg = marg[:2*n]
    elif marg.size == 2*n + 2:
        # Two-sided D(1): rows are  -f256 . w <= -floor  and  f256 . w <= cap.  Their duals
        # act on the same column quantity with opposite sign, so the NET coefficient on
        # F(256) is their difference. The overall sign convention is settled by the
        # eight-combination calibration below, exactly as for the one-sided case.
        lam0 = float(marg[2*n]) - float(marg[2*n + 1])
        marg = marg[:2*n]
    q0 = marg[:n] - marg[n:]
    basic = np.where(r.x > 1e-11)[0]
    if basic.size == 0:
        raise RuntimeError("master has empty support; cannot calibrate duals")

    best = None
    report = {}
    # The D(1)/F256 row is written as -f256 . w <= -floor, so its marginal carries an
    # INDEPENDENT sign relative to the band rows.  Tying lam's sign to q's tried only 4 of
    # the 8 combinations and raised RuntimeError whenever that row was ACTIVE (max|rc| 0.65).
    # It never fired here because lam = 0 while the row is slack -- a latent bug found by an
    # adversarial review that ran with that row binding.  Enumerate all 8.
    combos = []
    for nm, qv, muv in (("q,mu", q0, mu0), ("-q,mu", -q0, mu0),
                        ("q,-mu", q0, -mu0), ("-q,-mu", -q0, -mu0)):
        combos.append((nm, qv, muv, lam0))
        combos.append((nm + ",-lam", qv, muv, -lam0))
    for name, qv, muv, lamv in combos:
        rcs = cost[basic] + qv @ A[:, basic] + muv
        if f256 is not None:
            rcs = rcs + lamv*np.asarray(f256)[basic]
        worst = float(np.abs(rcs).max())
        report[name] = worst
        if best is None or worst < best[0]:
            best = (worst, name, qv, muv, lamv)

    worst, name, qv, muv, lamv = best
    report["lam"] = lamv
    if worst > tol:
        raise RuntimeError(
            f"no sign convention satisfies complementary slackness "
            f"(best {name}: max|rc| on basic columns = {worst:.3e} > {tol:.0e}); "
            f"all conventions: {report}"
        )
    return qv, muv, {"convention": name, "basic_max_rc": worst,
                     "lam": report.get("lam", 0.0), "all": report}
