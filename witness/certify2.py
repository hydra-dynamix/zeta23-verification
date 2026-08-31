"""Rigorous certification: dyadic snap + exact rational weights + interval enclosures.

Pipeline:
 1. load cg_cols.pkl; solve float master (tau_solve = 0.009, tightened)
 2. snap cg-config positions to dyadic k/2^20 (exact in float64)
 3. rebuild snapped rows; re-solve master; rationalize weights exactly (sum = 1)
 4. rigorous enclosures: mpmath.iv rows for cg configs (dyadic angles), exact Fraction
    rows for v4 atoms; S(j) interval = sum of exact-weight * row-interval
 5. verify S(j) in [j - tau_decl, j + tau_decl], tau_decl = 0.0101; sigma exact
 6. emit RowCertData-style integer enclosures at K = 2^60 -> law_certified_cg.txt
"""
import decimal, pickle, sys, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from mpmath import iv, mpf as _mpf
from colgen import N, ROWS, TAU, rows_of_config, sigma_of

iv.prec = int(__import__('os').environ.get('CERT_PREC', 80))
P0 = Fraction(10909258999421303588095230195816054408197, 16 * 10**39)
# Snap resolution.  2^-20 moves a position by up to 256/2^21 = 1.22e-4, which is LARGER
# than tau = 1e-4 -- that is why every certificate emitted so far had to declare
# tau_decl = 1/5000 rather than the tau the law was actually solved at.  Finer snaps make a
# tau=1e-4 certificate reachable; the cost is only integer size in the exact arithmetic.
SNAP = 1 << int(__import__('os').environ.get('CERT_SNAP_BITS', 20))

def snap(x):
    return np.round(x * SNAP / 256.0) * (256.0 / SNAP)

def iv_rows(xk, m):
    """rigorous interval rows for dyadic positions xk (integers, x = 256*xk/SNAP)."""
    out = []
    two_pi = 2 * iv.pi
    for j in range(1, N + 1):
        zr = iv.mpf(0); zi = iv.mpf(0)
        for kpos, mk in zip(xk, m):
            ang = two_pi * (int(j) * int(kpos) % SNAP) / SNAP
            zr += int(mk) * iv.cos(ang)
            zi += int(mk) * iv.sin(ang)
        out.append(zr * zr + zi * zi)
    return out

def main():
    t0 = time.time()
    cols = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else "cg_cols.pkl", "rb"))
    print(f"{len(cols)} columns loaded", flush=True)
    ts = float(sys.argv[4]) if len(sys.argv) > 4 else 0.009
    # The D(1) condition is a CONSTRAINT of the theorem, not a post-hoc check.  Solving without
    # it lets the optimiser drift above the published cap -- E[F(256)] is a free direction for
    # reducing p, so it WILL drift -- and the violation is then only discovered by the row-256
    # gate after an hour of interval arithmetic.  Carry it in the master instead.
    D1_CAP = Fraction(82395317, 100000000)
    F256_HI_TRUE = float(Fraction(128) + D1_CAP * 65536)  # E[F(256)] <= 54126.59494912
    # Solve strictly inside the cap: snapping positions perturbs F(256), and the certificate
    # gate checks the SNAPPED law against the true cap.  Without margin the solve sits exactly
    # on the boundary and the snapped master goes infeasible.
    F256_MARGIN = float(__import__('os').environ.get('CERT_F256_MARGIN', 150.0))
    F256_HI = F256_HI_TRUE - F256_MARGIN
    F256_LO = 128.0                                       # D(1) >= 0
    def f256_of(payload, kind, dyadic=False):
        """E[F(256)] contribution of one column.

        NOTE the two position representations: pre-snap columns carry float positions in
        [0, 256), while SNAPPED columns carry dyadic integer indices xk with x = 256*xk/SNAP.
        Evaluating cos(2*pi*xk) on the indices gives 1 for every atom and hence F(256) = 65536
        for every column, which silently makes the cap infeasible.  Convert first.
        """
        if kind == "v4":
            return 0.0                                    # no positions: bound BELOW (a floor)
        x, m = payload
        x = np.asarray(x, float)
        if dyadic:
            x = x * (256.0 / SNAP)
        a = 2 * np.pi * x
        m = np.asarray(m, float)
        return float((np.cos(a) * m).sum() ** 2 + (np.sin(a) * m).sum() ** 2)
    b_solve = np.concatenate([ROWS + ts, -(ROWS - ts), [-F256_LO], [F256_HI]])
    # solve on original columns to find support, then snap support and re-solve on it
    Am = np.array([c[0] for c in cols]).T
    cm = np.array([c[1] for c in cols])
    fm = np.array([f256_of(c[3], c[2]) for c in cols])
    r = linprog(cm, A_ub=np.vstack([Am, -Am, -fm[None, :], fm[None, :]]), b_ub=b_solve,
                A_eq=np.ones((1, len(cols))), b_eq=[1.0], bounds=(0, None), method="highs")
    print(f"master (tau=0.009): status {r.status} p = {r.fun:.9f}", flush=True)
    if r.status != 0: return 2
    supp = [i for i in range(len(cols)) if r.x[i] > 1e-11]
    print(f"support {len(supp)} (float master)", flush=True)
    # Snap the WHOLE pool, not just the support.  Snapping to 2^-20 moves a position by up
    # to 1.22e-4, which perturbs the rows enough to break a tau=1e-4 band; re-solving over
    # only the 256 support columns then has no freedom left and goes infeasible.  Snapping
    # everything and re-solving over the full pool lets the LP pick a different, snap-stable
    # support.
    supp = list(range(len(cols)))
    # snap support cg columns
    snapped = []      # (rows_float, sigma, kind, payload) kind: 'v4' exact | 'cg' dyadic
    for i in supp:
        Fv, sg, tag, payload = cols[i]
        # Classify by PAYLOAD, not by tag string.  Tags proliferated across the campaign
        # (cg, v4, sfw, sfw2, sfw3, hop, loc, js, swp, rsh, s209..s219); keying on
        # tag == "cg" silently routed every newer column into the exact-integer branch,
        # which is wrong for continuous positions.  A column is dyadic-certifiable iff it
        # actually carries a (positions, marks) pair.
        has_pos = (isinstance(payload, tuple) and len(payload) == 2
                   and all(isinstance(z, np.ndarray) for z in payload))
        if has_pos:
            x, m = payload
            xs = snap(x)
            xk = np.round(xs * SNAP / 256.0).astype(np.int64) % SNAP
            snapped.append((rows_of_config(xs, m), sg, "cg", (xk, m.astype(int))))
        else:
            snapped.append((Fv, sg, "v4", payload))
    A2 = np.array([s[0] for s in snapped]).T
    c2 = np.array([s[1] for s in snapped])
    f2v = np.array([f256_of(s_[3], s_[2], dyadic=True) for s_ in snapped])
    r2 = linprog(c2, A_ub=np.vstack([A2, -A2, -f2v[None, :], f2v[None, :]]), b_ub=b_solve,
                 A_eq=np.ones((1, len(snapped))), b_eq=[1.0], bounds=(0, None), method="highs")
    print(f"snapped master: status {r2.status} p = "
          f"{r2.fun if r2.fun is None else format(r2.fun, '.9f')}", flush=True)
    if r2.status != 0:
        print(f"snapped master infeasible (status {r2.status}); "
              f"raise CERT_F256_MARGIN (now {F256_MARGIN}) or loosen tau_solve"); return 3
    keep = [k for k in range(len(snapped)) if r2.x[k] > 1e-12]
    w = [Fraction(r2.x[k]).limit_denominator(10**25) for k in keep]
    tot = sum(w); w = [v / tot for v in w]
    print(f"final support {len(keep)}; enclosing rows rigorously...", flush=True)
    # exact rational rows for v4 columns come from integer row caches — they ARE ints
    # tau_decl as an EXACT decimal. Historically this took an integer in units of 1e-6, which
    # cannot express the tighter bands the witness now supports. A decimal string is parsed
    # exactly by Fraction, so "1e-7" or "0.0000001" both give the exact rational.
    if len(sys.argv) > 3:
        raw = sys.argv[3]
        tau_decl = (Fraction(int(raw), 10**6) if raw.isdigit()
                    else Fraction(decimal.Decimal(raw)))
    else:
        tau_decl = Fraction(101, 10000)
    # --- bulk rigorous pass (restructured for speed) ---
    print("bulk interval evaluation...", flush=True)
    Slo = [Fraction(0)] * 257
    Shi = [Fraction(0)] * 257
    for k, sk in enumerate(keep):
        Fv, sg, kind, payload = snapped[sk]
        if kind == "v4":
            for j in range(1, 256):
                v = w[k] * Fraction(int(round(Fv[j - 1])))
                Slo[j] += v; Shi[j] += v
            # Row 256 is NOT stored for v4 columns (Fv has only 255 entries, j=1..255), so
            # the old loop silently contributed ZERO to Slo[256]/Shi[256].  Row 256 feeds
            # D(1) = E[F(256)]/65536 - 1/512, which the ceiling theorem needs bounded, so a
            # missing contribution is a wrong certificate, not a loose one.  Without the
            # positions we cannot evaluate F(256); enclose it by its full attainable range
            # [0, 65536] instead.  That keeps the certificate VALID (merely conservative):
            # at the current v4 weight ~7e-4 this widens D(1) by only +-7e-4.
            Slo[256] += w[k] * Fraction(0)
            Shi[256] += w[k] * Fraction(65536)
        else:
            xk, m = payload
            rows_iv = iv_rows(xk, m)
            for j in range(1, 257):
                fj = rows_iv[j - 1]
                lo_t, hi_t = fj._mpi_
                def tup2frac(t):
                    s, man, exp, bc = t
                    fr = Fraction(man, 1) * (Fraction(2) ** exp)
                    return -fr if s else fr
                flo = tup2frac(lo_t)
                fhi = tup2frac(hi_t)
                Slo[j] += w[k] * flo
                Shi[j] += w[k] * fhi
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(keep)} cols enclosed [{time.time()-t0:.0f}s]", flush=True)
    ok = True
    worst = Fraction(0)
    for j in range(1, 256):
        lo_dev = j - Slo[j]; hi_dev = Shi[j] - j
        worst = max(worst, abs(lo_dev), abs(hi_dev))
        if Slo[j] < j - tau_decl or Shi[j] > j + tau_decl:
            print(f"row {j}: enclosure [{float(Slo[j]):.6f},{float(Shi[j]):.6f}] outside band")
            ok = False
    # ---- row 256 / D(1) gate -------------------------------------------------------
    # AUDIT FIX (bug 14): the band loop above runs j = 1..255, so row 256 was WRITTEN into
    # every emitted certificate but never validated -- no `ok = False` could ever fire for
    # it, meaning no certificate this script produced actually certified the D(1) condition
    # it claims to.  D(1) = (E[F(256)] - 128) / 65536.
    # We gate BOTH sides, which is strictly more conservative than the one-sided floor the
    # LP has been carrying: the published hypothesis is the two-sided |D(1)| <= D1_CAP.
    D1_CAP = Fraction(82395317, 100000000)          # published bound, rounded up
    EF_HI  = Fraction(128) + D1_CAP * 65536
    EF_LO  = Fraction(128)                           # our floor: D(1) >= 0
    d1_lo = (Slo[256] - 128) / 65536
    d1_hi = (Shi[256] - 128) / 65536
    print(f"row 256: E[F(256)] in [{float(Slo[256]):.3f}, {float(Shi[256]):.3f}]  "
          f"=> D(1) in [{float(d1_lo):.9f}, {float(d1_hi):.9f}]", flush=True)
    if Slo[256] < EF_LO:
        print(f"row 256: E[F(256)] lower enclosure {float(Slo[256]):.3f} < {float(EF_LO):.3f} "
              f"(D(1) >= 0 NOT established)")
        ok = False
    if Shi[256] > EF_HI:
        print(f"row 256: E[F(256)] upper enclosure {float(Shi[256]):.3f} > {float(EF_HI):.3f} "
              f"(|D(1)| <= {float(D1_CAP)} NOT established)")
        ok = False
    # --------------------------------------------------------------------------------
    p_exact = sum(Fraction(int(round(snapped[sk][1] * N)), N) * w[k] for k, sk in enumerate(keep))
    print(f"bands OK: {ok}; worst dev {float(worst):.4e}; p = {float(p_exact):.9f}", flush=True)
    print(f"gap to p0 = {float(p_exact - P0):+.6f}", flush=True)
    if ok:
        K64 = 1 << 60
        with open(sys.argv[2] if len(sys.argv) > 2 else "law_certified_cg.txt", "w") as f:
            f.write(f"p = {p_exact}\np_float = {float(p_exact):.9f}\ngap = {float(p_exact-P0):+.6f}\n"
                    f"tau_decl = {tau_decl}\nworst_dev = {float(worst):.4e}\nsupport = {len(keep)}\n"
                    f"K = 2^60\nenclosures: integer (lo, hi) per row j=1..256\n\n")
            for j in range(1, 257):
                lo_i = int(Slo[j] * K64) - 1
                hi_i = int(Shi[j] * K64) + 1
                f.write(f"({lo_i}, {hi_i}),\n")
            f.write("\nsupport configurations:\n")
            for k, sk in enumerate(keep):
                Fv, sg, kind, payload = snapped[sk]
                if kind == "cg":
                    xk, m = payload
                    f.write(f"w={w[k]} [cg] doubles={int((np.array(m)==2).sum())} n={len(m)}\n")
                    f.write(f"  xk(2^-20 units)={list(map(int, xk))}\n  m={list(map(int, m))}\n")
                else:
                    f.write(f"w={w[k]} [v4] atoms={payload}\n")
        print("written law_certified_cg.txt", flush=True)
    print(f"done [{time.time()-t0:.0f}s]")
    return 0 if ok else 4

if __name__ == "__main__":
    sys.exit(main())
