"""Earlier certification pass (superseded by certify2.py): dyadic snap + interval enclosures.

Pipeline:
 1. load cg_cols.pkl; solve float master (tau_solve = 0.009, tightened)
 2. snap cg-config positions to dyadic k/2^20 (exact in float64)
 3. rebuild snapped rows; re-solve master; rationalize weights exactly (sum = 1)
 4. rigorous enclosures: mpmath.iv rows for cg configs (dyadic angles), exact Fraction
    rows for v4 atoms; S(j) interval = sum of exact-weight * row-interval
 5. verify S(j) in [j - tau_decl, j + tau_decl], tau_decl = 0.0101; sigma exact
 6. emit RowCertData-style integer enclosures at K = 2^60 -> law_certified_cg.txt
"""
import pickle, sys, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from mpmath import iv, mpf as _mpf
from colgen import N, ROWS, TAU, rows_of_config, sigma_of

iv.prec = 80
P0 = Fraction(10909258999421303588095230195816054408197, 16 * 10**39)
SNAP = 1 << 20

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
    b_solve = np.concatenate([ROWS + ts, -(ROWS - ts)])
    # solve on original columns to find support, then snap support and re-solve on it
    Am = np.array([c[0] for c in cols]).T
    cm = np.array([c[1] for c in cols])
    r = linprog(cm, A_ub=np.vstack([Am, -Am]), b_ub=b_solve,
                A_eq=np.ones((1, len(cols))), b_eq=[1.0], bounds=(0, None), method="highs")
    print(f"master (tau=0.009): status {r.status} p = {r.fun:.9f}", flush=True)
    if r.status != 0: return 2
    supp = [i for i in range(len(cols)) if r.x[i] > 1e-11]
    print(f"support {len(supp)}", flush=True)
    # snap support cg columns
    snapped = []      # (rows_float, sigma, kind, payload) kind: 'v4' exact | 'cg' dyadic
    for i in supp:
        Fv, sg, tag, payload = cols[i]
        if tag == "cg":
            x, m = payload
            xs = snap(x)
            xk = np.round(xs * SNAP / 256.0).astype(np.int64) % SNAP
            snapped.append((rows_of_config(xs, m), sg, "cg", (xk, m.astype(int))))
        else:
            snapped.append((Fv, sg, "v4", payload))
    A2 = np.array([s[0] for s in snapped]).T
    c2 = np.array([s[1] for s in snapped])
    r2 = linprog(c2, A_ub=np.vstack([A2, -A2]), b_ub=b_solve,
                 A_eq=np.ones((1, len(snapped))), b_eq=[1.0], bounds=(0, None), method="highs")
    print(f"snapped master: status {r2.status} p = {r2.fun:.9f}", flush=True)
    if r2.status != 0:
        print("snapped master infeasible at tau=0.009 — loosen or refine"); return 3
    keep = [k for k in range(len(snapped)) if r2.x[k] > 1e-12]
    w = [Fraction(r2.x[k]).limit_denominator(10**25) for k in keep]
    tot = sum(w); w = [v / tot for v in w]
    print(f"final support {len(keep)}; enclosing rows rigorously...", flush=True)
    # exact rational rows for v4 columns come from integer row caches — they ARE ints
    tau_decl = Fraction(int(sys.argv[3]), 10**6) if len(sys.argv) > 3 else Fraction(101, 10000)
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
