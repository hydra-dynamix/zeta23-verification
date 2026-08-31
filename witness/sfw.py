"""Sliding Frank-Wolfe on the witness LP: let the optimiser MOVE atoms, not just add them.

Classical CG never relocates a column once inserted, so near-miss configurations
accumulate instead of being repaired -- the standard diagnosis for a CG plateau
(Denoyelle-Duval-Peyre-Soubies, Inverse Problems 36 (2020) 014001).

LP sensitivity gives the slide direction in closed form.  With duals q (rc = sigma_c
+ q.A[:,c] + mu) and optimal weights w,

    d(obj)/dF_c(j) = w_c q_j        =>     d(obj)/dx_i^c = w_c * dE_c/dx_i

where E_c = sum_j q_j F_c(j) and, with A_j = sum_i m_i cos(2pi j x_i/256),
B_j = sum_i m_i sin(2pi j x_i/256),

    dE_c/dx_i = (4 pi m_i / 256) * sum_j q_j j [ B_j cos(2pi j x_i/256) - A_j sin(2pi j x_i/256) ]
"""
import sys, time, pickle
import numpy as np
from scipy.optimize import linprog
from colgen import ROWS, rows_of_config, sigma_of
from duals import extract

TWOPI = 2*np.pi

def spectrum(x, m):
    ang = TWOPI * np.outer(ROWS, x) / 256.0
    return (np.cos(ang)*m).sum(axis=1), (np.sin(ang)*m).sum(axis=1)

def dE_dx(x, m, q):
    """gradient of E = sum_j q_j F(j) w.r.t. each position; direct sum, no FFT."""
    A, B = spectrum(x, m)
    ang = TWOPI * np.outer(ROWS, x) / 256.0          # (255, n)
    qj = q * ROWS                                     # (255,)
    term = (qj[:, None] * (B[:, None]*np.cos(ang) - A[:, None]*np.sin(ang))).sum(axis=0)
    return (2*TWOPI/256.0) * m * term

def _validate():
    rng = np.random.default_rng(7)
    for trial in range(3):
        n = 40
        x = np.sort(rng.uniform(0, 256, n)); m = np.ones(n); m[rng.choice(n, 9, False)] = 2
        q = rng.normal(0, 1e-3, 255)
        g = dE_dx(x, m, q)
        E = lambda xx: float(q @ rows_of_config(xx, m))
        h = 1e-6; err = []
        for i in rng.choice(n, 6, False):
            xp = x.copy(); xp[i] += h; xm = x.copy(); xm[i] -= h
            fd = (E(xp) - E(xm)) / (2*h)
            err.append(abs(fd - g[i]) / max(1e-12, abs(fd)))
        print(f"  trial {trial}: max rel err vs finite diff = {max(err):.3e}")
        assert max(err) < 1e-5, "GRADIENT WRONG"
    print("  gradient VALIDATED")


def slide(x, m, q, mu, steps=25):
    """Descend E = sum_j q_j F(j) from an existing support atom; backtracking line search."""
    x = x.copy(); best = float(q @ rows_of_config(x, m)); eta = 1e-3
    for s in range(steps):
        g = dE_dx(x, m, q)
        gn = float(np.linalg.norm(g))
        if gn < 1e-14: break
        improved = False
        for trial in range(12):
            xn = np.mod(x - eta*g, 256.0)
            v = float(q @ rows_of_config(xn, m))
            if v < best - 1e-16:
                x, best, improved = xn, v, True
                eta *= 1.6; break
            eta *= 0.4
        if not improved: break
    return x, best

# --- D(1) constraint -------------------------------------------------------------
# D(1) = sum_j s_j - 1/2 = (E[F(256)] - 128)/65536 + O(tau).
#
# CORRECTED after reading Ceiling.lean / Signed.lean at source.  d1 = 82395317/1e8 is NOT a
# design constraint: it is the authors' own law's MEASURED D(1) = 0.82395316071283519754,
# rounded up to 8 decimals (verified: ceil(D(1)*1e8) = 82395317 exactly, and 82395316 would
# fail their checker).  In ceiling_numeric, d1 is universally quantified and enters only as
# d1*|r(1)| -- a larger value weakens the conclusion, it does not invalidate the law.  For
# band-limited windows r(+-1) = 0 kills the term outright, and ceiling_law256_signed DROPS
# it entirely given only D(1) >= 0.  Capping E[F(256)] from ABOVE was therefore
# over-constraining: it forbade exactly the cheap low-atom-count columns we want (measured
# corr(n, F(256)) = -0.94).  The genuine requirement is the SIGN condition D(1) >= 0.
F256_FLOOR = 128.0            # D(1) >= 0  <=>  E[F(256)] >= 128
F256_CACHE = {}

def _f256_of(col):
    key = id(col)
    v = F256_CACHE.get(key)
    if v is not None: return v
    pm = col[3]
    if isinstance(pm, tuple) and len(pm) == 2 and all(isinstance(z, np.ndarray) for z in pm):
        x, m = pm
        ang = 2*np.pi*np.asarray(x, float)
        zr = float((np.cos(ang)*m).sum()); zi = float((np.sin(ang)*m).sum())
        v = zr*zr + zi*zi
    else:
        # AUDIT FIX: the F(256) row is a FLOOR (E[F(256)] >= 128), so a position-less
        # column must be bounded BELOW, not above.  65536.0 is the maximum, i.e. the
        # anti-conservative direction: it would credit such a column with satisfying
        # D(1) >= 0 on no evidence.  0.0 is the only sound enclosure without positions.
        v = 0.0                          # no positions: bound BELOW (conservative)
    F256_CACHE[key] = v
    return v

def f256_vec(cols):
    return np.array([_f256_of(c) for c in cols])

def master(cols, tau):
    A = np.array([c[0] for c in cols]).T
    cost = np.array([c[1] for c in cols]); K = len(cols)
    f256 = f256_vec(cols)
    b = np.concatenate([ROWS + tau, -(ROWS - tau), [-F256_FLOOR]])
    r = linprog(cost, A_ub=np.vstack([A, -A, -f256[None, :]]), b_ub=b,
                A_eq=np.ones((1, K)), b_eq=[1.0], bounds=(0, None), method="highs",
                options={"dual_feasibility_tolerance": 1e-10,
                         "primal_feasibility_tolerance": 1e-10})
    return r, A, cost

def run(src, out, tau=1e-4, iters=40, wmin=1e-11):
    cols = pickle.load(open(src, "rb"))
    print(f"loaded {len(cols)} columns from {src}", flush=True)
    t0 = time.time(); hist = []
    for it in range(iters):
        r, A, cost = master(cols, tau)
        if not r.success:
            print("master infeasible:", r.message, flush=True); break
        q, mu, rep = extract(r, A, cost)
        obj = float(r.fun); hist.append(obj)
        sup = np.where(r.x > wmin)[0]
        added = 0; bestrc = 0.0
        for ci in sup:
            entry = cols[ci]
            pm = entry[3]
            if not (isinstance(pm, tuple) and len(pm) == 2
                    and all(isinstance(z, np.ndarray) for z in pm)):
                continue                      # provenance placeholder, no stored positions
            x, m = pm
            xn, E = slide(np.asarray(x, float), np.asarray(m, float), q, mu)
            rc = sigma_of(m) + E + mu
            if rc < -1e-12:
                cols.append((rows_of_config(xn, m), sigma_of(m), "sfw", (xn, m)))
                added += 1; bestrc = min(bestrc, rc)
        print(f"  [{it:3d}] obj={obj:.9f}  support={len(sup)}  slid+added={added}  "
              f"best_rc={bestrc:+.3e}  cols={len(cols)}  {time.time()-t0:6.1f}s", flush=True)
        if added == 0:
            print("  no support atom can be improved by sliding", flush=True); break
        pickle.dump(cols, open(out, "wb"))
    r, A, cost = master(cols, tau)
    print(f"FINAL obj={float(r.fun):.9f}  (start {hist[0]:.9f}, delta {float(r.fun)-hist[0]:+.3e})", flush=True)
    pickle.dump(cols, open(out, "wb"))

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        _validate()
    else:
        src = sys.argv[1] if len(sys.argv) > 1 else "cg_cols_final.pkl"
        out = sys.argv[2] if len(sys.argv) > 2 else "cg_cols_sfw.pkl"
        tau = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
        its = int(sys.argv[4]) if len(sys.argv) > 4 else 40
        run(src, out, tau, its)
