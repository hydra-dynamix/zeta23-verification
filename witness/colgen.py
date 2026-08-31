"""Column generation over unrestricted positions.

Config = (positions x_i in [0,256), marks m_i in {1,2}), sum m_i = 256.
Rows F(j) = |sum_i m_i e(j x_i/256)|^2, j=1..255. No rationality constraint during
search — final certification is by rigorous interval enclosures (their own method).

Master:  min sigma.w  s.t.  j-tau <= (F w)(j) <= j+tau,  sum w = 1,  w >= 0.
Pricing: reduced cost rc(c) = sigma_c + sum_j q_j F_c(j) + mu, with q from HiGHS duals
(sign convention validated empirically each round). Minimize rc by gradient descent on
positions (table-lookup kernel) + double-count moves, batched restarts.
"""
import sys, time, pickle
import numpy as np
from scipy.optimize import linprog

N = 256
ROWS = np.arange(1, 256)
P0f = 0.6818286874638315
TAU = 0.01
rng = np.random.default_rng(23)

def rows_of_config(x, m):
    """F(j) for j=1..255 (float). x: positions array, m: marks array."""
    ang = 2 * np.pi * np.outer(ROWS, x) / 256.0        # (255, n)
    zr = (np.cos(ang) * m).sum(axis=1)
    zi = (np.sin(ang) * m).sum(axis=1)
    return zr * zr + zi * zi

def sigma_of(m):
    return float((m == 1).sum()) / N

def random_config(n_doubles):
    s = 256 - 2 * n_doubles
    n = s + n_doubles
    x = np.sort(rng.uniform(0, 256, n))
    m = np.ones(n); m[rng.choice(n, n_doubles, replace=False)] = 2
    return x, m

def extremal_config(n_doubles, jitter=0.5):
    """paper's extremal shape: rigid equispaced simples + spread doubles."""
    s = 256 - 2 * n_doubles
    n = s + n_doubles
    x = np.sort((np.arange(n) * (256.0 / n) + rng.uniform(-jitter, jitter, n)) % 256)
    m = np.ones(n)
    m[rng.choice(n, n_doubles, replace=False)] = 2
    return x, m

class Kernel:
    """table lookup for K(d) = sum_j q_j cos(2 pi j d/256) and K'."""
    G = 1 << 17
    def __init__(self, q):
        # build via FFT-like direct evaluation on grid
        d = np.arange(self.G) * (256.0 / self.G)
        ang = 2 * np.pi * np.outer(ROWS, d) / 256.0     # 255 x G — 33M floats, ok
        self.K = (q[:, None] * np.cos(ang)).sum(axis=0)
        self.Kp = (-(q * ROWS * 2 * np.pi / 256.0)[:, None] * np.sin(ang)).sum(axis=0)
        self.Q0 = q.sum()
    def k(self, d):
        idx = np.mod(d, 256.0) * (self.G / 256.0)
        return self.K[idx.astype(np.int64) % self.G]
    def kp(self, d):
        idx = np.mod(d, 256.0) * (self.G / 256.0)
        return self.Kp[idx.astype(np.int64) % self.G]

def energy(x, m, ker):
    D = x[:, None] - x[None, :]
    MM = m[:, None] * m[None, :]
    return float((MM * ker.k(D)).sum())

def grad(x, m, ker):
    D = x[:, None] - x[None, :]
    MM = m[:, None] * m[None, :]
    Kp = ker.kp(D)
    np.fill_diagonal(Kp, 0.0)
    return 2.0 * (MM * Kp).sum(axis=1)

def local_opt(x, m, ker, steps=120, lr0=2.0):
    lr = lr0
    e = energy(x, m, ker)
    for t in range(steps):
        g = grad(x, m, ker)
        gn = np.abs(g).max() + 1e-12
        x2 = np.mod(x - lr * g / gn, 256.0)
        e2 = energy(x2, m, ker)
        if e2 < e:
            x, e = x2, e2
            lr *= 1.15
        else:
            lr *= 0.5
            if lr < 1e-4: break
    return x, e

def price_batch(qvec, mu, n_restarts=24, doubles_range=(28, 52), seeds=None):
    """return list of (rc, x, m) with rc < 0 found by local search. qvec scaled for j rows."""
    ker = Kernel(qvec)
    out = []
    cands = []
    for r in range(n_restarts):
        nd = int(rng.integers(*doubles_range))
        cands.append(extremal_config(nd) if r % 2 == 0 else random_config(nd))
    if seeds:
        cands += seeds
    for x, m in cands:
        x, e = local_opt(x.copy(), m.copy(), ker)
        rc = sigma_of(m) + e + mu
        out.append((rc, x, m))
    out.sort(key=lambda t: t[0])
    return out

def main():
    t0 = time.time()
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    # ---- initial columns: previous v4 support + structured seeds ----
    cols = []      # list of (F float array 255, sigma, tag, payload)
    rows_prev = np.load("rows_v4.npy")
    meta_prev, sig_prev = pickle.load(open("fam_v4.pkl", "rb"))
    A = rows_prev[:, 1:256].astype(float).T
    c = np.array(sig_prev, dtype=float) / N
    b_ub = np.concatenate([ROWS + TAU, -(ROWS - TAU)])
    res = linprog(c, A_ub=np.vstack([A, -A]), b_ub=b_ub,
                  A_eq=np.ones((1, len(sig_prev))), b_eq=[1.0], bounds=(0, None), method="highs")
    supp = [i for i in range(len(sig_prev)) if res.x[i] > 1e-12]
    for i in supp:
        cols.append((rows_prev[i, 1:256].astype(float), sig_prev[i] / N, "v4", meta_prev[i]))
    print(f"warm start: {len(cols)} cols from v4 support (p={res.fun:.6f}) [{time.time()-t0:.0f}s]", flush=True)

    best_p = res.fun
    for rnd in range(max_rounds):
        Am = np.array([col[0] for col in cols]).T          # 255 x K
        cm = np.array([col[1] for col in cols])
        K_ = len(cols)
        r = linprog(cm, A_ub=np.vstack([Am, -Am]), b_ub=b_ub,
                    A_eq=np.ones((1, K_)), b_eq=[1.0], bounds=(0, None), method="highs")
        if r.status != 0:
            print(f"round {rnd}: master infeasible?!"); return 2
        mu = float(r.eqlin.marginals[0])
        marg = r.ineqlin.marginals                          # <= 0
        qvec = (marg[:255] - marg[255:])                    # candidate net dual, sign checked below
        # empirical sign check: rc of an existing basic column should be ~0
        i0 = int(np.argmax(r.x))
        rc0 = cm[i0] + float(qvec @ Am[:, i0]) + mu
        if abs(rc0) > 1e-6:
            qvec = -qvec
            rc0b = cm[i0] + float(qvec @ Am[:, i0]) + mu
            if abs(rc0b) > 1e-6:
                mu = -mu
                rc0c = cm[i0] + float(qvec @ Am[:, i0]) + mu
                if abs(rc0c) > 1e-6:
                    qvec = -qvec
        found = price_batch(qvec, mu, seeds=None)
        neg = [(rc, x, m) for rc, x, m in found if rc < -1e-7]
        for rc, x, m in neg[:8]:
            cols.append((rows_of_config(x, m), sigma_of(m), "cg", (x.copy(), m.copy())))
        print(f"round {rnd}: p = {r.fun:.9f} (gap {r.fun-P0f:+.5f})  cols {K_}  "
              f"pricing best rc {found[0][0]:+.2e}  added {min(len(neg),8)} [{time.time()-t0:.0f}s]", flush=True)
        best_p = min(best_p, r.fun)
        if not neg:
            print("pricing dry — CG converged at this restart budget", flush=True)
            break
    # save state for certification stage
    pickle.dump(cols, open("cg_cols.pkl", "wb"))
    print(f"FINAL master p = {r.fun:.9f}  gap {r.fun-P0f:+.5f}  cols {len(cols)} [{time.time()-t0:.0f}s]", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
