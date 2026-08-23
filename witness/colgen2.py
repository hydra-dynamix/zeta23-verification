"""exp7 stage 2: CG continuation with mark-move pricing (merge/split) + seeds.

Resumes from cg_cols.pkl. Pricing now interleaves position flow with mark moves:
  merge: two mark-1 points -> one mark-2 at energy-optimal midpoint (sigma -= 2/256)
  split: one mark-2 -> two mark-1 at +/-delta                       (sigma += 2/256)
accepted when they lower sigma + energy. Seeds = perturbed high-weight support configs.
"""
import sys, time, pickle
import numpy as np
from scipy.optimize import linprog
from colgen import (N, ROWS, P0f, TAU, rows_of_config, sigma_of, random_config,
                    extremal_config, Kernel, energy, grad, rng)

def local_opt2(x, m, ker, phases=4, flow_steps=70, lr0=2.0):
    def flow(x, m, steps, lr):
        e = energy(x, m, ker)
        for _ in range(steps):
            g = grad(x, m, ker)
            gn = np.abs(g).max() + 1e-12
            x2 = np.mod(x - lr * g / gn, 256.0)
            e2 = energy(x2, m, ker)
            if e2 < e: x, e, lr = x2, e2, lr * 1.15
            else:
                lr *= 0.5
                if lr < 1e-4: break
        return x, e, lr
    lr = lr0
    x, e, lr = flow(x, m, flow_steps, lr)
    for ph in range(phases):
        improved = False
        # objective = sigma + energy; mark moves trade sigma vs energy
        # --- merges: nearest simple pairs ---
        simples = np.where(m == 1)[0]
        if len(simples) >= 2:
            xs = x[simples]
            order = np.argsort(xs)
            cand_pairs = [(simples[order[i]], simples[order[(i + 1) % len(order)]])
                          for i in range(len(order))]
            gaps = [min(abs(x[a] - x[b]), 256 - abs(x[a] - x[b])) for a, b in cand_pairs]
            best_idx = np.argsort(gaps)[:6]
            for bi in best_idx:
                a, b = cand_pairs[bi]
                mid = x[a] + (x[b] - x[a]) / 2 if abs(x[b] - x[a]) < 128 else np.mod(x[a] - (256 - abs(x[b] - x[a])) / 2, 256)
                keep = np.ones(len(x), bool); keep[[a, b]] = False
                x2 = np.append(x[keep], mid); m2 = np.append(m[keep], 2.0)
                d = (energy(x2, m2, ker) - e) + (-2.0 / N)
                if d < 0:
                    x, m, e = x2, m2, e + (d + 2.0 / N)   # e tracks pure energy
                    e = energy(x, m, ker)
                    improved = True
                    break
        # --- splits ---
        doubles = np.where(m == 2)[0]
        for di in rng.permutation(doubles)[:6]:
            for delta in (0.15, 0.4):
                keep = np.ones(len(x), bool); keep[di] = False
                x2 = np.concatenate([x[keep], [np.mod(x[di] - delta, 256), np.mod(x[di] + delta, 256)]])
                m2 = np.concatenate([m[keep], [1.0, 1.0]])
                d = (energy(x2, m2, ker) - e) + (2.0 / N)
                if d < 0:
                    x, m, e = x2, m2, energy(x2, m2, ker)
                    improved = True
                    break
            else:
                continue
            break
        x, e, lr = flow(x, m, flow_steps, max(lr, 0.3))
        if not improved and ph > 0:
            break
    return x, m, e

def price2(qvec, mu, seeds, n_restarts=40, doubles_range=(24, 60)):
    ker = Kernel(qvec)
    cands = []
    for r in range(n_restarts):
        nd = int(rng.integers(*doubles_range))
        cands.append(extremal_config(nd) if r % 2 == 0 else random_config(nd))
    for (x, m) in seeds:
        xx = np.mod(x + rng.normal(0, 0.15, len(x)), 256.0)
        cands.append((xx, m.copy()))
    out = []
    for x, m in cands:
        x2, m2, e = local_opt2(x.copy(), m.copy(), ker)
        rc = sigma_of(m2) + e + mu
        out.append((rc, x2, m2))
    out.sort(key=lambda t: t[0])
    return out

def main():
    t0 = time.time()
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    cols = pickle.load(open("cg_cols.pkl", "rb"))
    print(f"resumed {len(cols)} cols", flush=True)
    b_ub = np.concatenate([ROWS + TAU, -(ROWS - TAU)])
    r = None
    for rnd in range(max_rounds):
        Am = np.array([col[0] for col in cols]).T
        cm = np.array([col[1] for col in cols])
        K_ = len(cols)
        r = linprog(cm, A_ub=np.vstack([Am, -Am]), b_ub=b_ub,
                    A_eq=np.ones((1, K_)), b_eq=[1.0], bounds=(0, None), method="highs")
        if r.status != 0: print("master infeasible"); return 2
        mu = float(r.eqlin.marginals[0])
        marg = r.ineqlin.marginals
        qvec = (marg[:255] - marg[255:])
        i0 = int(np.argmax(r.x))
        if abs(cm[i0] + float(qvec @ Am[:, i0]) + mu) > 1e-6:
            qvec = -qvec
            if abs(cm[i0] + float(qvec @ Am[:, i0]) + mu) > 1e-6:
                mu = -mu
                if abs(cm[i0] + float(qvec @ Am[:, i0]) + mu) > 1e-6:
                    qvec = -qvec
        # seeds: top-weight cg columns
        wi = np.argsort(-r.x)
        seeds = []
        for i in wi[:10]:
            if cols[i][2] == "cg":
                seeds.append(cols[i][3])
        found = price2(qvec, mu, seeds)
        neg = [(rc, x, m) for rc, x, m in found if rc < -1e-7]
        for rc, x, m in neg[:8]:
            cols.append((rows_of_config(x, m), sigma_of(m), "cg", (x.copy(), m.copy())))
        print(f"round {rnd}: p = {r.fun:.9f} (gap {r.fun-P0f:+.5f})  cols {K_}  "
              f"best rc {found[0][0]:+.2e}  added {min(len(neg),8)} [{time.time()-t0:.0f}s]", flush=True)
        if not neg:
            print("pricing dry", flush=True)
            break
        if rnd % 10 == 0:
            pickle.dump(cols, open("cg_cols.pkl", "wb"))
    pickle.dump(cols, open("cg_cols.pkl", "wb"))
    print(f"FINAL p = {r.fun:.9f}  gap {r.fun-P0f:+.5f}  cols {len(cols)} [{time.time()-t0:.0f}s]", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
