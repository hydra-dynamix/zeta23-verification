"""exp7 stage 5: fast CG — multiprocess pricing, FFT kernel, fat rounds, dual smoothing.
Resumes from cg_cols.pkl. Optional tight-tau mode for the endgame polish:
    python colgen3.py <rounds> [tau]
"""
import sys, time, pickle
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog
from colgen import N, ROWS, P0f, rows_of_config, sigma_of, random_config, extremal_config, rng
import colgen2

G = 1 << 17

class FastKernel:
    def __init__(self, q):
        spec = np.zeros(G, dtype=complex)
        js = np.arange(1, 256)
        spec[js * (G // 256) // 1] = 0  # placeholder; use direct irfft bins
        # K(d_g) = sum_j q_j cos(2 pi j g / G * (256/256)) with d = g*256/G:
        # angle = 2 pi j d /256 = 2 pi j g / G  -> bin j of length-G inverse FFT
        full = np.zeros(G, dtype=complex)
        full[js] = q / 2.0
        full[G - js] = q / 2.0
        self.K = np.fft.fft(full).real            # K[g] = sum_j q_j cos(2 pi j g/G)
        dspec = np.zeros(G, dtype=complex)
        coef = -(q * js * 2 * np.pi / 256.0)
        dspec[js] = 1j * coef / 2.0
        dspec[G - js] = -1j * coef / 2.0
        self.Kp = np.fft.fft(dspec).real
        self.Q0 = q.sum()
    def k(self, d):
        idx = np.mod(d, 256.0) * (G / 256.0)
        return self.K[idx.astype(np.int64) % G]
    def kp(self, d):
        idx = np.mod(d, 256.0) * (G / 256.0)
        return self.Kp[idx.astype(np.int64) % G]

_worker_state = {}
def _init_worker(qvec, mu):
    _worker_state["ker"] = FastKernel(qvec)
    _worker_state["mu"] = mu

def _price_one(args):
    x, m = args
    ker = _worker_state["ker"]
    x2, m2, e = colgen2.local_opt2(x, m, ker)
    return (sigma_of(m2) + e + _worker_state["mu"], x2, m2)

def main():
    t0 = time.time()
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    tau = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    cols = pickle.load(open("cg_cols.pkl", "rb"))
    print(f"resumed {len(cols)} cols; tau = {tau}", flush=True)
    b_ub = np.concatenate([ROWS + tau, -(ROWS - tau)])
    q_smooth = None
    r = None
    with Pool(processes=8) as dummy:
        pass  # warm import check
    for rnd in range(max_rounds):
        Am = np.array([c[0] for c in cols]).T
        cm = np.array([c[1] for c in cols])
        K_ = len(cols)
        r = linprog(cm, A_ub=np.vstack([Am, -Am]), b_ub=b_ub,
                    A_eq=np.ones((1, K_)), b_eq=[1.0], bounds=(0, None), method="highs")
        if r.status != 0:
            print(f"master infeasible at tau={tau}"); return 2
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
        q_smooth = qvec if q_smooth is None else 0.6 * qvec + 0.4 * q_smooth
        # candidates: seeds (top-weight cg cols perturbed) + fresh
        wi = np.argsort(-r.x)
        cands = []
        for i in wi[:16]:
            if cols[i][2] == "cg":
                x, m = cols[i][3]
                cands.append((np.mod(x + rng.normal(0, 0.12, len(x)), 256.0), m.copy()))
        for k in range(48):
            nd = int(rng.integers(24, 60))
            cands.append(extremal_config(nd) if k % 2 == 0 else random_config(nd))
        with Pool(processes=8, initializer=_init_worker, initargs=(q_smooth, mu)) as pool:
            found = pool.map(_price_one, cands)
        found.sort(key=lambda t: t[0])
        neg = [f for f in found if f[0] < -1e-8]
        for rc, x, m in neg[:32]:
            cols.append((rows_of_config(x, m), sigma_of(m), "cg", (x.copy(), m.copy())))
        print(f"round {rnd}: p = {r.fun:.9f} (gap {r.fun-P0f:+.5f}) cols {K_} "
              f"best rc {found[0][0]:+.2e} added {min(len(neg),32)} [{time.time()-t0:.0f}s]", flush=True)
        if not neg:
            print("pricing dry", flush=True)
            break
        if rnd % 5 == 0:
            pickle.dump(cols, open("cg_cols.pkl", "wb"))
    pickle.dump(cols, open("cg_cols.pkl", "wb"))
    print(f"FINAL p = {r.fun:.9f} gap {r.fun-P0f:+.5f} cols {len(cols)} [{time.time()-t0:.0f}s]", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
