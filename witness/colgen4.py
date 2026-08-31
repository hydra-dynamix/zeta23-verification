"""Pricer: field-based global coordinate descent.

For config (x, m) and kernel K(d) = sum_j q_j cos(2 pi j d/256):
  field phi(t) = sum_i m_i K(t - x_i) on the full grid via one inverse FFT of the
  config's Fourier coefficients (A_j - i B_j) weighted by q_j.
  Coordinate move: phi_minus_i = phi - m_i K(t - x_i)  (table shift), then relocate
  x_i to the global minimizer of 2 m_i phi_minus_i(t) over the grid. Marks via the
  same field: merge two simples at field-hot spots, split doubles at field-cold ones.
Usage: python colgen4.py <rounds> [tau]
"""
import sys, time, pickle
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog
from colgen import N, ROWS, P0f, rows_of_config, sigma_of, random_config, extremal_config, rng
from colgen3 import FastKernel, G

def config_spectrum(x, m):
    ang = 2 * np.pi * np.outer(np.arange(1, 256), x) / 256.0
    A = (np.cos(ang) * m).sum(axis=1)
    B = (np.sin(ang) * m).sum(axis=1)
    return A, B

def field_on_grid(q, A, B):
    """phi[g] = sum_j q_j (A_j cos(2pi j g/G) + B_j sin(2pi j g/G))"""
    spec = np.zeros(G, dtype=complex)
    js = np.arange(1, 256)
    spec[js] = q * (A + 1j * B) / 2.0
    spec[G - js] = q * (A - 1j * B) / 2.0
    return np.fft.fft(spec).real

def energy_from_spec(q, A, B, m):
    # E = sum_{i,i'} m m' K(xi - xi') = sum_j q_j (A_j^2 + B_j^2)
    return float((q * (A * A + B * B)).sum())

def global_sweep(x, m, q, ker, sweeps=6):
    """cyclic global coordinate moves + field-informed mark moves."""
    x = x.copy(); m = m.copy()
    A, B = config_spectrum(x, m)
    for sw in range(sweeps):
        moved = 0
        phi = field_on_grid(q, A, B)
        for i in rng.permutation(len(x)):
            # remove point i from spectrum
            angi = 2 * np.pi * np.arange(1, 256) * x[i] / 256.0
            A -= m[i] * np.cos(angi); B -= m[i] * np.sin(angi)
            # phi without i: recompute from spectrum (cheap: one FFT per few moves is
            # too costly per point; use shifted-kernel subtraction instead)
            idx0 = int(np.mod(x[i], 256.0) * (G / 256.0)) % G
            phi_wo = phi - m[i] * np.roll(ker.K, idx0)
            t = int(np.argmin(phi_wo))
            xn = t * (256.0 / G)
            if abs(xn - x[i]) > 1e-9:
                moved += 1
            x[i] = xn
            angn = 2 * np.pi * np.arange(1, 256) * xn / 256.0
            A += m[i] * np.cos(angn); B += m[i] * np.sin(angn)
            phi = phi_wo + m[i] * np.roll(ker.K, t)
        # mark moves: try merge hottest adjacent simple pair; split coldest double
        E = energy_from_spec(q, A, B, m)
        obj = sigma_of(m) + E
        simples = np.where(m == 1)[0]
        if len(simples) >= 2:
            xs = np.sort(x[simples])
            # try merging the closest pair
            order = np.argsort(x[simples])
            si = simples[order]
            gaps = np.diff(np.append(x[si], x[si][0] + 256))
            a = si[int(np.argmin(gaps))]; b = si[(int(np.argmin(gaps)) + 1) % len(si)]
            mid = np.mod(x[a] + np.mod(x[b] - x[a], 256.0) / 2, 256.0)
            keep = np.ones(len(x), bool); keep[[a, b]] = False
            x2 = np.append(x[keep], mid); m2 = np.append(m[keep], 2.0)
            A2, B2 = config_spectrum(x2, m2)
            if sigma_of(m2) + energy_from_spec(q, A2, B2, m2) < obj:
                x, m, A, B = x2, m2, A2, B2
                continue
        doubles = np.where(m == 2)[0]
        if len(doubles) >= 1:
            d = doubles[int(rng.integers(len(doubles)))]
            keep = np.ones(len(x), bool); keep[d] = False
            x2 = np.concatenate([x[keep], [np.mod(x[d] - 0.25, 256), np.mod(x[d] + 0.25, 256)]])
            m2 = np.concatenate([m[keep], [1.0, 1.0]])
            A2, B2 = config_spectrum(x2, m2)
            if sigma_of(m2) + energy_from_spec(q, A2, B2, m2) < obj:
                x, m, A, B = x2, m2, A2, B2
                continue
        if moved == 0:
            break
    return x, m, energy_from_spec(q, A, B, m)

_ws = {}
def _init_w(qvec, mu):
    _ws["q"] = qvec
    _ws["ker"] = FastKernel(qvec)
    _ws["mu"] = mu

def _price_one(args):
    x, m = args
    x2, m2, e = global_sweep(x, m, _ws["q"], _ws["ker"])
    return (sigma_of(m2) + e + _ws["mu"], x2, m2)

def main():
    t0 = time.time()
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    tau = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0001
    cols = pickle.load(open("cg_cols.pkl", "rb"))
    print(f"resumed {len(cols)} cols; tau = {tau}; pricer = global-field", flush=True)
    b_ub = np.concatenate([ROWS + tau, -(ROWS - tau)])
    q_smooth = None
    r = None
    for rnd in range(max_rounds):
        Am = np.array([c[0] for c in cols]).T
        cm = np.array([c[1] for c in cols])
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
        q_smooth = qvec if q_smooth is None else 0.6 * qvec + 0.4 * q_smooth
        wi = np.argsort(-r.x)
        cands = []
        for i in wi[:20]:
            if cols[i][2] == "cg":
                x, m = cols[i][3]
                cands.append((np.mod(x + rng.normal(0, 0.05, len(x)), 256.0), m.copy()))
        for k in range(44):
            nd = int(rng.integers(24, 60))
            cands.append(extremal_config(nd) if k % 2 == 0 else random_config(nd))
        with Pool(processes=8, initializer=_init_w, initargs=(q_smooth, mu)) as pool:
            found = pool.map(_price_one, cands)
        found.sort(key=lambda t: t[0])
        neg = [f for f in found if f[0] < -1e-9]
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
