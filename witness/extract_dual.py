"""record-exp8: extract the optimal bandwidth-one certificate (LP dual) from a CG snapshot.

The master LP:  min sigma.w  s.t.  A w <= j+tau, -A w <= -(j-tau), 1.w = 1, w >= 0
Dual feasibility for every configuration c:  sigma_c + sum_j q_j F_c(j) + mu >= 0,
i.e. the CERTIFICATE  p_cert(c) := -mu - sum_j q_j F_c(j)  satisfies p_cert(c) <= sigma_c
configuration-wise (over explored configs). Its value on near-CUE data (F(j) ~= j):
value = -mu - sum_j q_j * j  ~= master optimum p.

Paper normalization (ceiling_law256): certificate (c0, r) with
p_cert = c0 + sum_j (S_j/256) r(j/256) and S_j = F(j)/256, so
r(j/256) = -256^2 * q_j... up to sign bookkeeping; we report both raw q and the
scaled window samples R_j := -q_j * 256 * j-normalization-free form for shape analysis.
"""
import pickle, sys
import numpy as np
from scipy.optimize import linprog

N, ROWS = 256, np.arange(1, 256)
P0f = 0.6818286874638315

snap_file = sys.argv[1] if len(sys.argv) > 1 else "cg_cols.pkl"
out_file = sys.argv[2] if len(sys.argv) > 2 else "dual_certificate.txt"
tau = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0001

cols = pickle.load(open(snap_file, "rb"))
Am = np.array([c[0] for c in cols]).T
cm = np.array([c[1] for c in cols])
K = len(cols)
b_ub = np.concatenate([ROWS + tau, -(ROWS - tau)])
r = linprog(cm, A_ub=np.vstack([Am, -Am]), b_ub=b_ub, A_eq=np.ones((1, K)), b_eq=[1.0],
            bounds=(0, None), method="highs")
assert r.status == 0
mu = float(r.eqlin.marginals[0])
marg = r.ineqlin.marginals
qv = (marg[:255] - marg[255:])
# sign fix via basic column
i0 = int(np.argmax(r.x))
if abs(cm[i0] + float(qv @ Am[:, i0]) + mu) > 1e-6:
    qv = -qv
    if abs(cm[i0] + float(qv @ Am[:, i0]) + mu) > 1e-6:
        mu = -mu
        if abs(cm[i0] + float(qv @ Am[:, i0]) + mu) > 1e-6:
            qv = -qv
value = -mu - float(qv @ ROWS)
print(f"master p = {r.fun:.9f}; certificate value on ramp = {value:.9f} (match {abs(value-r.fun)<1e-6})")
print(f"c0 = -mu = {-mu:.9f}")
# window samples: paper form p = c0 + sum_j (S_j/256) r(j/256); S_j ~ F/256; our form
# p = -mu - sum q_j F_j = -mu + sum_j (F_j/256/256) * (-q_j*256^2)
Rsamp = -qv * 256.0 * 256.0          # r(j/256) samples
# Montgomery-Taylor comparison: the achieved certificate's effective r is quadratic-ish
# report shape stats
with open(out_file, "w") as f:
    f.write(f"snapshot = {snap_file}  tau = {tau}\nmaster_p = {r.fun:.9f}\n"
            f"certificate_value_on_ramp = {value:.9f}\nc0 = {-mu:.9f}\np0_ref = {P0f}\n\n")
    f.write("j, q_j, r_sample(j/256)\n")
    for j in range(1, 256):
        f.write(f"{j}, {qv[j-1]:+.6e}, {Rsamp[j-1]:+.6e}\n")
print(f"written {out_file}")
# quick shape summary
nz = np.abs(Rsamp) > np.abs(Rsamp).max() * 1e-3
print(f"support of r: {nz.sum()}/255 rows above 0.1% of max; max |r| at j = {int(ROWS[np.argmax(np.abs(Rsamp))])}")
print("r samples at j = 1,2,3,5,10,32,64,128,192,255:")
for j in (1,2,3,5,10,32,64,128,192,255):
    print(f"  r({j}/256) = {Rsamp[j-1]:+.4e}")
