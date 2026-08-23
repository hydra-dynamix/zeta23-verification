"""Independent re-implementation of zeta-23's RowCert checker (exp3, LDGR program riemann).

Re-verifies, from the enclosure integers in Zeta23/PairCeiling/LawN256.lean alone and
with exact integer arithmetic (no floats), everything the Lean `decide +kernel` call
`LawN256_check` certifies:
  (1) row inequalities  |N*lo_j - j*K| * td <= tn * K  and same for hi_j,  1 <= j < N
  (2) edge bound        |2*sum(lo) - K*N| * dd <= dn * (2*K*N)  and same for sum(hi)
plus structural facts and the p0 decimal claim from the file header.

This is an independent implementation (Python) of the checker semantics read from
RowCert.lean; it does NOT verify EnclOK itself (the law JSON cert_N256_blk_b128m.json,
sha256 cc3de991..., is not in the public repo or arXiv ancillary files).
"""
import re, sys
from fractions import Fraction

LEAN = sys.argv[1] if len(sys.argv) > 1 else "../zeta23-exp1/zeta-23-lean/Zeta23/PairCeiling/LawN256.lean"
src = open(LEAN, encoding="utf-8").read()

# ---- parse certificate data from the Lean source ----
N   = int(re.search(r"N := (\d+)", src).group(1))
K   = int(re.search(r"K := (\d+)", src).group(1))
tn  = int(re.search(r"tn := (\d+)", src).group(1))
td  = int(re.search(r"td := (\d+)", src).group(1))
dn  = int(re.search(r"dn := (\d+)", src).group(1))
dd  = int(re.search(r"dd := (\d+)", src).group(1))
encl = [(int(a), int(b)) for a, b in re.findall(r"\((\d+), (\d+)\)", src)]

print(f"parsed: N={N}  K=2^140? {K == 2**140}  tn/td = {tn}/{td}  dn/dd = {dn}/{dd}")
print(f"enclosures parsed: {len(encl)}  (checker requires len = N: {len(encl) == N})")

fail = 0

# ---- structural facts ----
widths = {hi - lo for lo, hi in encl}
print(f"enclosure widths: {widths} (all width 1 in units of 1/K = 2^-140: {widths == {1}})")
grid = all(encl[j - 1] == (j * 2**132 - 1, j * 2**132) for j in range(1, N))
print(f"rows 1..{N-1} are exactly (j*2^132 - 1, j*2^132): {grid}")

# ---- (1) row inequalities, exact integers, both endpoints ----
row_fail = []
for j in range(1, N):  # j = row index; rowsOK skips j >= N
    lo, hi = encl[j - 1]
    for name, e in (("lo", lo), ("hi", hi)):
        if abs(N * e - j * K) * td > tn * K:
            row_fail.append((j, name))
print(f"row inequalities |N*e - j*K|*td <= tn*K for j=1..{N-1}, both endpoints: "
      f"{'ALL PASS' if not row_fail else f'FAILURES: {row_fail}'}")
fail += len(row_fail)

# implied deviation bound vs tau (informative)
max_dev = max(Fraction(abs(N * e - j * K), K) for j in range(1, N) for e in encl[j - 1])
print(f"max implied |N*S(j) - j| over enclosure endpoints = {float(max_dev):.3e}  "
      f"(tau = {tn}/{td} = {tn/td:.1e}; margin factor ~{float(Fraction(tn, td)/max_dev):.2f}x)")

# ---- (2) edge bound (D(1)) ----
sum_lo = sum(lo for lo, _ in encl)
sum_hi = sum(hi for _, hi in encl)
edge_ok_lo = abs(2 * sum_lo - K * N) * dd <= dn * (2 * K * N)
edge_ok_hi = abs(2 * sum_hi - K * N) * dd <= dn * (2 * K * N)
print(f"edge bound |2*sumLo - K*N|*dd <= dn*2*K*N: {edge_ok_lo}; sumHi version: {edge_ok_hi}")
fail += (not edge_ok_lo) + (not edge_ok_hi)
D1_lo = Fraction(2 * sum_lo - K * N, 2 * K * N)   # implied D(1) = T/N - 1/2 endpoints
D1_hi = Fraction(2 * sum_hi - K * N, 2 * K * N)
print(f"implied D(1) in [{float(D1_lo):.9f}, {float(D1_hi):.9f}]  (claimed |D(1)| <= {dn/dd})")

# ---- S(256) and p0 header claims ----
S_N = Fraction(encl[N - 1][0], K)
print(f"S({N}) ~= {float(S_N):.6f} (from last enclosure; row check correctly skips j = N)")
p0 = Fraction(10909258999421303588095230195816054408197, 16 * 10**39)
print(f"p0 exact = {float(p0):.11f}; header decimal 0.68182868746... "
      f"{'OK' if str(float(p0)).startswith('0.6818286874') else 'MISMATCH'}; "
      f"paper's rounded-up bound p0 <= 0.6818287: {p0 <= Fraction(6818287, 10**7)}")

print(f"\nRESULT: {'PASS - independent checker agrees with LawN256_check' if fail == 0 else f'{fail} FAILURES'}")
print("NOT verified here (data non-public): EnclOK itself - that the true form factor S of the")
print("law in cert_N256_blk_b128m.json lies in these enclosures. Requires the JSON from the authors.")
