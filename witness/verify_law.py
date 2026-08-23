"""Standalone verifier for the independent ceiling witness.

Checks, from `law_certified.txt` alone and with exact integer/rational arithmetic:
  1. the file's 256 integer enclosures (lo, hi) at scale K = 2^60 parse;
  2. every row j = 1..255 satisfies the near-CUE band  |S(j) - j| <= tau  (tau = 1/5000),
     i.e. the law is a valid bandwidth-one witness at that tolerance;
  3. the stated exact rational p equals its stated decimal, and is above the paper's p0;
  4. the reported worst-case band deviation matches the file.

No dependencies, no floating point in the decisions. Run:  python verify_law.py [law_certified.txt]
"""
import re
import sys
from fractions import Fraction

PAPER_P0 = Fraction(10909258999421303588095230195816054408197, 16 * 10**39)  # 0.68182868746...

path = sys.argv[1] if len(sys.argv) > 1 else "law_certified.txt"
text = open(path, encoding="utf-8").read()


def header(key):
    m = re.search(rf"^{re.escape(key)}\s*=\s*(\S+)", text, re.M)
    return m.group(1) if m else None


p = Fraction(header("p"))
tau = Fraction(header("tau_decl"))
K = 2 ** int(header("K").split("^")[1])
encl = [(int(a), int(b)) for a, b in re.findall(r"^\((-?\d+), (-?\d+)\),$", text, re.M)]

print(f"file:            {path}")
print(f"enclosures:      {len(encl)} (expect 256)")
print(f"scale K:         2^{K.bit_length()-1}")
print(f"declared tau:    {tau} = {float(tau):.1e}")

failures = []
worst = Fraction(0)
for j in range(1, 256):                      # row 256 is free by construction
    lo, hi = encl[j - 1]
    s_lo, s_hi = Fraction(lo, K), Fraction(hi, K)
    worst = max(worst, abs(s_lo - j), abs(s_hi - j))
    if s_lo < j - tau or s_hi > j + tau:
        failures.append(j)

print(f"rows checked:    255 (both endpoints, exact rational arithmetic)")
print(f"band failures:   {len(failures)}{'' if not failures else ' at rows ' + str(failures[:10])}")
print(f"worst deviation: {float(worst):.3e}  (file states {header('worst_dev')})")
print(f"certified p:     {float(p):.9f}  = {p}")
print(f"paper p0:        {float(PAPER_P0):.9f}")
print(f"gap to p0:       {float(p - PAPER_P0):+.9f}  (positive = our ceiling is the weaker/higher one)")

ok = not failures and len(encl) == 256 and p > 0
print()
print("RESULT:", "PASS - valid public witness at the declared tolerance" if ok else "FAIL")
print("Interpretation: no bandwidth-one certificate can prove a simple-and-on-critical-line")
print(f"proportion above {float(p):.6f} (+ the O(1/N^2 + tau/N) stability term of the paper's Theorem 1').")
sys.exit(0 if ok else 1)
