# Independent verification of *zeta-23*, and a public witness for its bandwidth-one ceiling

This repository contains an independent third-party verification of the results of

> L. Alpöge and R. Furman, *More than two thirds of the zeta zeros are simple and on the critical line*, arXiv:2608.13637v1 (submitted 13 August 2026; a v2 dated 19 August 2026 adds provenance and extension material without altering the main theorems — this report audits v1), with Lean 4 formalization at [`anthropics/zeta-23-lean`](https://github.com/anthropics/zeta-23-lean),

together with an independently constructed, fully public **witness law** for the paper's bandwidth-one ceiling.

Everything here is reproducible from public materials. We introduce no new claims about the Riemann zeta function.

**Prior work.** We are not the first to check parts of this. An independent kernel rebuild and axiom audit
was published by [`stevemoraco/qs` PR #363](https://github.com/stevemoraco/qs/pull/363) on 13 August 2026,
and an independent reproduction of the ceiling data from the published enclosures by
[`teal-sea/zeta-lab` PR #12](https://github.com/teal-sea/zeta-lab/pull/12) on 11 August 2026 — the latter
obtaining the same digits we did. The `EnclOK` availability gap was raised in repository issues
[#8](https://github.com/anthropics/zeta-23-lean/issues/8) and
[#14](https://github.com/anthropics/zeta-23-lean/issues/14), both still open. Findings 1, 3 and 4 below
therefore corroborate work others published first. What is new here is finding 2 — a statement-fidelity
audit, which no other effort has attempted — and finding 5, an independently constructed witness.

**Author:** Bryan Carson, Hydra Dynamix — <bcarson@hydradynamix.com>
**Date:** 21 August 2026

## Summary of findings

| # | Finding | Status |
|---|---|---|
| 1 | The Lean formalization builds from source on independent hardware and all **34** headline theorems depend on exactly the three standard axioms, sorry-free | **Confirmed** |
| 2 | All **15** trusted definitions and **17** challenge statements faithfully encode the paper's informal claims | **Confirmed, 0 discrepancies** |
| 3 | The numerical certificate behind the ~0.682 ceiling (`cert_N256_blk_b128m.json`) is **not public**, so `EnclOK` cannot be verified by anyone outside the authors | **Data-availability gap** |
| 4 | The integer certificate layer downstream of `EnclOK` (255 row inequalities, edge bound, p₀ arithmetic) reproduces exactly under an independent implementation | **Confirmed** |
| 5 | An independently constructed witness law certifies a bandwidth-one ceiling of **p = 0.682434924**, within **6.1 × 10⁻⁴** of the paper's p₀ = 0.681828687 — with every link publicly reproducible | **New artifact** |

Findings 1, 2 and 4 support the paper (1 and 4 independently corroborating earlier third-party checks; see *Prior work* above). Finding 3 concerns reproducibility of one auxiliary section, not the correctness of the main theorems — **Theorems A–E are provably independent of the certificate in question**, as our axiom audit confirms. The repository's README and the paper both document the certificate's status accurately; the gap is availability, not disclosure.


## Why finding 5 matters

The paper's §7.2 shows that no *bandwidth-one certificate* — the entire method class its main theorem belongs to — can prove more than ≈ 0.6819 of zeros simple and on the critical line. That upper bound is established by exhibiting a witness law: a configuration law matching the pair-correlation data zeta is known to have, but containing only p₀ = 68.18% simple zeros.

That witness rests on a certificate file which is not distributed (finding 3). Our witness is different: it is built from scratch, its positions and weights are in this repository, and its form-factor enclosures are certified by interval arithmetic that anyone can re-run. It is very slightly weaker (0.682435 vs 0.681829 — a higher ceiling is a less informative one), but it is, to our knowledge, **the only independently constructed witness of this kind**: the prior
reproduction cited above re-checked the authors' own published enclosures, whereas this law was built from
scratch and certified without reference to their certificate.

## Layout

```
audit/
  build-summary.txt              build + axiom-audit transcript (finding 1)
  statement-fidelity-audit.md    per-definition and per-statement verdicts (finding 2)
  paper-reference-extraction.md  the paper-side reference the audit was diffed against
  independent_rowcheck.py        independent re-implementation of the RowCert checker (finding 4)
  rowcheck-output.txt            its output
witness/
  law_certified.txt              the certified witness: exact rational p, integer enclosures, support
  verify_law.py                  standalone exact-arithmetic verifier for the above (no dependencies)
  colgen.py colgen2.py           column generation (stage 1: positions; stage 2: mark moves)
  colgen3.py colgen4.py          faster pricers (FFT kernel; global-field coordinate descent)
  certify_cg.py                  dyadic snap + rigorous interval-enclosure certification
  extract_dual.py                extracts the optimal certificate (LP dual) from a snapshot
NOTE.md                          the technical note: methods, results, limitations
```

## Reproducing the audit (findings 1, 2, 4)

```sh
# 1. build the formalization and audit its axioms (~30 min with the Mathlib cache)
git clone https://github.com/anthropics/zeta-23-lean && cd zeta-23-lean
lake exe cache get && lake build && lake build Solution Solution.XiPrime
lake env lean comparator/PrintAxioms.lean            # expect: the three standard axioms only
lake env lean comparator/PrintAxioms/XiPrime.lean
lake env lean comparator/PrintAxioms/PairCeiling.lean

# 2. re-run the independent integer-certificate check (seconds, no dependencies)
python audit/independent_rowcheck.py path/to/zeta-23-lean/Zeta23/PairCeiling/LawN256.lean
```

Statement fidelity (finding 2) is a reading task: `audit/statement-fidelity-audit.md` gives per-item verdicts for
`comparator/ChallengeDeps.lean` (15 definitions) and `comparator/Challenge.lean` (17 statements) against the paper.

## Reproducing the witness (finding 5)

**Verifying the witness takes seconds and needs nothing but Python** — the certified law is
self-contained, so a skeptical reader never has to re-run the search:

```sh
cd witness && python verify_law.py        # exact rational arithmetic, no dependencies
# -> 255 rows checked, 0 band failures, certified p = 0.682434924, gap to p0 = +0.000606237
```

Regenerating the law from scratch (hours, stochastic):

```sh
pip install numpy scipy mpmath
cd witness
python colgen.py 150            # stage 1: column generation over free positions
python colgen2.py 200           # stage 2: adds mark merge/split moves
python colgen4.py 80 0.0001     # global-field pricer at tau = 1e-4
python certify_cg.py cg_cols.pkl law_certified.txt 200 0.0001   # rigorous certification
```

Column generation is stochastic; re-runs converge to nearby values rather than the identical law.
`witness/law_certified.txt` is the exact artifact our run produced, and its certification is deterministic:
the enclosures and the exact rational p can be re-verified from that file alone.

## Provenance

The verification was conducted as a pre-registered research programme in
[LDGR](https://ldgr.run) (`ldgr-research`), an evidence ledger for bounded, falsifiable experiments.
Every experiment in this report was registered with its pass/fail criteria *before* it ran, and every
artifact was checksummed at creation. (Later experiments in the programme additionally freeze a
cryptographic hash of the registered definition so that post-hoc edits are detectable; the experiments
behind this report predate that feature and are pre-registered without it.) The
findings above correspond to experiments `zeta23-exp1` (build and axioms), `zeta23-exp2` (statement
fidelity), `zeta23-exp3` (the `EnclOK` gap and the independent integer check), and `zeta23-exp7` (the
witness). Ledger records are available on request.

## Disclosure

This verification was carried out with substantial AI assistance (Claude, Anthropic), including the
column-generation search, the certification pipeline, and the drafting of this report. All quantitative
claims are machine-checkable from the artifacts here, and all reproduction commands in this README were
executed. Adversarial review of our own intermediate conclusions is described in `NOTE.md` §6.

Per the [Leiden Declaration on AI and Mathematics](https://leidendeclaration.ai/) (June 2026), we state
this involvement explicitly rather than in a footnote.

## License

Verification artifacts and code: Apache-2.0, matching the licence of the audited repository.
