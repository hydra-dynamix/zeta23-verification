# Statement-Fidelity Audit of the zeta-23 Lean Formalization

**Date:** 2026-08-19
**Auditor:** Bryan Carson (Hydra Dynamix), with AI assistance (Claude, Anthropic), independent of the zeta-23 team and its artifacts. Disclosure: the assisting model is an Anthropic product; every check below is mechanical or reproducible from the cited materials.
**Subject:** github.com/anthropics/zeta-23-lean @ commit cloned 2026-08-19 (paper: arXiv:2608.13637v1, Alpöge–Furman)
**Inputs:** comparator/ChallengeDeps.lean (109 lines), comparator/Challenge.lean (231 lines), comparator/Solution.lean, comparator/PrintAxioms*.lean, paper PDF (local copy; extraction reference in paper-reference.md)
**Companion:** exp1 (build + axiom reproduction, run 1, PASS) — this report covers WHAT is claimed; exp1 covered whether the proofs check.

## Verdict

**PASS — 0 discrepancies.** All 15 trusted definitions and all 17 challenge statements faithfully encode the paper's informal claims, subject to the declared trust base (Mathlib + Lean kernel) and the residual items in §5. The headline informal claim "≥ 67.25% of nontrivial zeta zeros are simple and on the critical line" is exactly `montgomery_taylor_simple_on_critical_line_mult`: ∀ ε > 0, eventually (2 − 1/cMT − ε)·N(T,2T) ≤ N₀ˢ(T,2T), with 2 − 1/cMT = 0.6725007… verified numerically and algebraically below.

## 1. Mechanical pre-checks (this audit, on top of exp1)

| Check | Result |
|---|---|
| Defs inlined in Challenge.lean vs ChallengeDeps.lean | **Character-identical** (mechanical diff) |
| 17 theorem statements: Solution.lean vs trusted Challenge.lean | **Byte-identical mod whitespace** (mechanical extraction + diff); so exp1's PrintAxioms audit (which imports Solution) transfers to the trusted statements |
| Count | Exactly 15 defs, exactly 17 statements, as advertised |
| Kernel-level statement agreement | Build success of Solution.lean is itself evidence: each delegation `exact Zeta23.thmX` typechecks only if the trusted root-namespace definitions unfold identically to the Zeta23-library ones in the kernel |
| cMT numeric | cMT = 0.7532960679; 1/cMT = 1.3274992963 = 1/2 + (1/√2)·cot(1/√2) (agreement < 1e-14); proportions 2 − 1/cMT = 0.6725007, (3 − 1/cMT)/2 = 0.8362504 — match the paper's 0.67250…/0.83625… |

Not run: the external leanprover/comparator tool itself (elaboration-level statement equality). The byte-identity diff plus the kernel-unfolding argument above make this a small residual (§5.2).

## 2. Definition-by-definition verdicts (15)

| # | Definition | Verdict | Notes |
|---|---|---|---|
| 1 | `IsNontrivialZero ρ := riemannZeta ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1` | **Faithful** | Open critical strip. Matches the paper's convention; that all non-trivial zeros lie in the open strip (boundary zero-freeness) is classical and not needed for the claims — both numerator and denominator use the same set, and the docstring is honest about this. |
| 2 | `zeroMult ρ := (analyticOrderAt riemannZeta ρ).toNat` | **Faithful** (note A) | Order of vanishing = multiplicity m_ρ. `.toNat` maps the ⊤ (locally-identically-zero) case to 0; vacuous since ζ is nowhere locally ≡ 0 (provable), and `zeroMult ρ = 1` in N0simple is unaffected. |
| 3 | `zerosIn T₁ T₂` | **Faithful** | Half-open window T₁ < Im ρ ≤ T₂, positive-ordinate (Im, not \|Im\|) — exactly the paper's §1.1 convention. Total in T₁,T₂ ∈ ℝ (paper restricts to 0 ≤ T₁ < T₂), but every theorem instantiates only (T,2T] with large T or (0,T] — no semantic drift. |
| 4 | `Ncount := ∑ᶠ ρ ∈ zerosIn, zeroMult ρ` | **Faithful** (note B) | Paper's N(T₁,T₂), with multiplicity. `∑ᶠ` (finsum) returns 0 on infinite support; finiteness of zeros in a horizontal strip window is classical, provable, and proved solution-side, so the definition denotes the true count. Since Ncount is the *denominator* (LHS multiplier) in all 17 statements, an infinite-support reading would have trivialized them — the reason this note matters and was checked. |
| 5 | `Ndist := (zerosIn).ncard` | **Faithful** | Paper's N_d (set count). ncard's infinite→0 convention could only shrink a *numerator*, i.e. strengthen the claim — no vacuous-truth risk. |
| 6 | `N0star := (zerosIn ∩ {re = 1/2}).ncard` | **Faithful** | Paper's N₀*: distinct on-line zeros, exact equality Re ρ = 1/2. |
| 7 | `N0simple := (zerosIn ∩ {re = 1/2} ∩ {zeroMult = 1}).ncard` | **Faithful** | Paper's N₀ˢ: simple (m_ρ = 1) and on the line. |
| 8–14 | `IsNontrivialZeroL, zeroMultL, zerosInL, NcountL, NdistL, N0starL, N0simpleL` | **Faithful** | Verbatim mirrors against Mathlib's `DirichletCharacter.LFunction χ` (analytic continuation of L(s,χ)). Same notes A/B apply. |
| 15 | `cMT := √2·tan(1/√2) / (1 + (1/√2)·tan(1/√2))` | **Faithful** | Equals the Montgomery–Taylor constant: algebraically, (1+θt)/(√2·t) = θ/t + 1/2 = 1/2 + (1/√2)cot(1/√2) at θ = 1/√2, t = tan θ; numerically verified to 1e-14. Docstring correctly states the decimals are not part of the formal statements. Solution.lean proves it equal to the library's division-safe `cStar 1` form. |

## 3. Statement-by-statement verdicts (17)

Encoding convention: "liminf_{T→∞} X(T)/N(T) ≥ c" as "∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, (c − ε)·N ≤ X". This is the standard and correct encoding: it is equivalent to the liminf claim whenever N(T) → ∞ (true: N(T,2T) ≍ T log T), and for small T where a window is empty the inequality is trivially true, matching the "as T → ∞" reading. All denominators are multiplicity-weighted (Ncount/NcountL), exactly as in the paper.

| Statement(s) | Paper claim | Verdict |
|---|---|---|
| `two_thirds_on_critical_line` (+`_cumulative`) | Thm A "a fortiori" N₀*/N ≥ 2/3, windows (T,2T] and (0,T] | **Faithful** |
| `two_thirds_simple_on_critical_line` (+`_cumulative`) | Thm A(i): N₀ˢ/N ≥ 2/3 | **Faithful** |
| `five_sixths_distinct` (+`_cumulative`) | Thm A(ii): N_d/N ≥ 5/6 | **Faithful** |
| `montgomery_taylor_on_critical_line` | MT window: N₀*/N ≥ 2 − 1/c₁* (dyadic) | **Faithful** |
| `montgomery_taylor_simple_on_critical_line_mult` (+`_cumulative`) | **The headline 67.25%**: N₀ˢ/N ≥ 2 − 1/c₁* | **Faithful** |
| `montgomery_taylor_distinct_mult` (+`_cumulative`) | N_d/N ≥ (3 − 1/c₁*)/2 = 0.83625… | **Faithful** (3/2 − cMT⁻¹/2 ≡ (3 − 1/cMT)/2) |
| `dirichlet_*` (6 statements) | Thm B: same bounds for L(s,χ), χ primitive mod q, dyadic | **Faithful** — hypotheses `1 < q` and `χ.IsPrimitive` match "fixed primitive χ"; q = 1 (χ trivial, L = ζ) is excluded and covered by the ζ theorems; constants identical; no q-uniformity claimed in paper or Lean |

**Coverage direction:** every Lean statement is a claim the paper makes (no overstatement). The paper makes a few claims *not* formalized (cumulative MT bound for N₀*, cumulative Dirichlet variants, Nˢ/N₀ a-fortiori variants) — omissions, which cannot inflate the result.

## 4. The §7.2 ceiling machinery — separate trust status (paper's own caveat, confirmed in source)

The bandwidth-one ceiling theorems (`Zeta23.PairCeiling.*`, incl. `ceiling_law256`) take the hypothesis `EnclOK LawN256.K S 0 LawN256.encl`. Confirmed directly in source (`Zeta23/PairCeiling/LawN256.lean` header): the enclosures are "**verified outside Lean by interval arithmetic**" — matching the paper's p. 14 caveat and the repository README, which states the same scoping accurately ("obtained outside Lean by interval arithmetic … available from the authors", with "everything downstream of the enclosures … checked in the kernel by `decide`"). Kernel-checked *given* EnclOK: the 255 near-CUE row inequalities, the edge bound |D(1)| ≤ 0.82395317, and the stability inequality. **Theorems A–E do not depend on any of this** — verified independently by exp1's axiom audit passing on all 17 statements. Independent recomputation of the enclosures is queued as exp3 (`zeta23-enclosure-recompute`); it is the one place where a certified claim of the paper rests on non-kernel numerics.

## 5. Residual trust base (what a reader must still accept)

1. **Mathlib semantics** (declared trust base): `riemannZeta` is the analytic continuation of ζ; `DirichletCharacter.LFunction` is L(s,χ); `analyticOrderAt` is the order of vanishing. These are mature, heavily-reviewed Mathlib definitions.
2. **Comparator tool not re-run here**: statement equality Challenge↔Solution was established by textual byte-identity plus the kernel-unfolding argument (§1), not by the leanprover/comparator elaboration check. Low residual; closable by running the tool (config present in repo).
3. **This audit's paper-side input** came from a lossy PDF text extraction with reconstructed math glyphs; every reconstructed item was resolved against the repo source or verified numerically (see paper-reference.md ambiguity ledger). The two ambiguities that could not be resolved from the PDF (`ceiling_law256` display form; exact radicals) were resolved by reading the Lean source directly, which is the ground truth being audited.
4. **Not audited**: the ξ′ comparator statements (`xiPrime_*`, 6 theorems) and PairCeiling statements (11) — exp1 confirmed their axiom hygiene but their statement fidelity vs Remark 7.1/§7.2 was out of exp2's registered scope (15 defs + 17 statements). Extendable in a follow-up if desired.

## 6. Conclusions

1. The zeta-23 headline claims are **faithfully formalized**: given Mathlib + the Lean kernel, "≥ 2/3 (indeed ≥ 0.67250…) of nontrivial zeta zeros are simple and on the critical line" (dyadic and cumulative liminf senses), the 5/6 distinctness bound, and the primitive-Dirichlet analogues are theorems whose formal statements mean what the paper says.
2. Combined with exp1 (independent build + axiom audit), **both legs of trust — "do the proofs check?" and "do the statements say the right thing?" — now hold on independent examination.**
3. The single remaining numerical-trust item in the *paper* (not in Theorems A–E) is the EnclOK enclosure certificate for the ~0.682 ceiling → exp3.
