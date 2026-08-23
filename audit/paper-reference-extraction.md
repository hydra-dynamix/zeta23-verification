# Paper reference extraction — arXiv:2608.13637 (Alpöge–Furman)

*Extracted 2026-08-19 from the local PDF via pdftotext + reconstruction. Math glyphs (√, ≤, ∧, ∩, ∑ᶠ, stars) were dropped by extraction and reconstructed; items marked [AMBIGUOUS] must be diffed against the repo/arXiv source, not against this reconstruction. This file is the paper-side input to the exp2 statement-fidelity audit.*

## Lettering note (critical for the audit)

The paper letters only **Theorem A** (zeta) and **Theorem B** (primitive Dirichlet L). The repo's A/B/C/D/E lettering is internal: `two_thirds_on_critical_line` (distinct-on-line 2/3), `thmB0_mult` (simple-on-line 2/3 = paper Thm A(i)), `thmC0_mult` (distinct 5/6 = paper Thm A(ii)), ThmD = Montgomery–Taylor constants, ThmE/ThmDE = Dirichlet.

## Theorem A (§1.1, pp. 1–2)

As T → ∞:
(i) N₀ˢ(T,2T) ≥ (2/3 − o(1)) · N(T,2T)
(ii) N_d(T,2T) ≥ (5/6 − o(1)) · N(T,2T)

Montgomery–Taylor window λ_MT in place of indicator λ₀: constants improve to 2 − c_MT⁻¹ = 0.67250… and ½(3 − c_MT⁻¹) = 0.83625…, where c_MT⁻¹ := 1/2 + (1/√2)·cot(1/√2) ≈ 1.3274988 [AMBIGUOUS: radicals reconstructed, numerically verified]. Optimal among windows/certificates of form 2 − R(λ) [CCLM17, Cor. 14]. Ceiling over all bandwidth-one certificates ≈ 0.682 (§7.2). A fortiori N₀*, N₀, Nˢ ≥ (2/3 − o(1))N on (T,2T].
Cumulative: same for (0,T], rate O(log log T / log T) (Remark 6.1; proof = dyadic summation).
Prior records: 5/12 for N₀ˢ/N [PRZZ20]; 0.6603 for N_d/N [Wu15]. Under RH: 2/3, 5/6, 0.6725 are Montgomery/CGG/Montgomery–Taylor; 0.6792 via SDP [CGdL20].
Lower bounds only (§1.4); insensitive to o(N) off-line zeros; argument also applies to Davenport–Heilbronn/Epstein where RH-analogue fails.

## Theorem B (§1.1, p. 2)

Theorem A holds verbatim for L(s,χ), any **fixed primitive** Dirichlet character χ; all errors O_q(·); no q-uniformity claimed. Paper gives NO verbatim Lean text for it ("the statement of Theorem B is analogous") — fidelity must be checked purely against the repo (ThmE/ThmDE).

## Counting definitions (§1.1, p. 1) — for 0 ≤ T₁ < T₂

ρ = β + iγ nontrivial zero, m_ρ ≥ 1 multiplicity. Interval half-open: T₁ < γ ≤ T₂ (γ = Im ρ, upper half-plane only).

| Symbol | Definition | Weighting |
|---|---|---|
| N(T₁,T₂) | Σ_{T₁<γ≤T₂} m_ρ | with multiplicity |
| N_d | #{ρ : T₁<γ≤T₂} | set count |
| N₀* | #{ρ : …, β = 1/2} | set count on line |
| N₀ | on-line zeros with multiplicity | multiplicity |
| N₀ˢ | #{ρ : …, β = 1/2, m_ρ = 1} | set count |
| Nˢ | simple zeros anywhere in strip | set count |

Chains: N₀ˢ ≤ N₀* ≤ N₀ ≤ N; N₀ˢ ≤ Nˢ ≤ N_d ≤ N. N(T) := N(0,T).

## "Proportion" semantics

Windowed dyadic (T,2T] lower bound with o(1) loss ⇔ liminf of ratio ≥ const. Lean encodes ε–T₀ form: ∀ ε > 0, ∃ T0 : ℝ, ∀ T ≥ T0, (c − ε)·Ncount ≤ numerator — over ALL real T ≥ T0. Cumulative variants use (0,T). ξ′ results (Remark 7.1) stated directly as liminf: ≥ 0.85838 (simple-on-line) and 0.92919 (distinct), quartic window v(s)=1−(7/100)(2s)²−(51/200)(2s)⁴ → 0.86864 / 0.93432; N′ counts ξ′ zeros with multiplicity in (T₁,T₂].

## Montgomery–Taylor window (§2.2 (2.6), Lemma 5.6 (5.13))

λ₀ := 1_[−1/2,1/2]; λ_MT(s) := cos(√2·s)·1_[−1/2,1/2](s) [AMBIGUOUS: √2 reconstructed, forced by R(λ_MT) value].
R(λ) := [∫λ² + ∫∫|u−v|λ(u)λ(v)] / (∫λ)². R(λ₀) = 4/3; R(λ_MT) = 1/2 + (1/√2)cot(1/√2) = c_MT⁻¹.
Certificate yields 2 − R(λ) (simple-on-line) and ½(3 − R(λ)) (distinct): λ₀ → 2/3, 5/6; λ_MT → 0.67250…, 0.83625….
Lean cMT := √2·tan(1/√2)/(1 + (1/√2)·tan(1/√2)) ≈ 0.75333 = 1/c_MT⁻¹ [AMBIGUOUS: radicals reconstructed, numerically verified equal].

## Lean listings printed in Appendix A (p. 15) — as transcribed

```lean
def IsNontrivialZero (ρ : ℂ) : Prop := riemannZeta ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1
def zeroMult (ρ : ℂ) : ℕ := (analyticOrderAt riemannZeta ρ).toNat
def zerosIn (T1 T2 : ℝ) : Set ℂ := {ρ | IsNontrivialZero ρ ∧ T1 < ρ.im ∧ ρ.im ≤ T2}
def Ncount (T1 T2 : ℝ) : ℕ := ∑ᶠ ρ ∈ zerosIn T1 T2, zeroMult ρ
def Ndist (T1 T2 : ℝ) : ℕ := (zerosIn T1 T2).ncard
def N0 (T1 T2 : ℝ) : ℕ := ∑ᶠ ρ ∈ zerosIn T1 T2 ∩ {ρ | ρ.re = 1 / 2}, zeroMult ρ
def N0star (T1 T2 : ℝ) : ℕ := (zerosIn T1 T2 ∩ {ρ | ρ.re = 1 / 2}).ncard
def N0simple (T1 T2 : ℝ) : ℕ := (zerosIn T1 T2 ∩ {ρ | ρ.re = 1/2} ∩ {ρ | zeroMult ρ = 1}).ncard
def Nsimple (T1 T2 : ℝ) : ℕ := (zerosIn T1 T2 ∩ {ρ | zeroMult ρ = 1}).ncard
```

Theorems (ε–T₀ form): two_thirds_on_critical_line[_cumulative] with N0star; thmB0_mult[_cumulative] with N0simple; thmC0_mult[_cumulative] with Ndist; montgomery_taylor_simple_on_critical_line_mult with (2 − 1/cMT − ε); montgomery_taylor_distinct_mult with (3/2 − cMT⁻¹/2 − ε). All against Ncount denominator. [Glyphs reconstructed — diff against repo files.]

## §7.2 bandwidth-one ceiling + numerical certificate (pp. 13–14) — KEY TRUST CAVEAT

- ceiling_law256: for any bandwidth-one certificate p and window r with Fourier support in [−1,1], an inequality bounding p; display heavily garbled in extraction [AMBIGUOUS — diff against repo only].
- Hypotheses: `hvalid` (sampled form-factor inequality for the 256-periodic law in Zeta23/PairCeiling/LawN256.lean) and **`EnclOK`** (form-factor enclosures LawN256.encl).
- p₀ := 1 − a_N, exact rational, rounded up **p₀ ≤ 0.6818287**; final window bound **below 0.6819**.
- **Paper's own caveat (p. 14, verbatim in substance):** at tag v1.0 the EnclOK enclosures are certified by interval arithmetic and are **NOT checked by the Lean kernel**; kernel-checked GIVEN EnclOK are: the 255 near-CUE row inequalities |256·S(j) − j| ≤ 3×10⁻⁴⁰, the edge bound |D(1)| ≤ 0.82395317 and its sign, and the analytic stability inequality. This is the only numerical certification in the paper; Theorems A and B are independent of it.
- The averaging step E_L[p] ≤ E_L[N₀ˢ/N] = 1 − a_N is done in the paper, not in Lean.
- NOTE: the repository README states the same scoping as the paper — enclosures obtained outside Lean by interval arithmetic, everything downstream kernel-checked by `decide` — and names the certificate file as available from the authors.

## Appendix A claims to cross-check

Toolchain v4.33.0-rc2; Mathlib 51e6992efd06; tag v1.0; top-level Theorem A declarations in Zeta23/Unconditional.lean and Zeta23/FinalMult.lean, "types carry no hypotheses"; comparator challenge file comparator/Challenge/Multiplicity.lean; counting functions against Mathlib riemannZeta; some inputs ported from PrimeNumberTheoremAnd [PNT+]; von Neumann trace inequality + Sylvester inertia contributed to Mathlib; #print axioms on all headline declarations = three standard axioms, no sorry [CONFIRMED independently by our exp1 on 2026-08-19].

## Ambiguity summary (do not diff against reconstructions)

1. All √ radicals (c_MT⁻¹, λ_MT, cMT) — numerically verified reconstructions only.
2. ceiling_law256 display — unrecoverable from extraction.
3. N₀ vs N₀* star glyphs — disambiguated via Lean listing.
4. ∧/∩/∑ᶠ/quantifier glyphs in Lean listings — reconstructed; inequality strictness legible and as reported.
5. Paper letters only Theorems A and B; C/D/E is repo lettering.
6. Paper prints "≈0.682", p₀ ≤ 0.6818287, "below 0.6819", LawN256, 255 kernel-checked rows (not "0.68185"/"256 enclosures" as circulated).
7. Theorem B Lean statement not printed in paper — repo-only check.
