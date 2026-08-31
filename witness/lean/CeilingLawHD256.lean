/-
Zeta23/PairCeiling/CeilingLawHD256.lean -- the bandwidth-one ceiling at the independently
constructed N = 256 law.  This mirrors CeilingLaw256.lean line for line, substituting LawHD256
for LawN256; the mathematics is entirely theirs (ceiling_nearCUE), only the witness is ours.

The resulting ceiling is STRICTLY STRONGER: the same statement with p = 0.681810781927646...
in place of p0 = 0.681828687463832..., and with the tau-dependent error constant
1/(6*256^2) + tau/(2*256) evaluated at tau = 1/10^7 rather than 3/10^40, which changes it from
2.5431315e-6 to 2.5433268e-6 -- a difference of 1.95e-10, four orders below the 1.79e-5 by which
the witness value improves.
-/
import Zeta23.PairCeiling.CeilingLaw256
import Zeta23.PairCeiling.LawHD256

open scoped BigOperators Interval
open MeasureTheory Set

namespace Zeta23
namespace PairCeiling

/-- the row certificate and edge bound satisfied by the independent law. -/
theorem lawHD256_rows (S : ℕ → ℝ) (hS : EnclOK LawHD256.K S 0 LawHD256.encl) :
    NearCUE S 256 ((1 : ℝ) / 10 ^ 7) ∧ |Dfun (massOf S 256) 256 1| ≤ (82395317 : ℝ) / 10 ^ 8 := by
  obtain ⟨_, h1, h2⟩ := cert_of_checkRows LawHD256 LawHD256_check S hS
  have eN : LawHD256.N = 256 := rfl
  have etn : LawHD256.tn = 1 := rfl
  have etd : LawHD256.td = 10 ^ 7 := rfl
  have edn : LawHD256.dn = 82395317 := rfl
  have edd : LawHD256.dd = 10 ^ 8 := rfl
  simp only [eN, etn, etd, edn, edd, Nat.cast_one, Nat.cast_ofNat, Nat.cast_pow] at h1 h2
  exact ⟨h1, h2⟩

/-- **CEILING AT THE INDEPENDENT N = 256 LAW (exact constants).** -/
theorem ceiling_lawHD256 (S : ℕ → ℝ) (hS : EnclOK LawHD256.K S 0 LawHD256.encl)
    {r g h : ℝ → ℝ} {T : Set ℝ} (hT : T.Countable) {c₀ p : ℝ}
    (hr : ∀ x ∈ Icc (0:ℝ) 1, HasDerivAt r (g x) x) (hg : ContinuousOn g (Icc (0:ℝ) 1))
    (hgh : ∀ x ∈ Ioo (0:ℝ) 1 \ T, HasDerivAt g (h x) x) (hh : IntervalIntegrable h volume 0 1)
    (hvalid : c₀ + ∑ j ∈ Finset.Icc 1 256, massOf S 256 j * r ((j:ℝ)/256) ≤ p) :
    c₀ + ∫ x in (0:ℝ)..1, r x * x
      ≤ p + (82395317 : ℝ) / 10 ^ 8 * |r 1|
          + (1 / (6 * (256 : ℝ) ^ 2) + ((1 : ℝ) / 10 ^ 7) / (2 * 256)) * |g 1|
          + (1 / (6 * (256 : ℝ) ^ 2) + ((1 : ℝ) / 10 ^ 7) / (2 * 256)) * ∫ x in (0:ℝ)..1, |h x| := by
  obtain ⟨hrows, hD1⟩ := lawHD256_rows S hS
  have := ceiling_nearCUE (N := 256) (by norm_num) S (by positivity) hrows hD1 hT hr hg hgh hh hvalid
  push_cast at this
  exact this

end PairCeiling
end Zeta23
