/-
Axiom audit for the independently constructed witness.  Mirrors comparator/PrintAxioms/PairCeiling.lean.
Expect: the three standard axioms only (propext, Classical.choice, Quot.sound) on the ceiling theorem,
and propext alone on the kernel row check -- exactly as for LawN256.
-/
import Zeta23.PairCeiling.CeilingLawHD256

#print axioms Zeta23.PairCeiling.LawHD256_check
#print axioms Zeta23.PairCeiling.lawHD256_rows
#print axioms Zeta23.PairCeiling.ceiling_lawHD256
