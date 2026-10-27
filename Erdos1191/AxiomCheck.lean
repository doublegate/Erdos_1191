/-
Copyright 2026 Lucas Parobek.
Licensed under the Apache License, Version 2.0.

Kernel axiom audit for all theorems in `Erdos1191.Proofs`.
Running `lean Erdos1191/AxiomCheck.lean` or `lake build` verifies that no theorem
depends on `sorryAx` or any custom axiom.
-/
import Erdos1191.Proofs

#print axioms Erdos1191.erdos_1191.variants.harmonic_conservation_numerator
#print axioms Erdos1191.erdos_1191.variants.two_shell_difference_localization
#print axioms Erdos1191.erdos_1191.variants.discrete_legendre_fenchel_excess
#print axioms Erdos1191.erdos_1191.variants.legendre_fenchel_squared_identity
#print axioms Erdos1191.erdos_1191.variants.robin_hood_pair_smoothing
#print axioms Erdos1191.erdos_1191.variants.scale_invariant_block_window_exponents
#print axioms Erdos1191.erdos_1191.variants.scale_invariant_block_window_product
#print axioms Erdos1191.erdos_1191.variants.hyperfactorial_sawtooth_abel_step
#print axioms Erdos1191.erdos_1191.variants.diophantine_haar_octave_divergence
#print axioms Erdos1191.erdos_1191.variants.bitshift_sidon_difference_index
#print axioms Erdos1191.maximinWitnessN32_is_sidon
#print axioms Erdos1191.mianChowlaPrefix10_is_sidon
#print axioms Erdos1191.erdos_1191.variants.cesaro_fejer_order2_step
#print axioms Erdos1191.erdos_1191.variants.zero_liminf_excludes_pos
