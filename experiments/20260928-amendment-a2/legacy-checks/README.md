# Legacy checks (step 1 of protocol.md)

Each producer was run with BLCA_A2=0 on the final code, in the environment that produced its
committed file, and its output deep-diffed against the committed file (runtime and host fields
ignored, floats compared exactly as stored). Results committed on or before 2026-09-11 were produced
with Python 3.8.20, NumPy 1.24.4, SciPy 1.10.1 (conda env "ml"); later ones with Python 3.9.6,
NumPy 2.0.2, SciPy 1.13.1.

- driver-results.jsonl: the producers run by run_jobs.py (differences = number of differing paths).
  s27's three differences are metadata keys the current script adds (slide_encoder,
  slide_encoder_is_the_primary, slides_available); every value is identical.
- s8 and s10 first differed only because their committed files had been completed by hand after
  the script ran (the BH q-values of 2026-08-18; the F3 margins and published gaps of 2026-08-17
  and 2026-09-11). Both scripts now write those fields themselves and reproduce the committed files
  with zero differences.
- s12, s14, s19, s7b: the individual .txt files here.
- The GEO chain (s24, s35 geo, s36 analyze) was checked on sysu with its system python3 (Python
  3.8.10, NumPy 1.24.2), the interpreter of the committed GEO files: zero differences for all three.
- Run by hand in the "ml" environment (commands in ../run_jobs.py, jobs s13, s22, s19b, s7abl;
  s6, s19c, s23, s23b also reran there): s6, s19b, s19c, s22, s23, s23b and s11 reproduce with zero
  differences. s19 reproduces with zero differences in "ml"; under Python 3.9.6 it differs in five
  auxiliary fields (one gate-log entry and one stacking fold's weights, fourth decimal), an
  environment effect.
- s13 (selection null, 200 permutations, 76 min): every committed field reproduces exactly. The one
  difference is a field the script gained after the committed file was written
  (A_selection_null/bars/as_frozen_in_the_protocol, read from the frozen protocol).
- s29 (five-study floors, Python 3.9.6): zero differences.
- s7_ablation_and_generalisation: differs in three fields under both interpreters. Two are the
  timing and memory of the cost measurement. The third, brca/best_published_verified_folds (0.759
  committed, 0.736 now), is a correction made to the script on 2026-09-11 (DIMAF's split files
  cover bladder only, so its breast value is not on verified folds) that was never written back to
  the committed file; Supplementary Fig. 1 recomputed it from the published values. The A2 run
  writes the corrected value. It is not an A2 effect.
- s19, s19b and s32 gained code after these checks (the tuned-ridge stacking comparator, and s19b's
  known-answer switch). Every addition runs only under BLCA_A2=1; s19b was rechecked after its edit
  (zero differences in "ml"), and s19's A2 output with the comparator equals its earlier A2 output
  in every field but the two it adds.
