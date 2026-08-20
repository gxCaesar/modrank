# Figures

Eight figures, five main and three supplementary. Every one is built by the script beside it, at its
final size, and loads its numbers from a results file rather than carrying them in the script.

| file | builder | size (in) | in the paper |
|---|---|---|---|
| `fig1_overview.pdf` | `fig1_overview.py` | 6.30 × 8.60 | Figure 1 — the defect, the architecture, the protocol |
| `fig2_missing_variable.pdf` | `fig2_missing_variable.py` | 6.30 × 4.70 | Figure 2 — the missing variable, quantified |
| `fig3_method_result.pdf` | `fig3_method_result.py` | 6.30 × 8.90 | Figure 3 — the result, and what fails Holm |
| `fig4_regime.pdf` | `fig4_regime.py` | 6.30 × 6.90 | Figure 4 — what limits a method at n=359 |
| `fig5_biology.pdf` | `fig5_biology.py` | 6.30 × 7.60 | Figure 5 — the biology |
| `figS1_components.pdf` | `figS1_components.py` | 6.30 × 3.30 | Supplementary S1 — seven falsified components |
| `figS2_generalisation.pdf` | `figS2_generalisation.py` | 6.30 × 3.10 | Supplementary S2 — five-cohort generalisation |
| `figS3_encoders.pdf` | `figS3_encoders.py` | 6.30 × 2.90 | Supplementary S3 — seven slide encoders |

**Superseded, kept for provenance, referenced by nothing:** `fig1_clinical_baseline.{py,pdf}` and
`fig1_source_data.json` are an earlier two-panel draft of Figure 1, written before the figure set was
designed. `fig1_overview` replaces them. They are left in place rather than deleted so that the
earlier state is recoverable, but nothing in either manuscript includes them — if you are looking
for Figure 1, it is `fig1_overview.pdf`.

## The rules these are built under

- **Authored at final size.** `figsize` is the exact inserted width, insertion is at `\textwidth`,
  so the scale is 1.0000 and source pt equals rendered pt. No `\resizebox`.
- **Saved through `style.save`**, which is the only thing that actually disables the tight bounding
  box — passing `bbox_inches=None` does not, because `savefig` reads `None` as "not given".
- **No numeral is typed into a figure script.** Every quantity is loaded from a results JSON through
  a lookup that raises on a missing key, so a renamed field fails the build instead of drawing a
  stale value. `fig1_overview.py` has the explicit `dig()`; the rest index directly and fail the
  same way.
- **Asserts that can fail.** Several panels exist to show something unflattering and assert that it
  is still true: `fig4` fails if the gated component's control stops beating it, `figS1` fails if no
  control beats its component, `figS3` fails if the seven-encoder ensemble stops losing
  significantly, `fig5` fails if the case grouping drops a case.
- **Palettes are checked at import.** `style.assert_palettes()` runs three checks: the ordered
  ladder survives greyscale, the nominal set does not smuggle an ordering into an unordered
  variable, and the risk triple comes from the ladder and runs dark-to-pale monotonically.
- **sRGB.** OUP's general guidance says CMYK; Briefings in Bioinformatics has been online-only open
  access since January 2024 and states no colour rule of its own. Recorded as a conflict in
  `briefings_in_bioinformatics_lock.json`, not resolved here — converting a luminance-ordered ladder
  to CMYK moves the gaps the ladder exists for.
