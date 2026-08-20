# Figure design — decide here, draw after

Nothing is drawn until this document is agreed. Every panel below names the file its numbers come
from, so drawing is mechanical and no panel can quietly invent a value.

## Geometry, fixed now because authoring at final size requires it

The manuscript is `a4paper, margin=2.5cm`, so **textwidth = 16.0 cm = 6.30 in**. Every figure is
authored at exactly that width in `figsize` and inserted at `\textwidth`, giving insertion scale
1.0000 and **source pt = rendered pt**. No `\resizebox`, no `bbox_inches="tight"` — the save is
wrapped in an `rc_context` that actually disables the tight box, because passing `bbox_inches=None`
does not.

| figure | panels | size (in) | page share |
|---|---|---|---|
| 1 · the missing variable | 4 | 6.30 × 6.60 | ~2/3 |
| 2 · the method and its result | 6 | 6.30 × 8.20 | ~4/5 |
| 3 · what limits a method at n = 359 | 5 | 6.30 × 6.60 | ~2/3 |
| 4 · the biology | 6 | 6.30 × 8.20 | ~4/5 |
| S1 · seven falsified components + generalisation | 3 | 6.30 × 4.20 | supplement |

## Colour, and the one rule that governs it

**Lightness is an ordered channel, so it is spent only on ordered quantities.** Two palettes, never
mixed — a pale colour must not mean "a category" in one panel and "high on the scale" in the next.

**Nominal** (arm, method, modality — unordered), chroma fixed, hues spread:

| role | hex |
|---|---|
| **ours** | `#6D98D3` |
| clinical / stage | `#89A467` |
| slide (morphology) | `#C98E66` |
| omics (transcriptome) | `#3CAC9C` |
| competitor (SurvPath, PIBD) | `#BC6A79` |
| published field | `#8C8C8C` (the one neutral) |

**Ordered ladder** (risk tertile, capacity *d*, tie threshold — anything with a direction):
`#667DB8 #CB7F9F #A4A8B6 #CAAEDF #9CD0F3 #F3CFED #DEEFFF`, dark → pale. For three risk groups the
build **computes** sRGB luminance and picks the three entries with the largest minimum gap rather
than my guessing, and asserts the gap ≥ 0.05 so the panel survives greyscale.

Ink `#28303F`. Fills solid. Arial everywhere **including mathtext** (which has its own font stack
and would otherwise render every `$…$` in DejaVu). `pdf.fonttype 42` — required by this venue
family, and it would violate an AAAI submission, so the module says so where it is set.

Two assertions run at import: the ladder needs a minimum adjacent luminance gap of 0.05; the
nominal set needs a luminance spread ≤ 0.22 and no colour shared with the ladder.

---

## Figure 1 — the variable this benchmark left out

*Exists as a 2-panel draft; this replaces it.* One claim: the field's clinical baseline is built on
a covariate that cannot work in this disease, and the variable that can was in the released files.

| panel | content | encoding | source |
|---|---|---|---|
| **a** | Nine published entries, our arm, and the two clinical constructions | horizontal bars; ours `#6D98D3`, published `#8C8C8C`, stage-clinical `#89A467`, grade-clinical `#BC6A79`; dashed rule at the stage arm so the six published values below it are countable | `reviewer-gaps.json`, `Table 2` values |
| **b** | **Why grade cannot work**: level composition of grade vs stage in this cohort | two stacked bars, ordered ladder by level; annotate "94.4% one level, entropy 0.311" against "5 levels, entropy 0.665" | `why-grade-fails.json` |
| **c** | The one-column swap across all five TCGA studies | dumbbell: grade `#BC6A79` → stage `#89A467`, one row per study, gain annotated | `why-grade-fails.json` |
| **d** | Grade entropy vs the swap's benefit, five studies | scatter, the three studies with a usable grade filled and the two without ringed and labelled "no grade in file"; Pearson −0.976 in the corner **with "n = 3, consistency check not an estimate"** | `why-grade-fails.json` |

Panel **d** is the one that must not overstate. Three points cannot support a correlation, and the
annotation says so inside the panel rather than only in the caption.

---

## Figure 2 — the method and what it is worth

Six panels, the paper's spine.

| panel | content | encoding | source |
|---|---|---|---|
| **a** | **Forest plot**, **six** prespecified comparisons | point + 95% CI, cases resampled; zero line; raw *p* and **Holm *q*** printed at right; **the four that fail Holm drawn open, the two that survive filled** — and the two input-parity comparisons are among the four that fail, which the panel must not hide | `pibd-parity.json → holm_family_of_six` |
| **b** | **Kaplan–Meier by risk tertile — ours** | the three luminance-selected ladder colours (`#667DB8` high / `#CAAEDF` middle / `#DEEFFF` low, min gap **0.272**), **risk table beneath at 0–84 mo**, log-rank χ²=**62.4**, p=2.8e−14; high-risk median **18.3 mo** marked, and low and middle annotated **"median not reached"** rather than left blank | `figure-source-data.json → F1` |
| **c** | **Kaplan–Meier — the incumbent as released, and the same method given stage** | two mini-panels sharing a y-axis. Left: SurvPath **as the benchmark publishes it**, χ²=**7.7**, tertile medians 106.1 / 65.7 / **88.0** — the high-risk third's median exceeds the middle third's. Right: the **same predictions plus the same clinical block**, χ²=**31.5**, medians 106.1 / 57.3 / 24.0, monotone. **The panel is the paper's thesis in one image: the missing variable is what restores the ordering** | `figure-source-data.json → F1` |
| **d** | **Time-dependent AUC vs horizon**, 20-point grid over 3-29 mo | one line per arm in nominal colours, marker + linestyle so the panel survives greyscale; **n-at-risk printed along the top axis**, since the right-hand end rests on 114 patients and a curve that does not say so is misleading. Ours runs 0.768 to 0.732 and is above every other arm at every horizon | `figure-source-data.json` -> F2 |
| **e** | **Ablation**, all seven modality subsets | horizontal bars, coloured by which modalities are present; leave-one-out deltas and their paired *p* annotated | `ablation-generalisation.json → B_ablation_bladder` |
| **f** | **Cost against benefit** | log-x parameters against concordance. **Three real points now**: ours (1,049, 0.7225), SurvPath (24,702,532, 0.6897) and PIBD (26,930,702, 0.6982) — each competitor's count printed by **its own trainer** during the runs reported here. The other eight entrants stay on a labelled "not reported" band at their published concordance | `parameter-counts.json` + Table 2 |

**Two corrections this document earned by existing.** (i) Panel **c** as first written attributed the
non-monotone medians 106.1 / 65.7 / 88.0 to *SurvPath + clinical*. They belong to **SurvPath alone,
as released** — the parity comparator is monotone. Drawing before checking would have put a wrong
attribution into a figure. (ii) The panel must not be read as "the incumbent is broken": its *event
counts* order correctly (31 / 39 / 43 of 120 / 119 / 120). What fails is the **strength** of the
separation — a 26 / 33 / 36 % event fraction across tertiles, against **15 / 29 / 51 %** for ours. The
caption states it in those terms, because event fractions are robust where a median crossing is not.

Panel **f** was the one real design risk and **the PIBD run removed most of it.** Both competitors we
could rerun print their own parameter count at initialisation, so two of the three points are
measured rather than reconstructed, and the band now covers only the eight entrants we genuinely
cannot place. Four orders of magnitude on the x-axis against 113 events is the whole argument of the
method section in one panel.

---

## Figure 3 — what limits a method at 359 patients and 113 events

The regime evidence. This figure is why the method is parameter-free, and it is the one a
methods-minded referee will read first.

| panel | content | encoding | source |
|---|---|---|---|
| **a** | **In-sample fit vs the corrected null vs out-of-fold**, by capacity *d* ∈ {2,4,8,16} | three lines; shade the region between real and null — *that shaded gap is all the evidence there is*, and it shrinks from d=8 to d=16 | `baseline-ledger.md`, oracle-correction table |
| **b** | **Weight sweep**, 21 points | curve with the maximum marked at exactly w = 0.50; horizontal rules for the in-sample linear-fusion oracle (0.7268) and the gated oracle (0.7492) | `fusion.json` |
| **c** | **Gated fusion against its permuted-stratum control** | two bars, control drawn **above** the real value, annotated "control wins → component void" | `c2.json → C5` |
| **d** | **Power curve**: detectable margin vs power at this cohort's SE (0.0210) | curve over 101 margins; vertical band covering ALL EIGHT observed gaps between consecutive published entries (0.002–0.029, median 0.010) — **the median gap sits at power 0.076 and even the largest reaches only 0.282**; three margins marked: ours-vs-PIBD 0.0249 (power **0.220**), ours-vs-SurvPath 0.0340 (**0.367**), ours-vs-table-incumbent 0.0422 (**0.520**); 80% power at **0.0590** | `figure-source-data.json` -> F3 |
| **e** | **Two unreported sources of variation, both larger than the benchmark's own resolution** | left half: SurvPath's five seed values as a strip, SD 0.0094, published 0.625 marked. Right half: PIBD at its best-validation checkpoint (0.6609) against its final epoch (0.5900), the **+0.0709** bracket labelled "checkpoint chosen on the fold it is scored on". The published-gap band from (d) — 0.002–0.029, median 0.010 — runs across both halves as a shared reference, which is what makes the comparison land | `survpath-multiseed.json`, `pibd-parity.json → published` |

Panels **d** and **e** make the same point by three independent routes — a power calculation, a
measured seed spread, and a measured checkpoint-selection effect — and putting them adjacent is
deliberate. **The checkpoint effect (+0.0709) is the largest single number in this figure and it is
about the benchmark, not about any method**, so it must not be drawn in a nominal colour that reads
as "a competitor's result".

---

## Figure 4 — the biology

Six panels. Built so the central claim could have failed: if the slide arm were a molecular proxy,
panel **c** would show it.

| panel | content | encoding | source |
|---|---|---|---|
| **a** | **Are the eight signatures prognostic** | horizontal *z* bars, signed; BH q < 0.05 filled, above it open; dashed rules at z = ±1.96; direction annotated so a reader sees all four significant ones point as the literature says | `biology.json → B2` |
| **b** | **Signature × arm correlation heatmap** | 8 signatures × 3 arms, **diverging** map built from the ordered ladder, `✱` where BH q < 0.05; row order by the slide arm's \|ρ\| | `biology.json → B3_B4` |
| **c** | **The conditioning test** | grouped bars: slide and omics arms on all pairs, on omics-tied pairs, and on pairs neither omics nor stage separates; chance line at 0.5 | `biology.json → B5` |
| **d** | **Modalities swap by molecular subtype** | slope graph, basal → luminal, one line per arm in nominal colours; the clinical and slide lines cross, the ours line is flat | `biology.json → B6` |
| **e** | **Volcano**: slide arm against all 275 pathways | ρ on x, −log₁₀ q on y; BH threshold as a horizontal rule; the four strongest metabolic hits labelled; "59 of 275 at q<0.05" in the corner | `biology.json → B3b` |
| **f** | **Swimmer plot**, the ten disagreement cases | one bar per patient to its observed time, event as a filled marker and censoring as an open one; up-revised above, down-revised below; stage printed at the left so the II–III / IV split is visible at a glance | `biology.json → case_studies_rule2` |

Panel **f** must show the failures, since the rule admitted them: the up-revised patient censored at
2.0 months and the down-revised one with an event at 20.0. A swimmer plot shows censoring natively,
which is exactly why it is the right form here and a bar chart of medians is not.

---

## Supplementary figures — expanded, because J8 established they are standard here

The earlier plan had one three-panel supplement, written when it was unclear how much supplementary
material Briefings in Bioinformatics accepts. **J8 settled it**: all three sampled BiB papers carry
supplementary material, and SurvBoard (bbaf521) carries supplementary methods, Figures S1–S5 and
Tables S1–S3. So the supplement is a normal part of the submission, not a concession, and the three
panels that were being compressed into one figure get their own.

**One constraint governs all of them**: BiB does not copyedit or typeset supplementary material, so
what we hand over is exactly what a reader sees. These are built to the same standard as the main
figures — same palette, same Arial, same authored-at-final-size rule — not to a lower one.

| figure | panels | size (in) | content | source |
|---|---|---|---|---|
| **S1** | 2 | 6.30 × 3.60 | Seven falsified components: effect against the 0.0145 stability bar, each beside its matched control; the one **void on its own control** marked. Second panel: the same seven as a timeline, so the order they were tried in is visible | `contribution-design.md` |
| **S2** | 2 | 6.30 × 3.60 | Five-cohort generalisation: ours against the best verified-fold method per study, **four of five lost**. Second panel: the stage-vs-grade swap in the same five studies, which is the finding that *does* replicate | `ablation-generalisation.json`, `why-grade-fails.json` |
| **S3** | 2 | 6.30 × 3.60 | Six further slide encoders against TITAN; and the seven-encoder ensemble being significantly worse than TITAN alone | `multi.json` |

**S2 is the honest one and it must not be softened.** Putting the losses and the replication in the
same figure is what makes the paper's own distinction visible: the five-cohort evidence supports the
*stage* finding and not the bladder number, and a supplement that showed only the second half would
be making a claim the paper explicitly declines to make.

## Data still to compute before drawing — **all closed**

Everything the five figures need now exists as a file. One source-data file per figure set, which is
also what a journal Source Data upload requires, so it is built now rather than assembled at
submission from numbers that live only inside a PDF.

| was | now | where |
|---|---|---|
| KM step functions + risk tables | **done**, 7 arms × 3 tertiles, at-risk at 0–84 mo, medians incl. "not reached" | `figure-source-data.json → F1` |
| dense time-dependent AUC | **done**, 20 horizons over 3–29 mo, 6 arms; the 3 published horizons reproduce exactly (asserted in code) | `→ F2` |
| power curve | **done**, 101 margins at SE 0.0210; 80% power at **0.0590**, our two parity margins at **0.367** and **0.220** | `→ F3` |
| ladder triple, luminance-selected | **done**, `#667DB8` / `#CAAEDF` / `#DEEFFF`, min adjacent gap **0.272** — five times the 0.05 floor | `→ F4` |
| parameter counts | **done**, and better than expected: both rerunnable competitors print their own | `parameter-counts.json` |

The four assertions that make the file self-checking, all passing: the three published td-AUC
horizons must equal the confirmatory run's; the power curve must pass through the already-reported
0.0589 at 80%; the already-reported power 0.355 must land on it; the ladder triple's gap must clear
0.05. A source-data file that cannot contradict the manuscript is not evidence, so each of these was
written to be able to fail.

## What this set deliberately does not include

- **No representative slide images or attention heatmaps.** We hold precomputed embeddings, not the
  whole-slide images, so any such panel would be fabricated. A computational-pathology paper without
  a picture of tissue is unusual and this one has a reason.
- **No decorative raster schematic.** A method diagram, if one is wanted, is drawn as vector
  primitives.
- **No graphical abstract yet** — the venue lock says Briefings in Bioinformatics requires one, but
  the venue is not fixed, and it is derived from Figure 2 once it is.
