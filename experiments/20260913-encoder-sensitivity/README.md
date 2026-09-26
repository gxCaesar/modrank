# 20260913-encoder-sensitivity: is the five-study validation a property of the method or of TITAN?

The validation rests on a slide arm built from one pretrained encoder. A reader is entitled to ask
whether three of five studies separate because the recipe works or because that encoder is good, and
the paper already concedes that the encoder is worth about 0.02 on bladder. So the whole thing was
rerun with a different encoder and nothing else changed.

**GigaSSL-GigaPath** is the alternative, chosen because it is the only encoder on disk that covers
all five studies: 359 bladder, 863 breast, 298 colorectal, 394 head and neck and 318 stomach cases.
It is a different architecture, trained by a different group under a different objective, so
agreement between the two is not agreement between two versions of the same thing.

## Observed: the same three studies, under both encoders

ModRank minus the corrected clinical arm, 95% interval from 6,000 patient resamples.

| study | TITAN | GigaSSL-GigaPath | separates under both |
|---|---|---|---|
| bladder | +0.0774 [+0.0284, +0.1264] | **+0.0569 [+0.0034, +0.1099]** | **yes** |
| breast | +0.1208 [+0.0525, +0.1904] | **+0.0847 [+0.0164, +0.1527]** | **yes** |
| head and neck | +0.0726 [+0.0206, +0.1231] | **+0.0580 [+0.0078, +0.1101]** | **yes** |
| stomach | +0.0495 [−0.0127, +0.1122] | +0.0501 [−0.0093, +0.1138] | no, under both |
| colorectal | +0.0236 [−0.0686, +0.1197] | +0.0301 [−0.0665, +0.1283] | no, under both |

**Three of five under TITAN, the same three of five under GigaSSL, and the same two that do not.**
The grade-control comparison agrees too: colorectal alone separates under either encoder.

## What this does and does not establish

It establishes that the conclusion is not an artefact of the encoder: the two agree on which studies
separate from a correctly specified clinical reference.

They disagree on magnitude in a way worth stating precisely rather than in summary, because the
first draft of this file said "every ModRank value is lower with GigaSSL" and that is false.
Concordance falls in the three studies that separate — bladder 0.6992 to 0.6788, breast 0.6892 to
0.6533, head and neck 0.6179 to 0.6032 — and **rises** in the two that do not, stomach 0.6311 to
0.6366 and colorectal 0.7141 to 0.7210. A weaker slide representation lowering exactly the studies
where the slide was contributing, and not the ones where it was not, is a second reading of where
the slide arm is doing work.

It does not establish that either encoder is adequate, nor that a third would agree. It is one
alternative, chosen for coverage rather than for being the best available.

## Known answers

Skipped on purpose. The committed point estimates describe the TITAN arm, and this run replaces it;
checking a rebuilt arm against them would fail by construction. The protocol, folds, clinical block,
transcriptome arm, estimator, penalty grid and bootstrap are unchanged, which is what makes the
comparison a swap of one component rather than a different experiment.
