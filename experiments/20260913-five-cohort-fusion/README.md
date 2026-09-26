# 20260913-five-cohort-fusion: the comparison that only one of five studies had

The development cohort compares the equal-weight rank average against a concatenated Cox and a
stacked Cox, and beats both in 24 of 24 re-partitions. The benchmark's other four studies never had
that comparison at all: they scored the three arms and their average and stopped. So any statement
about whether fitting the fusion helps rested on one study out of five.

This ran both alternatives on all five, on the released folds, under the five-study protocol
(benchmark clinical file, seed 0). Concatenation gets the penalty grid extended by 32,768 in its
favour. Stacking learns three weights per fold from the arms' inner out-of-fold percentiles, copied
from s19 rather than re-derived, because weights learned on the scores they weight would be fitted
on the validation fold through the back door.

## Observed

| study | ModRank | concatenated | ModRank − concat | stacked | ModRank − stacked |
|---|---|---|---|---|---|
| bladder | 0.6992 | 0.6801 | +0.0191 [−0.0150, +0.0534] | 0.6932 | +0.0059 [−0.0143, +0.0267] |
| **breast** | 0.6892 | **0.7571** | **−0.0679 [−0.1384, +0.0029]** | 0.7180 | −0.0286 [−0.0754, +0.0202] |
| colorectal | 0.7141 | 0.6530 | +0.0611 [−0.0206, +0.1442] | 0.6919 | +0.0219 [−0.0329, +0.0776] |
| head and neck | 0.6179 | 0.5966 | +0.0217 [−0.0170, +0.0594] | 0.6171 | +0.0009 [−0.0206, +0.0228] |
| stomach | 0.6311 | 0.6073 | +0.0242 [−0.0281, +0.0779] | 0.5862 | +0.0458 [−0.0152, +0.1081] |

**Ten comparisons, and not one interval excludes zero in either direction.**

## What it settles, and what it takes away

It takes away a sentence. A draft title said "fitting the fusion adds nothing", and this run is why
that wording did not ship. The strict version survives — no fitted fusion beats the rank average
with an interval excluding zero in any study — but in breast a concatenated Cox is ahead by 0.0679
at p=0.061, and a title that says a fitted fusion adds nothing while the archive shows that number
is a title contradicted by its own data.

Note also that bladder's advantage is protocol-dependent. Under the main protocol, five seeds and
the incumbent's clinical block, the rank average beats concatenation by +0.0464 and stacking by
+0.0111 in 24 of 24 re-partitions. Under this protocol, one seed and the benchmark's clinical file,
the same comparisons are +0.0191 and +0.0059 with intervals covering zero. Both are reported.

What it settles is the paper's own thesis, now with ten comparisons behind it rather than an
assertion: **at these sample sizes this benchmark separates a correct clinical reference from an
incorrect one and does not separate one fusion rule from another.** The reference correction is
+0.0972 [+0.047, +0.149]. Every method-versus-method difference here covers zero. That contrast is
the finding, and it is stronger for the null half being measured rather than assumed.
