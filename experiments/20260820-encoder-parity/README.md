# 2026-08-20 — encoder parity on the comparators' own features

## Observed

```
$ bash experiments/20260820-encoder-parity/run.sh
known-answer  titan slide 0.6596 (want 0.6596)  full 0.7225 (want 0.7225)  REPRODUCED
chief         slide 0.5836  ModRank 0.7009  (359 of 359 cases covered)
wrote experiments/20260820-encoder-parity/results/encoder-parity.json
14.1 s wall, CPU only
```

The known-answer line is the reason the second line can be read at all. A new code path producing a
new number is not evidence; the same path reproducing the shipped sweep's TITAN values to four
decimals is what makes it one.

CHIEF's slide arm at 0.5836 also equals its own development ledger rows, ids 72 and 100, scored
during the campaign by different code. That is a second, independent agreement.

## Interpretation

The full method on CHIEF reaches 0.7009. Every published entry on this benchmark sits below it,
including the highest at 0.6910, while the comparator that consumed those very features reruns at
0.6118. The slide arm on CHIEF alone is 0.5836, below eight of the nine published entries.

So the lead is not an artefact of reading a better slide encoder than the methods it is compared
against. It is also not independent of the encoder: TITAN's 0.7225 against CHIEF's 0.7009 puts the
encoder's contribution at about 0.022, which is a fifth of the margin over the incumbent and not
nothing.

Single seed, seed 0, which is what every entry in the encoder sweep is. The headline 0.7212 is a
five-seed mean of the amended arm and is a different quantity.

## Not run

The other direction, the comparators on TITAN, and this is a matter of definition rather than of
budget. SurvPath and PIBD consume a bag of 4,096 patch tokens and attend over them; the public TITAN
release is one vector per slide, so feeding a bag of one destroys the premise of the architecture.
Doing it properly means obtaining the tile-level encoder, re-extracting features for every slide
from the whole-slide images, and retraining both models. That is a large download and GPU-days, and
neither was requested or approved.
