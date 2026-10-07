Translated from the Korean plan file committed at e579bc2d20ce5bae611f18aafe24ebb11f47d742; the original file is authoritative.

# POG570 confirmatory hypotheses

Fixed order. Each one-sided α = 0.05. One sample per patient. Each hypothesis resamples its analysis set by patient 2,000 times. Seed 20261001. Percentiles 2.5 and 97.5. The random-generator call order is H1, H2, H3. The difference is the HostMix-TOO statistic minus the baseline statistic.

| Order | Hypothesis | Test |
| --- | --- | --- |
| H1 | In the at-risk set, the host-attraction rate of HostMix-TOO is lower than that of the baseline | Exact one-sided McNemar test on the host-attraction indicator. A favorable discordant pair is a biopsy that is a host-attraction error for the baseline only. p is the probability that a Binomial(n_disc, 0.5) variable is at least the observed favorable count. If n_disc = 0, p = 1 |
| H2 | In the evaluation set, the organ top-1 accuracy of HostMix-TOO is higher than that of the baseline | Exact one-sided McNemar test. A favorable discordant pair is a biopsy correct for HostMix-TOO only. Tested only if H1 is met |
| H3 | In the native-truth set, the top-1 accuracy of HostMix-TOO is not lower than that of the baseline by more than 10 percentage points | The 5th percentile of the same H3 bootstrap differences is greater than −0.10. Judged only if H2 is met |

H1, H2 and H3 all report the difference and a two-sided 95% confidence interval. H3 also reports the one-sided 5th percentile.

A prediction of NA is incorrect.

## Analysis sets

- Evaluation set: the organ is neither excluded nor NA. n = 512.
- At-risk set: evaluation samples whose standard site has a non-empty native-organ set and whose true organ lies outside that set. n = 378.
- Native-truth set: the true organ lies inside that set. n = 91.

## Pre-specified descriptive outputs

No multiplicity adjustment. Point estimates are reported.

- Site groups: liver, lymph_node, lung, soft_tissue, remainder. Remainder is any standard site other than those four. For each method: n, top-1, n_at_risk, host_rate.
- Tumor-content tertiles T1, T2, T3, from qcut with 3 bins on the 512 evaluation values, duplicates dropped, closed on the right. The confirmatory run does not recompute the boundaries. For each method: n, top-1, n_at_risk, host_rate.

Analyses that are not in the preregistration are reported separately as post hoc.
