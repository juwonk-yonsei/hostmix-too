Translated from the Korean plan file committed at ba3c074e96e0cfd0124e54a24e08964849e7879e; the original file is authoritative.

# Auxiliary confirmatory hypotheses

## Primary family

RNA-seq layer, one sample per patient, fixed order, each one-sided α = 0.05.

| Order | Hypothesis | Test |
| --- | --- | --- |
| AH1 | In the at-risk set, the host-attraction rate of HostMix-TOO is lower than that of the baseline | Exact one-sided McNemar test |
| AH2 | In the evaluation set, the top-1 accuracy of HostMix-TOO is higher than that of the baseline | Exact one-sided McNemar test. Tested only if AH1 is met |
| AH3 | In the native-truth set, the top-1 accuracy of HostMix-TOO is not lower than that of the baseline by more than 10 percentage points | One-sided 95% lower bound of the patient-bootstrap difference (2,000 replicates) > −0.10. Judged only if AH2 is met |

## Secondary family

RNA-seq layer, one sample per patient, Holm correction, family-wise one-sided α = 0.05. The secondary family is tested regardless of the primary family.

| Hypothesis | Statement | Test |
| --- | --- | --- |
| AS1 | In the native-truth set, the top-1 accuracy of the gated model is higher than that of HostMix-TOO | Exact one-sided McNemar test |
| AS2 | In the pool-out at-risk set, the host-attraction rate of the 22-tissue model is lower than that of HostMix-TOO | Exact one-sided McNemar test |
| AS3 | In the evaluation set, the top-1 accuracy of the gated model is higher than that of HostMix-TOO | Exact one-sided McNemar test |

Error-rate control is separate for each family. The primary family uses a fixed order and the secondary family uses Holm, each at family-wise α = 0.05. The two families are not controlled together. The paper's primary claims rest only on the primary family.

Effect size for every hypothesis: the difference and a 95% patient-bootstrap interval (2,000 replicates, seed 20261001, percentile method), plus the discordant counts.

Exact one-sided McNemar p = P(Binomial(m, 0.5) ≥ favorable count). m = favorable + unfavorable. If m = 0, the test is not rejected.

## Descriptive analyses named before unlock

No multiplicity adjustment.

The microarray layer repeats the same comparisons as effect sizes and intervals only.

A cohort meta-analysis uses cohorts with at least 10 at-risk samples and one sample per patient. d_i = (c_i − b_i) / n_i, where b_i is the number of patients who are a host-attraction error for the baseline only, c_i is the number who are a host-attraction error for HostMix-TOO only, and n_i is the number of at-risk patients. The variance is v_i = [(b_i + c_i) − (b_i − c_i)^2 / n_i] / n_i^2. If b_i or c_i is 0, 0.5 is added to both counts for the variance calculation only. The summary is a DerSimonian–Laird random-effects estimate with a 95% interval and I^2, separately by layer.

Analyses that are not in the preregistration are reported separately as post hoc.
