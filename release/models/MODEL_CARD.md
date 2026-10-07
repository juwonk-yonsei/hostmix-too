# Model card: BASE-Z, SA-Z and a post hoc model

File names: BASE_Z = baseline; SA_Z = HostMix-TOO; MIX_Z0_posthoc = the unstandardized mixture model (post hoc).

Name: HostMix-TOO is the locked SA-Z procedure.

## Intended use

Research use only. These models are not for clinical decisions.

## Training data

Primary-tumor expression from the UCSC Toil recomputation of TCGA, 7486 rows, 32 labels. GTEx-ref profiles from ten host tissues, 1103 samples, are the host pool. TARGET samples were counted and not used for training. This package does not contain the expression matrices.

## Features and classes

Z is a within-sample rank-normal score. Ranks use every gene in G, including zeros, and the model keeps the 5000 B0 columns in the stored order. The class list is in each portable JSON file. Colorectal is the COADREAD label. Lung combines LUAD and LUSC. Project labels are summed to those organ classes.

Both models are multinomial L2 logistic regression, class_weight balanced, solver lbfgs, max_iter 5000.

| Model | Scaler | C | How C was chosen | Rows |
| --- | --- | --- | --- | --- |
| BASE-Z | none | 0.1 | Fixed | 7486 |
| SA-Z | StandardScaler | 0.03 | 3-fold grouped cross-validation, macro-F1, grid 0.03, 0.1, 0.3, smallest C within 1e-12 of the best | 37430 |

## Augmentation

SA-Z uses one pure profile and four mixtures per tumor. The mixture weight is Uniform(0.15, 1.0), half-open. The host is a uniform draw of a host-pool tissue, then a uniform GTEx-ref sample of that tissue. The SA seed is 890338470.

## Performance summary

POG570. H1 host rate, n = 378: BASE-Z 0.228, SA-Z 0.024, difference -0.204 (95% CI -0.243 to -0.164). H2 top-1, n = 512: 0.682 and 0.756, difference 0.074 (95% CI 0.039 to 0.107). H3: non-inferiority (margin -0.10) was not shown: difference -0.066 (95% CI -0.132 to 0.000); one-sided 95% lower bound -0.121.

Auxiliary RNA-seq. AH1 host rate, n = 427: 0.204 and 0.056, difference -0.148. AH2 top-1, n = 729: 0.527 and 0.698, difference 0.171. AH3 non-inferiority was shown: difference 0.257, one-sided 95% lower bound 0.216. AS1-AS3 were not rejected after Holm correction.

Pre-specified description, not a hypothesis test. On MET500, POG570, and the auxiliary set, SCOPE and CUP-AI-Dx host rates were 0.122 to 0.304. SA-Z host rates were 0.024 to 0.056. SA-Z did not have the highest top-1 on every cohort.

Post hoc ablation: the mixture share of the host-rate gap was at least 0.810 in all three cohorts. That sentence is a description, not a preregistered result.

## Known limits

Hosts outside the training pool were a failed secondary test. Missing microarray genes move accuracy. Assay shift, including IMvigor210 esophagus calls, is a different failure from host attraction. POG570 sample-level predictions are not in this package.

## Post hoc model: the unstandardized mixture model

Files: `MIX_Z0_posthoc.joblib` (scikit-learn) and `MIX_Z0_posthoc.npz` with `MIX_Z0_posthoc.json` (portable, read without pickle by `scripts/predict_portable.py`). This model is not HostMix-TOO.

Status: post hoc ablation model; not evaluated in a preregistered confirmation.

Training data and settings: the 7486 TCGA primary tumors with the HostMix-TOO mixtures (the same mixture records as SA-Z, one pure profile and four mixtures per tumor, 37430 rows), the same Z features without feature standardization, and C = 0.1, fixed. Multinomial L2 logistic regression, class_weight balanced, solver lbfgs, max_iter 5000. Predictions rebuilt from these files match the stored ablation predictions on MET500 (437 biopsies, 0 top-1 mismatches).

Performance, post hoc (Additional file 1: Tables S8 and S9). Top-1 uses the evaluation set, the host-attraction rate the at-risk set and native-truth top-1 the native-truth set.

| Cohort | Model | Top-1 | Host-attraction rate | Native-truth top-1 |
| --- | --- | --- | --- | --- |
| MET500 (n = 437, 361, 76) | Baseline | 0.689 | 0.136 | 0.803 |
| | Unstandardized mixture model | 0.728 | 0.044 | 0.789 |
| | HostMix-TOO | 0.728 | 0.044 | 0.789 |
| POG570 (n = 512, 378, 91) | Baseline | 0.682 | 0.228 | 0.868 |
| | Unstandardized mixture model | 0.777 | 0.034 | 0.846 |
| | HostMix-TOO | 0.756 | 0.024 | 0.802 |
| Auxiliary RNA-seq (n = 729, 427, 315) | Baseline | 0.527 | 0.204 | 0.473 |
| | Unstandardized mixture model | 0.706 | 0.082 | 0.727 |
| | HostMix-TOO | 0.698 | 0.056 | 0.730 |
| Microarray (n = 536, 207, 341) | Baseline | 0.724 | 0.256 | 0.871 |
| | Unstandardized mixture model | 0.776 | 0.092 | 0.891 |
| | HostMix-TOO | 0.455 | 0.101 | 0.554 |

In the microarray evaluation set the unstandardized mixture model made 11 esophagus calls (baseline 30, HostMix-TOO 147). When as many randomly chosen classifier genes were removed from POG570 and MET500 profiles as each of four microarray platforms lacks, its mean top-1 loss was 0.006 (baseline 0.006, HostMix-TOO 0.031).

Recommended use: research use; confirm before application.

## Citation

Cite the archive (https://doi.org/10.5281/zenodo.23205336) and the GTEx sentence in `host_pool/CITATION.txt`. The article is under review; its DOI will be added on publication.
