# Additional file 1. Supplementary methods, tables and figures

## 1. Data processing

The matrix is distributed as log2(TPM + 0.001), and we converted it to linear transcripts per million (TPM) as max(2^x − 0.001, 0).

Gene identifiers were mapped to symbols with GENCODE version 23. The gene universe is the 17552 HGNC protein-coding symbols present in the Toil matrix, the MET500 matrix and the POG570 matrix. Previous HGNC symbols were updated to the current symbol. Aliases were not used.

Within each cohort, values were converted to a linear scale summing to 10^6 per sample. FPKM values were rescaled. Counts were divided by the GENCODE version 23 exon-union length and then rescaled. Log-scale values were exponentiated. PRINCE provides log2-transformed normalized values, which were converted by 2^v − 1 without a length correction.

Microarray probes matching more than one gene were removed. For each remaining gene the probe with the highest mean signal was kept. GSE71729 and GSE74685 are two-channel arrays. The sample channel was used and the reference channel was not.

## 2. Cohort search and selection

The search log has 537 rows and 502 study identifiers (`results/stage5/aux/screening_log.tsv`). The last row for each identifier is the decision used here. Of those 502 studies, 348 were excluded, 142 were held after a title-only screen, and 12 were taken forward. These three counts sum to 502. The 35 extra rows are later updates of an identifier already in the log.

Excluded studies, by the criterion in `config/aux_eligibility_A5.md`: 144 failed criterion 2 (not a human tumor-tissue sample; cell lines, xenografts, organoids, blood or single-cell profiles), 67 failed criterion 3 (no per-sample biopsy site), 4 failed criterion 4 (the primary organ did not map to an organ group), 8 failed criterion 6 (fewer than 10 at-risk samples in the metadata), and 125 failed criterion 7 (the study belongs to TCGA, GTEx, MET500, POG570 or CCLE). These five counts sum to 348.

The 142 held studies are GEO series whose title and summary were read and whose sample characteristics were not. Held means that the title and summary did not name a cell line, organoid, xenograft, single-cell or spatial profile, or a blood, plasma or serum sample, so the series was not excluded under criterion 2, and that it was not reviewed further against the other criteria; no reason was recorded for not reviewing these series further. The 12 studies taken forward are the cohorts in Table 1: five that passed the biopsy-site screen (IMvigor210, Anders, DFCI melanoma, PRINCE and SU2C/PCF) and seven marked included (GSE50760, GSE209998, GSE41258, GSE14018, GSE71729, GSE74685 and FHCRC).

## 3. Duplicate screening and other exclusions

A sample was marked as a duplicate when the maximum Spearman correlation of log2(TPM + 1) with an external MET500 or POG570 sample was at least 0.98, or at least 0.95 and greater than its maximum correlation with another patient in the same cohort (`config/prereg_A6.md`). That rule excluded 76 SU2C/PCF samples. No sample in any other cohort was excluded as a duplicate.

Other exclusions used the organ-label rules, not the correlation rule. Normal tissue, polyps and samples with an empty tissue field were excluded under rule 8 (GSE209998, 6; GSE41258, 125; GSE50760, 18; GSE71729, 134). Cell-line samples were excluded under rule 2 (GSE41258, 12; GSE71729, 17). Neuroendocrine or small-cell features were excluded under rule 1 (SU2C/PCF, 17; GSE74685, 20 chromogranin-A-positive samples). A further 93 samples kept the cohort organ label but had no standard biopsy site, so they were outside the evaluation set (IMvigor210, 71; PRINCE, 21; Anders, 1).

## 4. Labels and sites

Table S1 gives the 32 tumor labels and the 26 organ groups. Table S2 gives the standard site, the host tissue or proxy, and the native organ set.

Fifty-nine MET500 specimens were excluded from the evaluation set because the annotated diagnosis did not map to an organ group. Fifty-eight POG570 samples were excluded for the same reason, leaving 512. In IMvigor210, 67 samples were metastases biopsied from the kidney. Their true organ, bladder, lies outside the kidney native set, so they are at risk, and they form 67 of the 71 biopsies in the pool-out at-risk set.

### Table S1. Tumor labels and organ groups

Thirty-two tumor labels are summed to 26 organ groups after prediction. Lung combines lung adenocarcinoma and lung squamous cell carcinoma. Colorectal combines colon and rectal adenocarcinoma.

| Tumor label | Organ group |
|---|---|
| BRCA | Breast |
| PRAD | Prostate |
| COADREAD | Colorectal |
| LUAD | Lung |
| LUSC | Lung |
| LIHC | Liver |
| CHOL | Biliary |
| PAAD | Pancreas |
| STAD | Stomach |
| ESCA | Esophagus |
| HNSC | HeadNeck |
| THCA | Thyroid |
| KIRC | Kidney |
| KIRP | Kidney |
| KICH | Kidney |
| BLCA | Bladder |
| OV | Ovary |
| UCEC | Uterus |
| UCS | Uterus |
| CESC | Cervix |
| SKCM | Melanoma |
| UVM | UvealMelanoma |
| GBM | Brain |
| LGG | Brain |
| ACC | Adrenal |
| PCPG | Adrenal |
| TGCT | Testis |
| THYM | Thymus |
| SARC | Sarcoma |
| MESO | Mesothelioma |
| DLBC | Lymphoid |
| LAML | Myeloid |

### Table S2. Standard biopsy sites and native organ sets

A biopsy is at risk when its true organ is outside the native set of its standard site. A native-truth biopsy is one whose true organ is inside that set. Lymph node uses spleen as a proxy, and bone marrow and bone use whole blood. Tissues listed here that are not in the ten-tissue host pool are outside the pool.

| Standard site | Host tissue or proxy | Native organ set | In host pool |
|---|---|---|---|
| liver | Liver | Liver; Biliary | yes |
| lung | Lung | Lung; Mesothelioma | yes |
| lymph_node | Spleen (proxy) | Lymphoid | yes |
| bone_marrow | Whole Blood (proxy) | Myeloid; Lymphoid | yes |
| bone | Whole Blood (proxy) | Myeloid; Lymphoid | yes |
| brain | Brain - Cortex | Brain | yes |
| adrenal | Adrenal Gland | Adrenal | yes |
| skin | Skin - Not Sun Exposed (Suprapubic) | Melanoma | yes |
| subcutaneous | Skin - Not Sun Exposed (Suprapubic) | Melanoma | yes |
| soft_tissue | Adipose - Subcutaneous; Muscle - Skeletal | Sarcoma | yes |
| kidney | Kidney - Cortex | Kidney | no |
| pancreas | Pancreas | Pancreas | no |
| peritoneum | Adipose - Visceral (Omentum) | Mesothelioma | yes |
| omentum | Adipose - Visceral (Omentum) | Mesothelioma | yes |
| pleura |  | Mesothelioma | no |
| bladder | Bladder | Bladder | no |
| thyroid | Thyroid | Thyroid | no |
| ovary | Ovary | Ovary | no |
| breast | Breast - Mammary Tissue | Breast | no |
| stomach | Stomach | Stomach | no |
| colon_rectum | Colon - Transverse | Colorectal | no |
| prostate | Prostate | Prostate | no |
| esophagus | Esophagus - Mucosa | Esophagus | no |
| head_neck | Minor Salivary Gland | HeadNeck | no |
| cervix | Cervix - Ectocervix | Cervix | no |
| other |  |  | no |

## 5. Host pool and proxies

The host pool is ten GTEx tissues. Lymph node and bone marrow have no GTEx tissue of the same name. The proxy was the candidate tissue with the highest Spearman correlation of mean expression against the Human Protein Atlas tissue query. Table S3 lists the pool counts and the eight correlations.

### Table S3. Host pool and proxy correlations

GTEx-ref sample counts for the ten host-pool tissues, then the eight Spearman correlations used to choose proxies for lymph node and bone marrow. The selected proxies are spleen for lymph node and whole blood for bone marrow.

| Host-pool tissue | GTEx-ref samples |
|---|---|
| Adipose - Subcutaneous | 156 |
| Adipose - Visceral (Omentum) | 98 |
| Adrenal Gland | 64 |
| Brain - Cortex | 51 |
| Liver | 55 |
| Lung | 145 |
| Muscle - Skeletal | 204 |
| Skin - Not Sun Exposed (Suprapubic) | 116 |
| Spleen | 49 |
| Whole Blood | 165 |

| Query | Candidate tissue | Genes compared | GTEx-ref samples | Spearman |
|---|---|---|---|---|
| lymph node | Spleen | 17539 | 49 | 0.8813875885578315 |
| lymph node | Whole Blood | 17539 | 165 | 0.8424373159977525 |
| lymph node | Cells - Ebv-Transformed Lymphocytes | 17539 | 53 | 0.8705113793189775 |
| lymph node | Small Intestine - Terminal Ileum | 17539 | 46 | 0.8305221329877467 |
| bone marrow | Spleen | 17539 | 49 | 0.8199266753031522 |
| bone marrow | Whole Blood | 17539 | 165 | 0.8832743760976061 |
| bone marrow | Cells - Ebv-Transformed Lymphocytes | 17539 | 53 | 0.8325539228733679 |
| bone marrow | Small Intestine - Terminal Ileum | 17539 | 46 | 0.7461527036575726 |

## 6. Mixture generation and seeds

Each HostMix-TOO training tumor contributes one unchanged profile and four mixtures. The tumor fraction of RNA is drawn from the uniform distribution on [0.15, 1.0). The host is a uniform draw of one pool tissue and then a uniform GTEx-ref sample of that tissue. The mixture seed recorded for HostMix-TOO is 890338470. Bootstrap intervals use 2000 patient resamples and seed 20261001. The random gene-removal control uses seed 20261023. Cross-validation uses seed 20261001.

## 7. Development comparison and exploratory variants

Table S4 is the MET500 comparison used during development. Table S5 is the set of variants fit between the two confirmations. Neither table is a hypothesis test.

### Table S4. MET500 development comparison

MET500 was the development cohort. Intervals are the bootstrap intervals from the development table. Candidate rank is the development ordering among the ten candidates; a blank rank means the classifier was not in that candidate list. This comparison is exploratory.

| Classifier | Top-1 | Top-1 interval | Host-attraction rate | Host-rate interval | Native-truth top-1 | Candidate rank |
|---|---|---|---|---|---|---|
| Baseline | 0.689 | 0.645 to 0.730 | 0.136 | 0.101 to 0.172 | 0.803 |  |
| Native-organ masking | 0.600 | 0.554 to 0.645 | 0.000 | 0.000 to 0.000 | 0.000 |  |
| Normal classes | 0.675 | 0.634 to 0.719 | 0.119 | 0.086 to 0.153 | 0.776 |  |
| Linear deconvolution | 0.709 | 0.668 to 0.751 | 0.089 | 0.061 to 0.118 | 0.816 |  |
| HostMix-TOO | 0.728 | 0.686 to 0.767 | 0.044 | 0.023 to 0.067 | 0.789 | 1 |
| Gene removal 20% | 0.677 | 0.636 to 0.721 | 0.119 | 0.085 to 0.153 | 0.789 | 5 |
| Site-specific mixtures | 0.714 | 0.670 to 0.753 | 0.030 | 0.014 to 0.049 | 0.737 | 2 |
| Site-specific + removal 20% | 0.709 | 0.668 to 0.751 | 0.030 | 0.014 to 0.049 | 0.711 | 3 |
| Site-specific + removal 40% | 0.705 | 0.661 to 0.746 | 0.028 | 0.011 to 0.045 | 0.724 | 4 |
| Gene sets | 0.568 | 0.522 to 0.613 | 0.141 | 0.106 to 0.178 | 0.697 |  |
| Gene-set form of Native-organ masking | 0.481 | 0.435 to 0.529 | 0.000 | 0.000 to 0.000 | 0.000 |  |
| Gene sets + host correction | 0.579 | 0.531 to 0.622 | 0.053 | 0.030 to 0.077 | 0.566 |  |
| Gene-set form of normal classes | 0.590 | 0.547 to 0.634 | 0.119 | 0.086 to 0.155 | 0.711 |  |
| Gene-set form of Linear deconvolution | 0.549 | 0.503 to 0.595 | 0.066 | 0.040 to 0.093 | 0.618 |  |
| Gene sets + mixtures | 0.670 | 0.627 to 0.712 | 0.044 | 0.023 to 0.066 | 0.737 | 6 |
| Gene-set form of Gene removal 20% | 0.574 | 0.531 to 0.618 | 0.141 | 0.106 to 0.178 | 0.697 | 10 |
| Gene-set form of Site-specific mixtures | 0.609 | 0.565 to 0.652 | 0.022 | 0.008 to 0.039 | 0.684 | 9 |
| Gene-set form of Site-specific + removal 20% | 0.616 | 0.570 to 0.659 | 0.025 | 0.011 to 0.042 | 0.684 | 8 |
| Gene-set form of Site-specific + removal 40% | 0.629 | 0.584 to 0.673 | 0.025 | 0.011 to 0.042 | 0.724 | 7 |

### Table S5. Exploratory variants

These variants were fit between the two confirmations. Simulation rates are unweighted means across tissues at tumor RNA fraction 0.6. The leave-one-host-out model that omitted liver had a POG570 host-attraction rate of 0.183. Across the 10 left-out tissues, the POG570 host-attraction rate of the matching leave-one-host-out model ranged from 0.021 to 0.183. These rows are exploratory.

| Variant | TCGA-test top-1 | In-pool host rate at rho 0.6 | Pool-out host rate at rho 0.6 | MET500 top-1 | MET500 host rate | POG570 top-1 | POG570 host rate | POG570 native-truth top-1 |
|---|---|---|---|---|---|---|---|---|
| Baseline | 0.982 | 0.111 | 0.120 | 0.689 | 0.136 | 0.682 | 0.228 | 0.868 |
| HostMix-TOO | 0.979 | 0.002 | 0.084 | 0.728 | 0.044 | 0.756 | 0.024 | 0.802 |
| Standardized, pure tumors | 0.981 | 0.159 | 0.106 | 0.666 | 0.139 | 0.670 | 0.222 | 0.868 |
| One mixture per tumor | 0.980 | 0.002 | 0.085 | 0.728 | 0.044 | 0.752 | 0.029 | 0.835 |
| Eight mixtures per tumor | 0.982 | 0.002 | 0.086 | 0.735 | 0.036 | 0.748 | 0.019 | 0.802 |
| Mixture lower bound 0.3 | 0.979 | 0.002 | 0.076 | 0.730 | 0.036 | 0.766 | 0.029 | 0.813 |
| Mixture lower bound 0.5 | 0.980 | 0.003 | 0.077 | 0.721 | 0.033 | 0.758 | 0.034 | 0.824 |
| Mixtures, C = 0.01 | 0.980 | 0.002 | 0.076 | 0.737 | 0.039 | 0.781 | 0.026 | 0.835 |
| Mixtures, C = 0.1 | 0.980 | 0.002 | 0.089 | 0.741 | 0.033 | 0.770 | 0.026 | 0.813 |
| 3-tissue model | 0.981 | 0.084 | 0.101 | 0.737 | 0.030 | 0.760 | 0.037 | 0.846 |
| 22-tissue model | 0.979 | 0.003 | 0.002 | 0.705 | 0.047 | 0.746 | 0.048 | 0.835 |
| Perceptron | 0.978 | 0.179 | 0.110 | 0.654 | 0.147 | 0.652 | 0.230 | 0.890 |
| Perceptron + mixtures | 0.980 | 0.003 | 0.114 | 0.696 | 0.030 | 0.750 | 0.024 | 0.813 |

## 8. POG570 alternatives and strata

### Table S6. POG570 alternatives and pre-specified strata

These outputs were named before the POG570 labels were unlocked. They are pre-specified descriptive results, not hypothesis tests.

Overall.

| Classifier | n | Top-1 | n at risk | Host-attraction rate | n native truth | Native-truth top-1 |
|---|---|---|---|---|---|---|
| Baseline | 512 | 0.682 | 378 | 0.228 | 91 | 0.868 |
| Gene sets | 512 | 0.596 | 378 | 0.246 | 91 | 0.791 |
| HostMix-TOO | 512 | 0.756 | 378 | 0.024 | 91 | 0.802 |
| Site-specific mixtures | 512 | 0.758 | 378 | 0.024 | 91 | 0.802 |
| Linear deconvolution | 512 | 0.723 | 378 | 0.143 | 91 | 0.879 |
| Normal classes | 512 | 0.709 | 378 | 0.180 | 91 | 0.857 |
| Native-organ masking | 512 | 0.637 | 378 | 0.000 | 91 | 0.000 |
| Gene sets + mixtures | 512 | 0.727 | 378 | 0.037 | 91 | 0.802 |
| Gene sets + host correction | 512 | 0.678 | 378 | 0.074 | 91 | 0.714 |

By biopsy-site group.

| Classifier | Site group | n | Top-1 | n at risk | Host-attraction rate | n native truth | Native-truth top-1 |
|---|---|---|---|---|---|---|---|
| Baseline | liver | 186 | 0.478 | 175 | 0.451 | 11 | 0.545 |
| Baseline | lung | 52 | 0.827 | 30 | 0.133 | 22 | 0.909 |
| Baseline | lymph_node | 92 | 0.804 | 89 | 0.000 | 3 | 0.667 |
| Baseline | remainder | 133 | 0.759 | 45 | 0.044 | 45 | 0.911 |
| Baseline | soft_tissue | 49 | 0.857 | 39 | 0.026 | 10 | 1.000 |
| Gene sets | liver | 186 | 0.339 | 175 | 0.514 | 11 | 0.364 |
| Gene sets | lung | 52 | 0.788 | 30 | 0.033 | 22 | 0.864 |
| Gene sets | lymph_node | 92 | 0.772 | 89 | 0.000 | 3 | 0.667 |
| Gene sets | remainder | 133 | 0.699 | 45 | 0.022 | 45 | 0.844 |
| Gene sets | soft_tissue | 49 | 0.755 | 39 | 0.026 | 10 | 0.900 |
| HostMix-TOO | liver | 186 | 0.726 | 175 | 0.046 | 11 | 0.545 |
| HostMix-TOO | lung | 52 | 0.846 | 30 | 0.000 | 22 | 0.909 |
| HostMix-TOO | lymph_node | 92 | 0.783 | 89 | 0.000 | 3 | 0.667 |
| HostMix-TOO | remainder | 133 | 0.737 | 45 | 0.000 | 45 | 0.822 |
| HostMix-TOO | soft_tissue | 49 | 0.776 | 39 | 0.026 | 10 | 0.800 |
| Site-specific mixtures | liver | 186 | 0.726 | 175 | 0.046 | 11 | 0.545 |
| Site-specific mixtures | lung | 52 | 0.846 | 30 | 0.000 | 22 | 0.909 |
| Site-specific mixtures | lymph_node | 92 | 0.783 | 89 | 0.000 | 3 | 0.667 |
| Site-specific mixtures | remainder | 133 | 0.744 | 45 | 0.000 | 45 | 0.822 |
| Site-specific mixtures | soft_tissue | 49 | 0.776 | 39 | 0.026 | 10 | 0.800 |
| Linear deconvolution | liver | 186 | 0.597 | 175 | 0.291 | 11 | 0.636 |
| Linear deconvolution | lung | 52 | 0.827 | 30 | 0.033 | 22 | 0.909 |
| Linear deconvolution | lymph_node | 92 | 0.793 | 89 | 0.000 | 3 | 0.667 |
| Linear deconvolution | remainder | 133 | 0.767 | 45 | 0.022 | 45 | 0.911 |
| Linear deconvolution | soft_tissue | 49 | 0.837 | 39 | 0.026 | 10 | 1.000 |
| Normal classes | liver | 186 | 0.543 | 175 | 0.371 | 11 | 0.545 |
| Normal classes | lung | 52 | 0.846 | 30 | 0.000 | 22 | 0.909 |
| Normal classes | lymph_node | 92 | 0.804 | 89 | 0.000 | 3 | 0.667 |
| Normal classes | remainder | 133 | 0.789 | 45 | 0.022 | 45 | 0.911 |
| Normal classes | soft_tissue | 49 | 0.796 | 39 | 0.051 | 10 | 0.900 |
| Native-organ masking | liver | 186 | 0.731 | 175 | 0.000 | 11 | 0.000 |
| Native-organ masking | lung | 52 | 0.481 | 30 | 0.000 | 22 | 0.000 |
| Native-organ masking | lymph_node | 92 | 0.783 | 89 | 0.000 | 3 | 0.000 |
| Native-organ masking | remainder | 133 | 0.459 | 45 | 0.000 | 45 | 0.000 |
| Native-organ masking | soft_tissue | 49 | 0.653 | 39 | 0.000 | 10 | 0.000 |
| Gene sets + mixtures | liver | 186 | 0.672 | 175 | 0.069 | 11 | 0.364 |
| Gene sets + mixtures | lung | 52 | 0.846 | 30 | 0.000 | 22 | 0.909 |
| Gene sets + mixtures | lymph_node | 92 | 0.772 | 89 | 0.000 | 3 | 0.667 |
| Gene sets + mixtures | remainder | 133 | 0.707 | 45 | 0.022 | 45 | 0.844 |
| Gene sets + mixtures | soft_tissue | 49 | 0.776 | 39 | 0.026 | 10 | 0.900 |
| Gene sets + host correction | liver | 186 | 0.613 | 175 | 0.149 | 11 | 0.545 |
| Gene sets + host correction | lung | 52 | 0.788 | 30 | 0.000 | 22 | 0.818 |
| Gene sets + host correction | lymph_node | 92 | 0.772 | 89 | 0.000 | 3 | 0.667 |
| Gene sets + host correction | remainder | 133 | 0.669 | 45 | 0.022 | 45 | 0.756 |
| Gene sets + host correction | soft_tissue | 49 | 0.653 | 39 | 0.026 | 10 | 0.500 |

By tumor-content tertile.

| Classifier | Tumor-content tertile | n | Top-1 | n at risk | Host-attraction rate |
|---|---|---|---|---|---|
| Baseline | T1 | 173 | 0.636 | 128 | 0.273 |
| Baseline | T2 | 168 | 0.685 | 124 | 0.274 |
| Baseline | T3 | 171 | 0.725 | 126 | 0.135 |
| Gene sets | T1 | 173 | 0.561 | 128 | 0.266 |
| Gene sets | T2 | 168 | 0.565 | 124 | 0.306 |
| Gene sets | T3 | 171 | 0.661 | 126 | 0.167 |
| HostMix-TOO | T1 | 173 | 0.682 | 128 | 0.055 |
| HostMix-TOO | T2 | 168 | 0.786 | 124 | 0.008 |
| HostMix-TOO | T3 | 171 | 0.801 | 126 | 0.008 |
| Site-specific mixtures | T1 | 173 | 0.688 | 128 | 0.039 |
| Site-specific mixtures | T2 | 168 | 0.792 | 124 | 0.016 |
| Site-specific mixtures | T3 | 171 | 0.795 | 126 | 0.016 |
| Linear deconvolution | T1 | 173 | 0.676 | 128 | 0.180 |
| Linear deconvolution | T2 | 168 | 0.726 | 124 | 0.194 |
| Linear deconvolution | T3 | 171 | 0.766 | 126 | 0.056 |
| Normal classes | T1 | 173 | 0.642 | 128 | 0.234 |
| Normal classes | T2 | 168 | 0.744 | 124 | 0.202 |
| Normal classes | T3 | 171 | 0.743 | 126 | 0.103 |
| Native-organ masking | T1 | 173 | 0.601 | 128 | 0.000 |
| Native-organ masking | T2 | 168 | 0.667 | 124 | 0.000 |
| Native-organ masking | T3 | 171 | 0.643 | 126 | 0.000 |
| Gene sets + mixtures | T1 | 173 | 0.665 | 128 | 0.055 |
| Gene sets + mixtures | T2 | 168 | 0.732 | 124 | 0.048 |
| Gene sets + mixtures | T3 | 171 | 0.784 | 126 | 0.008 |
| Gene sets + host correction | T1 | 173 | 0.630 | 128 | 0.109 |
| Gene sets + host correction | T2 | 168 | 0.673 | 124 | 0.097 |
| Gene sets + host correction | T3 | 171 | 0.731 | 126 | 0.016 |

## 9. Registered hypotheses

The following English text is a translation. The Korean plan files are authoritative.

Translated from the Korean plan file committed at e579bc2d20ce5bae611f18aafe24ebb11f47d742; the original file is authoritative.

### POG570 confirmatory hypotheses

Fixed order. Each one-sided α = 0.05. One sample per patient. Each hypothesis resamples its analysis set by patient 2,000 times. Seed 20261001. Percentiles 2.5 and 97.5. The random-generator call order is H1, H2, H3. The difference is the HostMix-TOO statistic minus the baseline statistic.

| Order | Hypothesis | Test |
| --- | --- | --- |
| H1 | In the at-risk set, the host-attraction rate of HostMix-TOO is lower than that of the baseline | Exact one-sided McNemar test on the host-attraction indicator. A favorable discordant pair is a biopsy that is a host-attraction error for the baseline only. p is the probability that a Binomial(n_disc, 0.5) variable is at least the observed favorable count. If n_disc = 0, p = 1 |
| H2 | In the evaluation set, the organ top-1 accuracy of HostMix-TOO is higher than that of the baseline | Exact one-sided McNemar test. A favorable discordant pair is a biopsy correct for HostMix-TOO only. Tested only if H1 is met |
| H3 | In the native-truth set, the top-1 accuracy of HostMix-TOO is not lower than that of the baseline by more than 10 percentage points | The 5th percentile of the same H3 bootstrap differences is greater than −0.10. Judged only if H2 is met |

H1, H2 and H3 all report the difference and a two-sided 95% confidence interval. H3 also reports the one-sided 5th percentile.

A prediction of NA is incorrect.

**Analysis sets**

- Evaluation set: the organ is neither excluded nor NA. n = 512.
- At-risk set: evaluation samples whose standard site has a non-empty native-organ set and whose true organ lies outside that set. n = 378.
- Native-truth set: the true organ lies inside that set. n = 91.

**Pre-specified descriptive outputs**

No multiplicity adjustment. Point estimates are reported.

- Site groups: liver, lymph_node, lung, soft_tissue, remainder. Remainder is any standard site other than those four. For each method: n, top-1, n_at_risk, host_rate.
- Tumor-content tertiles T1, T2, T3, from qcut with 3 bins on the 512 evaluation values, duplicates dropped, closed on the right. The confirmatory run does not recompute the boundaries. For each method: n, top-1, n_at_risk, host_rate.

Analyses that are not in the preregistration are reported separately as post hoc.

Translated from the Korean plan file committed at ba3c074e96e0cfd0124e54a24e08964849e7879e; the original file is authoritative.

### Auxiliary confirmatory hypotheses

**Primary family**

RNA-seq layer, one sample per patient, fixed order, each one-sided α = 0.05.

| Order | Hypothesis | Test |
| --- | --- | --- |
| AH1 | In the at-risk set, the host-attraction rate of HostMix-TOO is lower than that of the baseline | Exact one-sided McNemar test |
| AH2 | In the evaluation set, the top-1 accuracy of HostMix-TOO is higher than that of the baseline | Exact one-sided McNemar test. Tested only if AH1 is met |
| AH3 | In the native-truth set, the top-1 accuracy of HostMix-TOO is not lower than that of the baseline by more than 10 percentage points | One-sided 95% lower bound of the patient-bootstrap difference (2,000 replicates) > −0.10. Judged only if AH2 is met |

**Secondary family**

RNA-seq layer, one sample per patient, Holm correction, family-wise one-sided α = 0.05. The secondary family is tested regardless of the primary family.

| Hypothesis | Statement | Test |
| --- | --- | --- |
| AS1 | In the native-truth set, the top-1 accuracy of the gated model is higher than that of HostMix-TOO | Exact one-sided McNemar test |
| AS2 | In the pool-out at-risk set, the host-attraction rate of the 22-tissue model is lower than that of HostMix-TOO | Exact one-sided McNemar test |
| AS3 | In the evaluation set, the top-1 accuracy of the gated model is higher than that of HostMix-TOO | Exact one-sided McNemar test |

Error-rate control is separate for each family. The primary family uses a fixed order and the secondary family uses Holm, each at family-wise α = 0.05. The two families are not controlled together. The paper's primary claims rest only on the primary family.

Effect size for every hypothesis: the difference and a 95% patient-bootstrap interval (2,000 replicates, seed 20261001, percentile method), plus the discordant counts.

Exact one-sided McNemar p = P(Binomial(m, 0.5) ≥ favorable count). m = favorable + unfavorable. If m = 0, the test is not rejected.

**Descriptive analyses named before unlock**

No multiplicity adjustment.

The microarray layer repeats the same comparisons as effect sizes and intervals only.

A cohort meta-analysis uses cohorts with at least 10 at-risk samples and one sample per patient. d_i = (c_i − b_i) / n_i, where b_i is the number of patients who are a host-attraction error for the baseline only, c_i is the number who are a host-attraction error for HostMix-TOO only, and n_i is the number of at-risk patients. The variance is v_i = [(b_i + c_i) − (b_i − c_i)^2 / n_i] / n_i^2. If b_i or c_i is 0, 0.5 is added to both counts for the variance calculation only. The summary is a DerSimonian–Laird random-effects estimate with a 95% interval and I^2, separately by layer.

Analyses that are not in the preregistration are reported separately as post hoc.

## 10. Design power, synthetic labels, planted predictions and the independent checker

Design power for POG570 was computed before that confirmation (`results/stage3/power.json`) and is not a result of the confirmatory data. In each simulated data set every biopsy gave a discordant pair favoring the hypothesis with probability p_favor, a discordant pair against it with probability p_against, and a concordant pair otherwise; the one-sided exact McNemar test was applied at α = 0.05, with 10000 simulations per scenario and seed 20261001. The full scenario used the MET500 development proportions of the baseline and HostMix-TOO: for H1, at-risk biopsies with a host-attraction error under the baseline only (34 of 361) or under HostMix-TOO only (1 of 361); for H2, evaluation biopsies correct under HostMix-TOO only (29 of 437) or under the baseline only (12 of 437). The half scenario halved p_favor and kept p_against. These assumptions were preregistered before the POG570 features were computed.

| Hypothesis | Scenario | n | p_favor | p_against | α | Simulations | Power |
|---|---|---|---|---|---|---|---|
| H1 | MET500 proportions | 378 | 0.094 (34/361) | 0.003 (1/361) | 0.05 | 10000 | 1.000 |
| H1 | Half effect | 378 | 0.047 (17/361) | 0.003 (1/361) | 0.05 | 10000 | 0.996 |
| H2 | MET500 proportions | 512 | 0.066 (29/437) | 0.027 (12/437) | 0.05 | 10000 | 0.873 |
| H2 | Half effect | 512 | 0.033 (29/874) | 0.027 (12/437) | 0.05 | 10000 | 0.099 |

Design power for the auxiliary family was computed before unlock. The table below is that calculation. It is not a result of the confirmatory data.

| hypothesis | scenario | n | p_favor | p_against | alpha | power |
|---|---|---|---|---|---|---|
| AH1 | 1 | 427 | 0.204 | 0.000 | 0.05 | 1.000 |
| AH1 | 2 | 427 | 0.094 | 0.003 | 0.05 | 1.000 |
| AH2 | 1 | 729 | 0.119 | 0.045 | 0.05 | 1.000 |
| AH2 | 2 | 729 | 0.066 | 0.027 | 0.05 | 0.960 |
| AS1 | 1 | 315 | 0.077 | 0.011 | 0.017 | 0.978 |
| AS1 | 2 | 315 | 0.060 | 0.027 | 0.017 | 0.360 |
| AS2 | 1 | 71 | 0.084 | 0.002 | 0.017 | 0.507 |
| AS2 | 2 | 71 | 0.043 | 0.002 | 0.017 | 0.076 |
| AS3 | 1 | 729 | 0.037 | 0.018 | 0.017 | 0.493 |
| AS3 | 2 | 729 | 0.018 | 0.011 | 0.017 | 0.103 |

Before unlock, the confirmatory code was run on synthetic labels and on planted predictions with known expected results. An independent script that does not import the analysis code recomputed the auxiliary test statistics. That comparison has 274 rows and 0 mismatches. No separate independent checker was recorded for the POG570 confirmation.

## 11. Published classifiers

SCOPE and CUP-AI-Dx were run from the authors' public code and trained models. The input rules were committed before any cohort was run with either tool. CUP-AI-Dx uses the 791 features in the repository file `features_791.csv`. Cohort expression was converted to log2(TPM + 1), features absent from a cohort were set to 0, and each profile was then standardized across the 791 features, as in the repository's external test script. SCOPE uses its 17688-gene list. Features absent from a cohort were set to 0. Exon-union length normalization was applied only to the PRINCE counts, which had not already been divided by length; the other cohorts were passed as TPM or FPKM. The pre-specified rule passed inputs that are a constant multiple of RPKM within each sample unchanged and ran SCOPE once more on the doubled matrix to check scale dependence; the reported SCOPE calls are those of the unchanged run, and the doubled run is not the result in any table. The fraction of samples whose SCOPE top-1 class changed in the doubled run was 0 in POG570, IMvigor210, Anders and PRINCE and 0.008 (1 of 121) in DFCI melanoma.

| Cohort | Values provided | SCOPE input | Length normalization for SCOPE | Scale changed for the reported SCOPE run | CUP-AI-Dx input before standardization | Evidence (file: line) |
|---|---|---|---|---|---|---|
| MET500 | FPKM | FPKM as provided | No | No; no doubled run | FPKM converted to a sum of 10^6 per sample, then log2(TPM + 1) | `config/external_inputs_A8.tsv`: 2 |
| POG570 | TPM | TPM as provided | No | No; doubled run for the scale check only | log2(TPM + 1) | `config/external_inputs_A8.tsv`: 3 |
| IMvigor210 | log2(TPM + 1) | TPM after 2^v − 1, negative values set to 0 | No | No; doubled run for the scale check only | log2(TPM + 1) | `config/external_inputs_A8.tsv`: 4 |
| Anders | log2(TPM + 1) | TPM after 2^v − 1, negative values set to 0 | No | No; doubled run for the scale check only | log2(TPM + 1) | `config/external_inputs_A8.tsv`: 5 |
| DFCI melanoma | TPM | TPM as provided | No | No; doubled run for the scale check only | log2(TPM + 1) | `config/external_inputs_A8.tsv`: 6 |
| PRINCE | log2(upper-quartile normalized counts + 1) | 2^v − 1, negative values set to 0, divided by GENCODE v23 exon-union length in kb | Yes | No; doubled run for the scale check only | Length-divided values converted to a sum of 10^6 per sample, then log2(TPM + 1) | `config/external_inputs_A8.tsv`: 7 |
| GSE50760 | FPKM | FPKM as provided | No | No; no doubled run | FPKM converted to a sum of 10^6 per sample, then log2(TPM + 1) | `config/external_inputs_A8.tsv`: 8 |
| SU2C/PCF | FPKM | FPKM as provided | No | No; no doubled run | FPKM converted to a sum of 10^6 per sample, then log2(TPM + 1) | `config/external_inputs_A8.tsv`: 9 |

The scale-check rule is line 64 of `config/analysis_plan_A8.yaml`, and the CUP-AI-Dx preprocessing is line 55 of that file. For GSE50760, 18 at-risk biopsies outside the evaluation matrix were run once more with both tools under the same input rules, so that the at-risk set was complete; the SCOPE and CUP-AI-Dx outputs of both runs were kept and combined for scoring. SCOPE calls of a normal-tissue class, marked by the authors' normal-sample suffix, were counted as incorrect in the primary comparison. A sensitivity count that mapped those calls to the corresponding organ is reported in the main text.

The CUP-AI-Dx example was reproduced before the study cohorts were scored. The wrapper compares each prediction with the tumor type in the authors' metastatic metadata (`ExternalDataMeta.csv`, column tumor.type): 363 of 394 were correct (92.1%). The possible overlap between POG570 and the SCOPE development samples could not be checked, because the SCOPE development identifiers were not available to us. MET500 was the development cohort for HostMix-TOO and is not a confirmation.

## 12. Post hoc analyses

The ablation plan was committed before the three extra logistic models were fit. For each cohort it defined the gap between the baseline and HostMix-TOO host-attraction rates and the share of that gap attributable to the mixtures. The reading rule was that the mixtures would be considered the main source if all nine shares were at least 0.5. That rule was met. The ablation is post hoc relative to the preregistered hypotheses.

The missing-gene simulation restricted POG570 and MET500 profiles to the genes measured on each microarray platform and compared the change in top-1 accuracy with ten random removals of the same number of classifier genes. Models were not retrained.

The gene contribution in Figure S4 is the coefficient difference, esophagus logit minus bladder logit, multiplied by the model input: the rank-normal score for the baseline and the standardized score for HostMix-TOO. Positive values favor esophagus. Values are not comparable between the two models. The ten largest positive and the five most negative mean contributions are shown.

Esophagus calls were assigned to the side with the greater Spearman correlation to the adenocarcinoma or squamous centroid of the TCGA esophageal carcinoma class. Figure S5 uses one cohort order in every panel.


### Table S8. Post hoc ablation

Five logistic models compared after the confirmations: the baseline (pure tumors, unstandardized, C = 0.1); two standardized models trained on pure tumors (C = 0.03 and C = 0.15); the unstandardized mixture model, trained with the HostMix-TOO mixtures (C = 0.1); and HostMix-TOO (mixtures, standardized, C = 0.03). Top-1 accuracy uses the evaluation set, the host-attraction rate uses the at-risk set, and native-truth top-1 uses the native-truth set. This table is post hoc.

| Cohort | Model | n | Top-1 | n at risk | Host-attraction rate | n native-truth | Native-truth top-1 |
|---|---|---|---|---|---|---|---|
| MET500 | Baseline | 437 | 0.689 | 361 | 0.136 | 76 | 0.803 |
| MET500 | Standardized, pure (C = 0.03) | 437 | 0.680 | 361 | 0.122 | 76 | 0.803 |
| MET500 | Standardized, pure (C = 0.15) | 437 | 0.666 | 361 | 0.139 | 76 | 0.803 |
| MET500 | Unstandardized, mixtures (C = 0.1) | 437 | 0.728 | 361 | 0.044 | 76 | 0.789 |
| MET500 | HostMix-TOO | 437 | 0.728 | 361 | 0.044 | 76 | 0.789 |
| POG570 | Baseline | 512 | 0.682 | 378 | 0.228 | 91 | 0.868 |
| POG570 | Standardized, pure (C = 0.03) | 512 | 0.705 | 378 | 0.196 | 91 | 0.857 |
| POG570 | Standardized, pure (C = 0.15) | 512 | 0.670 | 378 | 0.222 | 91 | 0.868 |
| POG570 | Unstandardized, mixtures (C = 0.1) | 512 | 0.777 | 378 | 0.034 | 91 | 0.846 |
| POG570 | HostMix-TOO | 512 | 0.756 | 378 | 0.024 | 91 | 0.802 |
| Auxiliary RNA-seq | Baseline | 729 | 0.527 | 427 | 0.204 | 315 | 0.473 |
| Auxiliary RNA-seq | Standardized, pure (C = 0.03) | 729 | 0.494 | 427 | 0.176 | 315 | 0.429 |
| Auxiliary RNA-seq | Standardized, pure (C = 0.15) | 729 | 0.491 | 427 | 0.201 | 315 | 0.454 |
| Auxiliary RNA-seq | Unstandardized, mixtures (C = 0.1) | 729 | 0.706 | 427 | 0.082 | 315 | 0.727 |
| Auxiliary RNA-seq | HostMix-TOO | 729 | 0.698 | 427 | 0.056 | 315 | 0.730 |

### Table S9. Post hoc ablation in the microarray layer and in the missing-gene simulation

Top-1 accuracy in the microarray evaluation set (n = 536), host-attraction rate in the microarray at-risk set (n = 207), top-1 accuracy in the microarray native-truth set (n = 341) and number of esophagus calls in the evaluation set for the five models of the ablation (the baseline, the two standardized models trained on pure tumors, the unstandardized mixture model and HostMix-TOO; Table S8), and the mean loss of top-1 accuracy in POG570 and MET500 when the classifier genes missing from each microarray platform were removed and when as many randomly chosen classifier genes were removed (ten draws; four platforms, two cohorts). Shares are the fraction of HostMix-TOO's loss reproduced by each model. The plan, including the reading rule, was committed before these models were applied to microarray or gene-restricted profiles. All values are post hoc.

| Model | Top-1 (n = 536) | Host-attraction rate (n = 207) | Native-truth top-1 (n = 341) | Esophagus calls | Mean random-removal loss |
|---|---|---|---|---|---|
| Baseline | 0.724 | 0.256 | 0.871 | 30 | 0.006 |
| Standardized, pure (C = 0.03) | 0.360 | 0.208 | 0.446 | 183 | 0.049 |
| Standardized, pure (C = 0.15) | 0.332 | 0.251 | 0.408 | 191 | 0.049 |
| Unstandardized, mixtures (C = 0.1) | 0.776 | 0.092 | 0.891 | 11 | 0.006 |
| HostMix-TOO | 0.455 | 0.101 | 0.554 | 147 | 0.031 |

The shares of HostMix-TOO's top-1 loss on the microarray evaluation set were −0.194 for the unstandardized mixture model, 1.354 for standardized pure tumors with C = 0.03 and 1.458 for C = 0.15. The reading rule assigned this loss to standardization (category S). The shares of the extra loss under random gene removal were −0.028, 1.711 and 1.704, also category S. The gap in that loss between HostMix-TOO and the baseline was 0.025.

## 13. Supplementary figures

**Figure S1. Accuracy and host attraction of the locked classifiers**

Each cell is the top-1 organ accuracy (a) or the host-attraction rate (b) of one locked classifier. Columns are MET500, POG570, the six confirmatory RNA-seq cohorts, and the microarray cohorts. The MET500 column is exploratory. In the POG570 column, the baseline and HostMix-TOO cells are the two arms of the preregistered comparisons H1 and H2, and the other cells are pre-specified descriptive or exploratory; the status of every cell is given in Additional file 2. The auxiliary column is a post hoc aggregate of the confirmatory predictions on the six RNA-seq cohorts. The microarray column is a post hoc aggregate of the one-patient descriptive layer (top-1 on 536 evaluation biopsies; host-attraction rate on 207 at-risk biopsies). Shading runs from light to dark gray over 0 to 1 in a and over 0 to the largest rate in b; labels give the value to three decimals. n/c, not computed (hatched cells). Native-organ masking has a host-attraction rate of 0 by construction.

**Figure S2. Accuracy and host attraction in three analysis sets**

Top-1 accuracy (upper row) and a second metric (lower row) of the baseline, HostMix-TOO, SCOPE and CUP-AI-Dx in the native-truth set (a), the at-risk set (b) and the pool-out at-risk set (c), shown separately for MET500, POG570 and the auxiliary RNA-seq cohorts; n under each group is the size of that set. In b and c the lower row is the host-attraction rate. Host attraction is defined only for at-risk biopsies, so in a the lower row ("called another native organ") is the fraction of native-truth biopsies called a native organ of the biopsy site other than the true organ. A bar of height zero is labeled 0. This decomposition was made after both confirmations and is post hoc.

**Figure S3. Esophagus calls and the effect of missing genes**

**a** Fraction of biopsies whose true organ is not esophagus that were called esophagus, in MET500, POG570, the six auxiliary RNA-seq cohorts and the microarray cohorts; n is that denominator. For the microarray cohorts, the set has one sample per patient across the layer, as in Table 1 (n = 536). n/c, not computed: the esophagus fraction of linear deconvolution was not computed for MET500. **b** HostMix-TOO calls for the 183 GSE41258 native-truth colorectal tumors; 102 were called esophagus. **c** Change in top-1 accuracy after restricting POG570 and MET500 profiles to the genes of GPL96, GPL20769, GPL15659 or the FHCRC Agilent array; solid bars are POG570 and hatched bars MET500, colored by classifier. The horizontal mark is the mean change across ten random removals of the same number of classifier genes, and the vertical line spans the minimum to maximum of those ten changes. The missing-gene simulation is post hoc. Models were not retrained.

**Figure S4. Genes contributing to esophagus calls in IMvigor210**

Mean contribution of each gene to the difference between the esophagus and bladder logits for IMvigor210 bladder tumors in the native-truth set that were called esophagus: **a** baseline (n = 145), **b** HostMix-TOO (n = 65). The contribution is the coefficient difference multiplied by the model input, the rank-normal score for the baseline and the standardized score for HostMix-TOO, so values are not comparable between panels. Positive values favor esophagus. The ten largest positive and the five most negative mean contributions are shown. This analysis is post hoc and does not identify a mechanism.

**Figure S5. Assignment of esophagus calls to adenocarcinoma or squamous centroids**

Each row is one cohort, grouped as the six auxiliary RNA-seq cohorts, the other RNA-seq cohort (AURORA US) and the five microarray cohorts, in the same order in every panel: **a** baseline, **b** HostMix-TOO, **c** linear deconvolution. Each esophagus call is assigned to the adenocarcinoma or the squamous centroid, whichever has the greater Spearman correlation with the profile. The n at the end of a bar is that panel's number of esophagus calls, including n = 0 when the classifier made none; GSE50760 had none under any classifier. Of 147 HostMix-TOO esophagus calls in GSE41258, 145 were closer to the adenocarcinoma centroid and 2 to the squamous centroid; in IMvigor210, 69 of 105 were closer to the squamous centroid. This assignment is post hoc.

**Figure S6. Cohort-level differences in host-attraction rate**

Difference in host-attraction rate, HostMix-TOO minus baseline, in each cohort with at least 10 at-risk biopsies (at-risk n in parentheses), with Wald 95% intervals from the recorded variances. Diamonds show the DerSimonian–Laird random-effects estimate and its 95% interval for each layer (RNA-seq: seven cohorts, including the descriptive cohort AURORA US; microarray: four cohorts); the columns on the right give each difference with its 95% interval and the I² of each layer as a percentage. Cohort-level differences use each cohort's own one-sample-per-patient set, so patients profiled in both GSE74685 and FHCRC are counted in each cohort, and the at-risk n of GSE74685 (57) is larger than in Table 1 (49). FHCRC is absent because its at-risk count was below 10. The cohort differences and layer estimates were specified in advance as descriptive analyses; they do not replace the preregistered tests.


## 14. Repository correspondence

### Table S7. Correspondence between manuscript names and project files

Internal script and file names appear in this table and in the source notes of Sections 2, 3, 9, 10 and 11. Files under `results/` are produced by the scripts and are not in the public repository; `config/prereg_A6.md` is withheld from it because it lists POG570 patient identifiers with per-patient results (see `WITHHELD.tsv` there). Its English translation, `config/prereg_A6_en.md`, is in the repository, and the copy of the git history available to editors and reviewers contains its hash and a copy with these identifiers masked.

| Manuscript name | Internal name | Script | Output |
|---|---|---|---|
| Baseline | BASE-Z | scripts/08_train.py | results/stage2/models/BASE_Z.joblib |
| HostMix-TOO | SA-Z | scripts/08_train.py | results/stage2/models/SA_Z.joblib |
| Normal classes | NC-Z | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Site-specific mixtures | SC-Z | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Linear deconvolution | LD-Z | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Native-organ masking | M1-Z | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Gene removal 20% | IF20-Z | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Gene sets | BASE-K | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Gene sets + mixtures | SA-K | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Gene sets + host correction | V0-K | scripts/08_train.py | results/stage2/met500_overall.tsv |
| Gated model | SA-G | scripts/24_aux_confirm.py | results/stage7/confirm/tables/metrics.tsv |
| 22-tissue model | SA-pool22 | scripts/14_stage5_train.py | results/stage5/posthoc/pog_overall.tsv |
| 3-tissue model | SA-pool3 | scripts/14_stage5_train.py | results/stage5/posthoc/pog_overall.tsv |
| POG570 hypotheses | H1 H2 H3 | scripts/12_confirm.py | results/stage4/confirm/primary.tsv |
| Auxiliary hypotheses | AH AS | scripts/24_aux_confirm.py | results/stage7/confirm/tables/hypothesis_primary.tsv |
| Published classifiers | SCOPE CUP-AI-Dx | scripts/29_cup_wrapper_A8.py | results/stage8/external/tables/metrics.tsv |
| Ablation | PURE-Zs PURE-Zs-w MIX-Z0 | scripts/45_ablation_A12.py | results/stage9/ablation/metrics.tsv |
| Figure source | display names | scripts/50_figures_A13.py | manuscript/bmc/figure_source/ |
