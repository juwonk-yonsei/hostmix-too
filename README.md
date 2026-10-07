# HostMix-TOO

Code, trained models and public-cohort predictions for the manuscript on host attraction in tissue-of-origin
classification of metastatic biopsies and its reduction by host-tissue mixture augmentation (HostMix-TOO).

Repository: https://github.com/juwonk-yonsei/hostmix-too. Archived version (Zenodo, all versions):
https://doi.org/10.5281/zenodo.23205336.

## Contents

| Folder | Contents |
| --- | --- |
| `release/` | The reuse package: trained models (baseline, HostMix-TOO and the post hoc unstandardized mixture model) in scikit-learn format and as portable NumPy archives, the host-pool profiles, the site dictionary, sample-level predictions for the public cohorts, environment pins and a data manifest with source URLs and hashes. Start with `release/README.md` and `release/models/MODEL_CARD.md`. |
| `scripts/` | Every analysis script of the project, numbered in the order they were written and run. |
| `config/` | Analysis plans, internal preregistration plans, mappings and split files. |
| `manuscript/` | Manuscript sources, figures, figure source data, additional files and the number-check files. |
| `logs/` | Environment and resource records of the runs. |

No expression data are included. TCGA, GTEx, MET500, the cBioPortal and GEO cohorts and POG570 are available from
their sources, listed with URLs and SHA-256 hashes in `release/data_manifest.tsv`.

## What this snapshot leaves out

This repository is a snapshot of the project repository without its git history. Files that list POG570
patient identifiers with their tumor labels or per-patient results, the internal work reports and the git history
bundle are left out; `WITHHELD.tsv` lists each of them with its SHA-256 and the reason. Manuscript drafts and
submission notes are not included either. A copy of the git history without the contents of the withheld files,
in which every commit hash, including those of the internal preregistration commits listed in Additional file 3
of the manuscript, can be verified, is available to editors and reviewers on request. `SNAPSHOT.txt` gives the
project commit the snapshot was made from.

## Licenses

Code: MIT (`LICENSE`). Trained models, host-pool profiles and predictions for the public cohorts other than
cBioPortal studies: CC BY 4.0. Predictions derived from cBioPortal studies: ODbL 1.0. Details are in
`LICENSE-DATA.md`.

## Citation

The article is under review; its DOI will be added here on publication. Until then, please cite this archive
(https://doi.org/10.5281/zenodo.23205336; `CITATION.cff`).
