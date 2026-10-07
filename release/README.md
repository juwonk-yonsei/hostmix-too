# Release skeleton

This directory is the reuse package of the public repository https://github.com/juwonk-yonsei/hostmix-too, archived at Zenodo (https://doi.org/10.5281/zenodo.23205336).

## Licenses

- Code: MIT (`LICENSE`)
- Trained models and the host-pool profile: CC BY 4.0
- Predictions derived from cBioPortal studies: ODbL, stated in `predictions/ODBL_NOTICE.txt`
- Predictions for the other public cohorts: CC BY 4.0

## Git history

The public repository is a snapshot of the project repository without its git history. The history, which contains the preregistration commits listed in Additional file 3 of the manuscript, also contains files that list POG570 patient identifiers with their tumor labels or per-patient results and is therefore not public. A copy of the history without the contents of these files, in which every commit hash can be verified, is available to editors and reviewers on request. Files left out of the snapshot are listed with their SHA-256 in `WITHHELD.tsv` at the top of the repository.

## Reproduce

1. Set `project_root` and `figure_dir` in `config.yaml`. The file has no thread-count key.
2. Run the scripts in `scripts/` in numeric order, from the project root, with the environment pins in `env_pins.txt`.
3. Figure scripts `34_figures_A8.py` and `38_figures_A9.py` write to `figure_dir`. They do not retrain models.
4. Locked confirmation models are `models/BASE_Z.joblib` and `models/SA_Z.joblib`. The unstandardized mixture model of the post hoc ablation is `models/MIX_Z0_posthoc.joblib`; it was not evaluated in a preregistered confirmation. See `models/MODEL_CARD.md`.
5. Host-pool profiles are `host_pool/`. Read `host_pool/CITATION.txt` before reuse.
6. The site dictionary is `site_dictionary/mapping_rules.py`.
7. Sample-level predictions for public cohorts are under `predictions/`. cBioPortal-derived folders include `ODBL_NOTICE.txt`. POG570 sample-level predictions are not included. Aggregate POG570 metrics are only in the manuscript and in Additional file 2.

## Not included

POG570 expression, POG570 sample-level predictions and the POG570 label file. TCGA and GTEx expression matrices. Those stay at their sources under the terms quoted in the manuscript. The git history is not included (see "Git history").
