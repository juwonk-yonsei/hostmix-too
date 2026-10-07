#!/usr/bin/env python3
"""Write config/analysis_plan_A1.yaml from the gene summary and mapping-file hashes.

Run this after scripts/03_splits_genes.py and before any accuracy calculation.
"""
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import json

from util import load_paths, now_iso, sha256_file


def main() -> None:
    paths = load_paths()
    summary = json.loads((paths["data_processed"] / "genes/G_summary.json").read_text())
    map_dir = paths["config"] / "mappings"
    hashes = []
    for path in sorted(map_dir.glob("*.tsv")):
        hashes.append(f"    {path.name}: {sha256_file(path)}")
    split_hashes = []
    for name in ("split_tcga.tsv", "split_gtex.tsv", "met500_eval_samples.tsv", "locked_hashes.tsv"):
        path = paths["config"] / name
        split_hashes.append(f"    {name}: {sha256_file(path)}")
    text = f"""# Analysis plan A1. Written {now_iso()} before section 10.
# Do not edit after MET500 or simulation metrics have been computed.
seed: 20261001
n_threads: 16

common_genes:
  n_G: {summary['n_G']}
  n_sets_kept: {summary['n_sets_kept']}
  n_gmt_raw: {summary['n_gmt_raw']}
  min_set_members_in_G: 5
  symbol_rule: >
    Toil and MET500 versioned Ensembl ids map through gencode v23 probemap symbols.
    A symbol is kept only if HGNC locus_group is protein-coding gene.
    POG570 unversioned Ensembl ids match the probemap after stripping the version.
    Unversioned ids with more than one symbol are dropped.
    Duplicate symbols in a matrix are summed on the linear scale.

classifiers:
  B0:
    features: within-sample rank-normal z on G
    gene_filter: variance ddof=1 on TCGA-train, top 5000, ties broken by smaller gene index
    gene_filter_scope: once on all TCGA-train, not re-fit inside CV folds
    model: multinomial logistic regression
    penalty: l2
    class_weight: balanced
    solver: lbfgs
    max_iter: 5000
    C_grid: [0.01, 0.1, 1.0, 10.0]
    C_selection: TCGA-train stratified 5-fold CV macro-F1, labels = classes present in the validation fold
    C_tie_break: smaller C
  B1:
    features: ONCOfind KS scores for gene sets with at least 5 members in G
    scaling: StandardScaler mean and sd fit inside each CV training fold; final scaler fit on all TCGA-train
    model: same logistic regression settings as B0

adjustments:
  M1: set organ probabilities in the biopsy site native organ set to 0 and renormalize. If the remainder sums to 0, the prediction is NA and counts as incorrect.
  V0: ES_corrected = (ES - 0.3 * GTEx-ref mean ES) / 0.7, then apply the fitted B1 scaler and model. Sites with no GTEx reference are not corrected.
  V0_multi_tissue: unweighted mean of the per-tissue GTEx-ref mean vectors.
  combinations: [B0, B0+M1, B1, B1+M1, B1+V0]

met500_eval_unit:
  rule: one library per sample_source
  library_field: token capt or poly inside Sample_id. The column sample_source is the specimen id, not the chemistry.
  preference: keep poly when both capt and poly exist; otherwise keep the only library
  bootstrap_unit: sample_source. No patient id is present, so patients are not grouped beyond sample_source.
  include_test_column: true
  unmapped_cohorts_excluded: [MISC, SECR]

metrics:
  top1: predicted label equals truth. NA predictions are incorrect.
  top3: truth is among the three highest probabilities. NA rows miss.
  macro_f1: sklearn f1_score average=macro, labels = sorted classes present in y_true, zero_division=0
  organ_probability: sum of TCGA project probabilities inside each organ group
  tie_break: numpy argmax, class order is sorted
  host_error: predicted organ is in the biopsy site native organ set and the true organ is not
  H: n_host_error / n_error. Undefined when n_error is 0.
  at_risk: native organ set is non-empty and the true organ is outside it
  host_rate: n_host_error / n_at_risk
  bootstrap: 2000 resamples of the evaluation unit with replacement
  ci: percentile 2.5 and 97.5 of the bootstrap statistics
  H_undefined_replicates: omitted from the H percentile and counted
  tc_tertiles: qcut on the mapped MET500 evaluation set, duplicates dropped
  library_strata: capt versus poly token of the kept library

simulation:
  n_tumors: 1000
  tumor_source: TCGA-test
  stratify: learning_label
  tissues: ["Liver", "Lung", "Adipose - Subcutaneous"]
  rho: [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2]
  rho_meaning: RNA fraction, not cell fraction
  mix: renormalize each profile on G to sum 1e6, then x = rho * tumor + (1-rho) * normal
  normal_draw: with replacement from GTEx-sim of that tissue
  rng_order: tissues in the list above, rho in the list above
  recompute: z and KS after mixing, no refit
  host_native_organs:
    Liver: [Liver, Biliary]
    Lung: [Lung, Mesothelioma]
    "Adipose - Subcutaneous": [Sarcoma]
  host_pull_rate_denominator: tumors whose true organ is outside that tissue native set

ks_implementation_rule:
  compare_n_pairs: 1000
  include: high-zero MET500 samples and TCGA samples
  accept_fast_if: max absolute score difference <= 1e-6 and sign mismatches = 0
  otherwise: original GetES with 16 processes

pog570:
  diagnosis_field: ANALYSIS_COHORT
  biopsy_field: BIOPSY_SITE
  expression: not used for risk counts
  unmapped_diagnosis: not counted as at risk

skcm:
  clinical_file: data/raw/xena/SKCM_clinicalMatrix
  join: exact 15-character sample id
  site_field: tumor_tissue_site

file_sha256:
{chr(10).join(hashes)}
{chr(10).join(split_hashes)}
"""
    dest = paths["config"] / "analysis_plan_A1.yaml"
    dest.write_text(text)
    print("wrote", dest)


if __name__ == "__main__":
    main()
