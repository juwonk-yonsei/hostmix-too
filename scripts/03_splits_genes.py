#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import json

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from genes import build_common_genes, load_gene_sets, matrix_gene_ids
from loader import read_met500_meta, read_phenotype
from mapping_rules import DISEASE_TO_PROJECT, LEARNING_LABEL
from util import SEED, assert_disjoint, load_paths


def learning_and_met(ph: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tcga = ph.loc[ph["_study"] == "TCGA"].copy()
    tcga["patient"] = tcga["sample"].str.slice(0, 12)
    tcga["code"] = tcga["sample"].str.split("-").str[3].str[:2]
    tcga["project"] = tcga["primary disease or tissue"].map(DISEASE_TO_PROJECT)
    tcga["learning_label"] = tcga["project"].map(LEARNING_LABEL)
    learn = tcga.loc[
        ((tcga["code"] == "01") & (tcga["_sample_type"] == "Primary Tumor"))
        | (
            (tcga["code"] == "03")
            & (tcga["_sample_type"] == "Primary Blood Derived Cancer - Peripheral Blood")
            & (tcga["project"] == "LAML")
        )
    ].copy()
    met = tcga.loc[tcga["code"] == "06"].copy()
    if learn["learning_label"].isna().any():
        missing = sorted(learn.loc[learn["learning_label"].isna(), "primary disease or tissue"].unique())
        raise SystemExit(f"Unmapped TCGA disease in the learning set: {missing}")
    if learn["patient"].duplicated().any():
        raise SystemExit("More than one learning sample for a patient; the dedup rule was expected to be a no-op.")
    return learn, met


def main() -> None:
    paths = load_paths()
    ph = read_phenotype()
    learn, met = learning_and_met(ph)
    train_patients, test_patients = train_test_split(
        learn["patient"].to_numpy(),
        test_size=0.2,
        stratify=learn["learning_label"].to_numpy(),
        random_state=SEED,
    )
    train_patients = set(train_patients)
    test_patients = set(test_patients)
    assert_disjoint(train_patients, test_patients, "TCGA train/test patients")

    learn["split"] = np.where(learn["patient"].isin(train_patients), "TCGA-train", "TCGA-test")
    met["split"] = "TCGA-met"
    met["exclude_from_eval"] = met["patient"].isin(train_patients)
    met["exclude_reason"] = np.where(met["exclude_from_eval"], "primary_in_TCGA-train", "")
    columns = [
        "sample",
        "patient",
        "code",
        "project",
        "learning_label",
        "_sample_type",
        "primary disease or tissue",
        "split",
    ]
    learn_out = learn[columns].copy()
    met_out = met[columns + ["exclude_from_eval", "exclude_reason"]].copy()
    learn_out["exclude_from_eval"] = False
    learn_out["exclude_reason"] = ""
    split = pd.concat([learn_out, met_out], ignore_index=True)
    split_path = paths["config"] / "split_tcga.tsv"
    split.to_csv(split_path, sep="\t", index=False)
    assert_disjoint(
        split.loc[split["split"] == "TCGA-train", "patient"],
        split.loc[split["split"] == "TCGA-test", "patient"],
        "written TCGA split",
    )

    gtex_all = ph.loc[ph["_study"] == "GTEX"].copy()
    # K-562 cell-line ids are not GTEX donor barcodes. They are excluded from the donor split
    # and listed, rather than assigned an invented donor id.
    non_donor = gtex_all.loc[~gtex_all["sample"].str.startswith("GTEX-"), ["sample", "primary disease or tissue", "_sample_type"]]
    non_donor.to_csv(paths["results"] / "audit/gtex_nondonor_ids.tsv", sep="\t", index=False)
    gtex = gtex_all.loc[gtex_all["sample"].str.startswith("GTEX-")].copy()
    gtex["donor"] = gtex["sample"].str.split("-").str[:2].str.join("-")
    rng = np.random.default_rng(SEED)
    rows = []
    for tissue in sorted(gtex["primary disease or tissue"].dropna().unique()):
        sub = gtex.loc[gtex["primary disease or tissue"] == tissue]
        donors = np.array(sorted(sub["donor"].unique()))
        rng.shuffle(donors)
        n_ref = len(donors) // 2
        ref = set(donors[:n_ref])
        for row in sub.itertuples(index=False):
            rows.append(
                {
                    "sample": row.sample,
                    "donor": row.donor,
                    "tissue": tissue,
                    "split": "GTEx-ref" if row.donor in ref else "GTEx-sim",
                }
            )
    gtex_split = pd.DataFrame(rows)
    gtex_path = paths["config"] / "split_gtex.tsv"
    gtex_split.to_csv(gtex_path, sep="\t", index=False)

    # MET500 evaluation-unit table. One library per sample_source; poly over capt.
    meta = read_met500_meta()
    meta["library"] = meta["Sample_id"].str.extract(r"-(capt|poly)-", expand=False)
    chosen = []
    for source, sub in meta.groupby("sample_source", sort=True):
        poly = sub.loc[sub["library"] == "poly"]
        pick = poly.iloc[0] if len(poly) else sub.sort_values("Sample_id").iloc[0]
        chosen.append(pick)
    met_eval = pd.DataFrame(chosen)
    met_eval.to_csv(paths["config"] / "met500_eval_samples.tsv", sep="\t", index=False)

    toil_path = paths["data_raw"] / "toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    met_path = paths["data_raw"] / "met500/M.mx.txt.gz"
    if not toil_path.exists() or not met_path.exists():
        print("splits written; expression files not both present, gene set skipped")
        return
    print("scanning gene ids")
    common = build_common_genes(matrix_gene_ids(toil_path), matrix_gene_ids(met_path))
    genes = common["genes"]
    gene_dir = paths["data_processed"] / "genes"
    gene_dir.mkdir(parents=True, exist_ok=True)
    pd.Series(genes, name="symbol").to_csv(gene_dir / "G_symbols.txt", index=False, header=False)
    names, members, n_raw = load_gene_sets(genes)
    payload = {k: v for k, v in common.items() if k not in {"genes", "id_to_symbol", "unique_base"}}
    payload["n_G"] = len(genes)
    payload["n_gmt_raw"] = n_raw
    payload["n_sets_kept"] = len(names)
    (gene_dir / "G_summary.json").write_text(json.dumps(payload, indent=2))
    # Member indices are stored with the set names for the KS step.
    packed = paths["data_processed"] / "genes/gene_sets.npz"
    offsets = np.zeros(len(members) + 1, dtype=np.int64)
    chunks = [np.asarray(m, dtype=np.int32) for m in members]
    for i, chunk in enumerate(chunks):
        offsets[i + 1] = offsets[i] + chunk.size
    indices = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int32)
    np.savez(packed, offsets=offsets, indices=indices, names=np.array(names, dtype=object))
    print(json.dumps(payload, indent=2))
    print("train", int((split["split"] == "TCGA-train").sum()), "test", int((split["split"] == "TCGA-test").sum()))
    print("met", int((split["split"] == "TCGA-met").sum()), "met excluded", int(met["exclude_from_eval"].sum()))
    print("met500 eval", len(met_eval))


if __name__ == "__main__":
    main()
