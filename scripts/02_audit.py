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
import zipfile

import numpy as np
import pandas as pd

from loader import (
    POG570Locked,
    load_pog570_expression,
    pog570_expression_structure,
    read_met500_meta,
    read_phenotype,
    read_pog570_table,
    read_skcm_clinical,
    verify_pog570_lock_hashes,
    write_pog570_lock_hashes,
)
from mapping_rules import (
    DISEASE_TO_PROJECT,
    LEARNING_LABEL,
    MET500_COHORT_TO_ORGAN,
    POG570_COHORT_TO_ORGAN,
    PROJECT_TO_ORGAN,
    SITE_DICTIONARY,
    SITE_NATIVE,
    map_site,
)
from util import load_paths, now_iso


def freq_frame(series: pd.Series, full_limit: int = 50, top: int = 30) -> pd.DataFrame:
    values = series.where(series.notna(), "<NA>").astype(str)
    counts = values.value_counts(dropna=False)
    frame = counts.rename_axis("value").reset_index(name="n")
    frame["n_unique_total"] = len(counts)
    frame["truncated"] = len(counts) > full_limit
    if len(counts) > full_limit:
        frame = frame.head(top)
    return frame


def write_freq(folder: Path, name: str, series: pd.Series) -> None:
    freq_frame(series).to_csv(folder / f"{name}.tsv", sep="\t", index=False)


def main() -> None:
    paths = load_paths()
    audit = paths["results"] / "audit"
    audit.mkdir(parents=True, exist_ok=True)
    map_dir = paths["config"] / "mappings"
    map_dir.mkdir(parents=True, exist_ok=True)

    ph = read_phenotype()
    pd.Series(ph.columns, name="column").to_csv(audit / "phenotype_columns.tsv", sep="\t", index=False)
    for column in ph.columns:
        if column == "sample":
            continue
        write_freq(audit, f"phenotype_{column.replace(' ', '_')}", ph[column])

    tcga = ph.loc[ph["_study"] == "TCGA"].copy()
    tcga["patient"] = tcga["sample"].str.slice(0, 12)
    tcga["code"] = tcga["sample"].str.split("-").str[3].str[:2]
    tcga["project"] = tcga["primary disease or tissue"].map(DISEASE_TO_PROJECT)
    tcga["learning_label"] = tcga["project"].map(LEARNING_LABEL)
    learn_mask = ((tcga["code"] == "01") & (tcga["_sample_type"] == "Primary Tumor")) | (
        (tcga["code"] == "03")
        & (tcga["_sample_type"] == "Primary Blood Derived Cancer - Peripheral Blood")
        & (tcga["project"] == "LAML")
    )
    learn = tcga.loc[learn_mask]
    met = tcga.loc[tcga["code"] == "06"]
    normal = tcga.loc[tcga["code"] == "11"]
    recurrent_01 = tcga.loc[(tcga["code"] == "01") & (tcga["_sample_type"] != "Primary Tumor")]
    recurrent_01.to_csv(audit / "tcga_code01_not_primary.tsv", sep="\t", index=False)

    def project_counts(frame: pd.DataFrame) -> pd.DataFrame:
        counts = frame["learning_label"].fillna(frame["project"]).value_counts(dropna=False)
        return counts.rename_axis("project").reset_index(name="n")

    project_counts(learn).to_csv(audit / "tcga_learning_by_project.tsv", sep="\t", index=False)
    project_counts(met).to_csv(audit / "tcga_met06_by_project.tsv", sep="\t", index=False)
    project_counts(normal).to_csv(audit / "tcga_normal11_by_project.tsv", sep="\t", index=False)
    pd.crosstab(tcga["code"].fillna("<NA>"), tcga["_sample_type"].fillna("<NA>")).to_csv(
        audit / "tcga_code_by_sample_type.tsv", sep="\t"
    )
    summary = {
        "n_phenotype": int(len(ph)),
        "n_tcga": int(len(tcga)),
        "n_gtex": int((ph["_study"] == "GTEX").sum()),
        "n_target": int((ph["_study"] == "TARGET").sum()),
        "n_learning": int(len(learn)),
        "n_met06": int(len(met)),
        "n_normal11": int(len(normal)),
        "n_code01_not_primary": int(len(recurrent_01)),
        "learning_min_class": int(learn["learning_label"].value_counts().min()),
        "patient_code_max": int(tcga.groupby(["patient", "code"]).size().max()),
        "datetime": now_iso(),
    }
    (audit / "tcga_summary.json").write_text(json.dumps(summary, indent=2))

    gtex = ph.loc[ph["_study"] == "GTEX"].copy()
    gtex["donor"] = gtex["sample"].str.split("-").str[:2].str.join("-")
    gtex_table = (
        gtex.groupby("primary disease or tissue", dropna=False)
        .agg(n_samples=("sample", "size"), n_donors=("donor", "nunique"))
        .reset_index()
        .sort_values("n_samples", ascending=False)
    )
    gtex_table.to_csv(audit / "gtex_tissue_counts.tsv", sep="\t", index=False)

    meta = read_met500_meta()
    pd.Series(meta.columns, name="column").to_csv(audit / "met500_columns.tsv", sep="\t", index=False)
    for column in meta.columns:
        write_freq(audit, f"met500_{column.replace('.', '_')}", meta[column])
    tc = pd.to_numeric(meta["tc"], errors="coerce")
    quantiles = tc.quantile([0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]).rename_axis("quantile").reset_index(name="tc")
    quantiles.to_csv(audit / "met500_tc_quantiles.tsv", sep="\t", index=False)
    meta = meta.copy()
    meta["library"] = meta["Sample_id"].str.extract(r"-(capt|poly)-", expand=False)
    lib_counts = meta.groupby("sample_source")["library"].apply(lambda s: tuple(sorted(s.dropna())))
    lib_counts.value_counts().rename_axis("libraries").reset_index(name="n_sample_source").to_csv(
        audit / "met500_library_pairs.tsv", sep="\t", index=False
    )
    (audit / "met500_structure.json").write_text(
        json.dumps(
            {
                "n_rows": int(len(meta)),
                "n_sample_id": int(meta["Sample_id"].nunique()),
                "n_sample_source": int(meta["sample_source"].nunique()),
                "n_run_id_equals_sample_id": int((meta["run.id"] == meta["Sample_id"]).sum()),
                "tc_missing": int(tc.isna().sum()),
                "library_capt": int((meta["library"] == "capt").sum()),
                "library_poly": int((meta["library"] == "poly").sum()),
                "library_other": int(meta["library"].isna().sum()),
            },
            indent=2,
        )
    )

    # Lock hashes before any POG expression read. Verify on later runs.
    lock_path = paths["config"] / "locked_hashes.tsv"
    if not lock_path.exists():
        write_pog570_lock_hashes()
    verify_pog570_lock_hashes()
    locked = False
    try:
        load_pog570_expression(unlock=False)
    except POG570Locked:
        locked = True
    structure = pog570_expression_structure()
    s1 = read_pog570_table("s1")
    # ID-format check only: column ids versus PATIENT_ID, no expression values.
    import gzip

    with gzip.open(paths["data_raw"] / "POG570/POG570_TPM_expression.txt.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    expr_ids = set(header[1:])
    structure["n_expr_ids_in_patient_id"] = len(expr_ids & set(s1["PATIENT_ID"]))
    structure["lock_refused_without_unlock"] = locked
    (audit / "pog570_structure.json").write_text(json.dumps(structure, indent=2))

    s2 = read_pog570_table("s2")
    pd.Series(s1.columns, name="column").to_csv(audit / "pog570_s1_columns.tsv", sep="\t", index=False)
    pd.Series(s2.columns, name="column").to_csv(audit / "pog570_s2_columns.tsv", sep="\t", index=False)
    for column in s1.columns:
        if column in {"PATIENT_ID", "SAMPLE_ID_DNA", "SAMPLE_ID_RNA", "EGAD_ID"}:
            write_freq(audit, f"pog570_s1_{column}", s1[column])
            continue
        write_freq(audit, f"pog570_s1_{column}", s1[column])
    for column in s2.columns:
        if column in {"Patient_ID", "Therapy_start_date", "Therapy_end_or_biopsy_date"}:
            write_freq(audit, f"pog570_s2_{column}", s2[column])
            continue
        write_freq(audit, f"pog570_s2_{column}", s2[column])
    content = pd.to_numeric(s1["TUMOUR_CONTENT"], errors="coerce")
    content.quantile([0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]).rename_axis("quantile").reset_index(name="tumour_content").to_csv(
        audit / "pog570_tumour_content_quantiles.tsv", sep="\t", index=False
    )
    pd.crosstab(s1["BIOPSY_SITE"].fillna("<NA>"), s1["ANALYSIS_COHORT"].fillna("<NA>")).to_csv(
        audit / "pog570_biopsy_site_by_analysis_cohort.tsv", sep="\t"
    )
    pd.crosstab(s1["BIOPSY_COHORT"].fillna("<NA>"), s1["ANALYSIS_COHORT"].fillna("<NA>")).to_csv(
        audit / "pog570_biopsy_cohort_by_analysis_cohort.tsv", sep="\t"
    )

    # Mapping tables.
    disease_rows = [
        {"disease": disease, "tcga_project": project, "learning_label": LEARNING_LABEL[project]}
        for disease, project in DISEASE_TO_PROJECT.items()
    ]
    pd.DataFrame(disease_rows).to_csv(map_dir / "tcga_disease_to_project.tsv", sep="\t", index=False)
    organ_rows = [{"project": project, "organ": organ} for project, organ in PROJECT_TO_ORGAN.items()]
    pd.DataFrame(organ_rows).to_csv(map_dir / "project_to_organ.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"cohort": k, "organ": v} for k, v in MET500_COHORT_TO_ORGAN.items()]
    ).to_csv(map_dir / "met500_cohort_to_organ.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"analysis_cohort": k, "organ": v} for k, v in POG570_COHORT_TO_ORGAN.items()]
    ).to_csv(map_dir / "pog570_cohort_to_organ.tsv", sep="\t", index=False)
    site_rows = [{"source": a, "raw_value": b, "standard_site": c} for (a, b), c in SITE_DICTIONARY.items()]
    pd.DataFrame(site_rows).to_csv(map_dir / "site_dictionary.tsv", sep="\t", index=False)
    native_rows = []
    for site, (tissues, projects, organs) in SITE_NATIVE.items():
        native_rows.append(
            {
                "standard_site": site,
                "gtex_tissues": "|".join(tissues),
                "native_projects": "|".join(projects),
                "native_organs": "|".join(organs),
            }
        )
    pd.DataFrame(native_rows).to_csv(map_dir / "site_native.tsv", sep="\t", index=False)

    # Values present in the metadata that the dictionary does not list.
    met_missing = sorted(set(meta["biopsy_tissue"].dropna()) - {b for (src, b) in SITE_DICTIONARY if src.startswith("met500")})
    pog_missing = sorted(set(s1["BIOPSY_SITE"].dropna()) - {b for (src, b) in SITE_DICTIONARY if src.startswith("pog570")})
    (audit / "site_dictionary_unlisted.json").write_text(
        json.dumps({"met500_biopsy_tissue": met_missing, "pog570_biopsy_site": pog_missing}, indent=2)
    )

    skcm = read_skcm_clinical()
    pd.Series(skcm.columns, name="column").to_csv(audit / "skcm_columns.tsv", sep="\t", index=False)
    for column in ["sample_type", "sample_type_id", "tumor_tissue_site", "melanoma_origin_skin_anatomic_site"]:
        write_freq(audit, f"skcm_{column}", skcm[column])

    hpa_report = audit_hpa(paths["data_raw"] / "hpa", audit)
    (audit / "hpa_summary.json").write_text(json.dumps(hpa_report, indent=2))
    print("audit done", summary["n_learning"], "learning", structure["n_genes"], "pog genes")


def audit_hpa(folder: Path, audit: Path) -> dict:
    report = {}
    for name in ["rna_tissue_consensus.tsv.zip", "rna_tissue_hpa.tsv.zip", "rna_tissue_hpa_samples.tsv.zip"]:
        path = folder / name
        report[name] = {"exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else 0}
        if not path.exists():
            continue
        with zipfile.ZipFile(path) as zf:
            inner = zf.namelist()[0]
            with zf.open(inner) as handle:
                header = handle.readline().decode().rstrip("\n").split("\t")
            report[name]["inner"] = inner
            report[name]["columns"] = header
        if name.endswith("samples.tsv.zip"):
            with zipfile.ZipFile(path) as zf:
                text = zf.read(zf.namelist()[0]).decode()
            frame = pd.read_csv(pd.io.common.StringIO(text), sep="\t", dtype=str)
            tissues = sorted(frame["Tissue"].dropna().unique())
            report[name]["n_rows"] = int(len(frame))
            report[name]["tissues"] = tissues
            report[name]["has_lymph_node"] = any("lymph" in t.lower() for t in tissues)
            report[name]["has_bone_marrow"] = any("bone marrow" in t.lower() for t in tissues)
        else:
            tissue_col = "Tissue"
            counts = {}
            with zipfile.ZipFile(path) as zf:
                with zf.open(zf.namelist()[0]) as handle:
                    header = handle.readline().decode().rstrip("\n").split("\t")
                    tissue_i = header.index(tissue_col)
                    for line in handle:
                        tissue = line.decode().rstrip("\n").split("\t")[tissue_i]
                        counts[tissue] = counts.get(tissue, 0) + 1
            tissues = sorted(counts)
            report[name]["n_rows"] = int(sum(counts.values()))
            report[name]["tissues"] = tissues
            report[name]["has_lymph_node"] = any(t.lower() == "lymph node" for t in tissues)
            report[name]["has_bone_marrow"] = any(t.lower() == "bone marrow" for t in tissues)
            pd.Series(counts, name="n").rename_axis("tissue").reset_index().to_csv(
                audit / f"hpa_{name.replace('.tsv.zip', '')}_tissues.tsv", sep="\t", index=False
            )
    return report


if __name__ == "__main__":
    main()
