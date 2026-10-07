#!/usr/bin/env python3
"""Keep prad_fhcrc Agilent samples whose IDs are not already in GSE74685.

The public file is gene-level, so the probe-with-highest-mean rule is not applied.
Shared IDs stay with GSE74685, which keeps the sample-channel intensities.
Does not apply a classifier.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("auxgeo", ROOT / "20_aux_geo.py")
geo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geo)
aux = geo.aux

sys.path.insert(0, str(ROOT))
from stage2_common import load_gene_pack
from stage5_lib import load_paths

SITES = {
    "LN": "lymph_node",
    "LUNG": "lung",
    "LIVER": "liver",
    "BONE": "bone",
    "PROSTATE": "prostate",
    "BLADDER": "bladder",
    "ADRENAL": "adrenal",
    "SKIN": "skin",
    "PERITONEUM": "peritoneum",
    "PERITONEAL": "peritoneum",
    "KIDNEY": "kidney",
    "RENAL": "kidney",
    "BLADDER_NECK": "other",
    "PELVIC_MASS": "other",
    "RETROPERITONEAL": "other",
    "SPLEEN": "other",
    "SCROTUM": "other",
    "APPENDIX": "other",
}


def gse74685_ids(paths) -> set[str]:
    header = geo.header_frame(Path(paths["results"]) / "stage5" / "aux" / "geo_headers" / "GSE74685_header.txt")
    return set(header["source_name_ch2"].astype(str))


def main() -> None:
    paths = load_paths()
    expr_path = Path(paths["data_processed"]) / "aux" / "downloads" / "prad_fhcrc_data_mrna_agilent_microarray.txt"
    clinical_path = Path(paths["results"]) / "stage5" / "aux" / "clinical" / "prad_fhcrc_sample.txt"
    matrix = aux.read_cbio_expression(expr_path)
    shared = gse74685_ids(paths)
    keep = [col for col in matrix.columns if col not in shared]
    clinical = pd.read_csv(clinical_path, sep="\t", comment="#", dtype=str)
    clinical = clinical.drop_duplicates("SAMPLE_ID").set_index("SAMPLE_ID")
    rows = []
    for sample in keep:
        meta = clinical.loc[sample] if sample in clinical.index else None
        raw = "" if meta is None else str(meta["TISSUE_SITE"])
        if raw in {"", "nan"}:
            organ, site, rule = "exclude", "", "3"
        else:
            organ, site, rule = "Prostate", SITES.get(raw, "other"), "cohort"
        patient = sample if meta is None else str(meta["PATIENT_ID"])
        rows.append({
            "sample_id": sample, "patient_id": patient, "raw_site": raw,
            "standard_site": site, "organ": organ, "organ_rule": rule,
            "library": "cBio data_mrna_agilent_microarray.txt", "lcm": False,
            "overlap_note": "sample id absent from GSE74685 source_name_ch2",
        })
    samples = pd.DataFrame(rows)
    genes, b0_idx, offsets, indices, _names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    info = aux.finish_cohort(
        "prad_fhcrc", samples, matrix[keep], "microarray", pd.Series(dtype=float),
        genes, b0, offsets, indices, paths, "microarray",
        "gene-level Agilent file; probe-max rule not applied because probe rows are not in the cBio download; values min>0;  GSE74685-shared sample ids removed",
        None,
    )
    table = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / "prad_fhcrc" / "samples.tsv", sep="\t")
    log_path = Path(paths["results"]) / "stage5" / "aux" / "screening_log.tsv"
    log = pd.read_csv(log_path, sep="\t", dtype=str)
    note = (
        f"최종 판정: Agilent gene-level file, not RNA-seq; {len(shared & set(matrix.columns))} sample ids also in GSE74685 were removed; "
        f"kept {len(table)}; n_at_risk={int(table['in_risk'].sum())}"
    )
    if not (log["memo"] == note).any():
        base = log.loc[log["id"] == "prad_fhcrc"].iloc[0]
        extra = pd.DataFrame([{
            "path": base["path"], "query": base["query"], "id": "prad_fhcrc", "title": base["title"],
            "n": str(len(table)), "platform": "microarray", "decision": "포함",
            "failed_criterion": "", "memo": note,
        }])
        pd.concat([log, extra], ignore_index=True).to_csv(log_path, sep="\t", index=False)
    site_path = Path(paths["config"]) / "mappings" / "site_dictionary_aux.tsv"
    sites = pd.read_csv(site_path, sep="\t", dtype=str)
    add = pd.DataFrame([{"source": "prad_fhcrc", "raw_value": raw, "standard_site": standard} for raw, standard in SITES.items()])
    pd.concat([sites, add], ignore_index=True).drop_duplicates(["source", "raw_value"], keep="last").to_csv(site_path, sep="\t", index=False)
    license_text = aux.get_bytes("https://raw.githubusercontent.com/cBioPortal/datahub/master/public/prad_fhcrc/LICENSE").decode().strip()
    downloads = pd.read_csv(Path(paths["config"]) / "aux_downloads_A5.tsv", sep="\t", dtype=str)
    row = pd.DataFrame([{
        "cohort": "prad_fhcrc", "file": expr_path.name,
        "url": "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/prad_fhcrc/data_mrna_agilent_microarray.txt",
        "bytes": str(expr_path.stat().st_size), "sha256": aux.sha256_file(expr_path),
        "license": license_text,
    }])
    pd.concat([downloads, row], ignore_index=True).drop_duplicates(["cohort", "file"], keep="last").to_csv(
        Path(paths["config"]) / "aux_downloads_A5.tsv", sep="\t", index=False,
    )
    features = pd.read_csv(Path(paths["config"]) / "aux_features_A5.tsv", sep="\t", dtype=str)
    pd.concat([features, pd.DataFrame([info]).astype(str)], ignore_index=True).drop_duplicates("cohort", keep="last").to_csv(
        Path(paths["config"]) / "aux_features_A5.tsv", sep="\t", index=False,
    )
    print("prad_fhcrc kept", len(table), "risk", int(table["in_risk"].sum()), flush=True)


if __name__ == "__main__":
    main()
