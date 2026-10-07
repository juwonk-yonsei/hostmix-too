#!/usr/bin/env python3
"""Choose LN-proxy and BM-proxy by Spearman correlation of HPA tissue nTPM with GTEx-ref means.

Labels are not used. The tie-break is the order of CANDIDATE_TISSUES.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from stage2_common import (
    CANDIDATE_TISSUES,
    HPA_TARGETS,
    gtex_ref_ids,
    load_gene_pack,
    load_tpm,
    sum_1e6,
)
from util import load_paths, now_iso


def hpa_profiles(path: Path, genes: list[str]) -> tuple[dict[str, np.ndarray], dict]:
    frame = pd.read_csv(path, sep="\t", compression="zip")
    tissues = sorted(frame["Tissue"].astype(str).unique())
    wanted = {}
    for target in HPA_TARGETS:
        hits = [tissue for tissue in tissues if tissue.lower() == target]
        wanted[target] = hits
    gene_index = {gene: i for i, gene in enumerate(genes)}
    profiles = {}
    used = {}
    for target, hits in wanted.items():
        if len(hits) != 1:
            raise SystemExit(f"HPA tissue match for {target!r} is {hits}")
        # Duplicate Gene name rows are summed on the nTPM scale, then log2(sum + 1).
        sub = frame.loc[frame["Tissue"] == hits[0], ["Gene name", "nTPM"]].copy()
        sub["nTPM"] = pd.to_numeric(sub["nTPM"], errors="coerce")
        summed = sub.groupby("Gene name", sort=False)["nTPM"].sum()
        vec = np.full(len(genes), np.nan, dtype=np.float64)
        n_hit = 0
        for name, value in summed.items():
            idx = gene_index.get(str(name))
            if idx is None or not np.isfinite(value):
                continue
            vec[idx] = np.log2(float(value) + 1.0)
            n_hit += 1
        profiles[target] = vec
        used[target] = {"hpa_tissue": hits[0], "n_genes_in_G": n_hit, "duplicate_rule": "sum nTPM"}
    return profiles, {"hpa_tissues_all": tissues, "used": used}


def main() -> None:
    paths = load_paths()
    genes, _, _, _, _ = load_gene_pack(paths)
    profiles, meta = hpa_profiles(paths["data_raw"] / "hpa/rna_tissue_consensus.tsv.zip", genes)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    ref_ids = gtex_ref_ids(paths, CANDIDATE_TISSUES)
    means = {}
    for tissue, ids in ref_ids.items():
        if not ids:
            raise SystemExit(f"No GTEx-ref samples for {tissue}")
        block = sum_1e6(tpm[[row[s] for s in ids if s in row]])
        missing = [s for s in ids if s not in row]
        if missing:
            raise SystemExit(f"{tissue} missing from Toil: {len(missing)}")
        means[tissue] = np.log2(block.astype(np.float64) + 1.0).mean(axis=0)
        print(tissue, "ref", len(ids), flush=True)
    rows = []
    for target, hpa in profiles.items():
        for tissue in CANDIDATE_TISSUES:
            mask = np.isfinite(hpa)
            corr = float(spearmanr(hpa[mask], means[tissue][mask]).statistic)
            rows.append({
                "hpa_query": target,
                "hpa_tissue": meta["used"][target]["hpa_tissue"],
                "gtex_tissue": tissue,
                "n_genes": int(mask.sum()),
                "n_gtex_ref": len(ref_ids[tissue]),
                "spearman": corr,
            })
    table = pd.DataFrame(rows)
    out_dir = paths["results"] / "stage2"
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / "proxy_correlation.tsv", sep="\t", index=False)
    chosen = {}
    for target, sub in table.groupby("hpa_query", sort=False):
        best = float(sub["spearman"].max())
        tied = sub.loc[np.isclose(sub["spearman"], best)]
        order = {tissue: i for i, tissue in enumerate(CANDIDATE_TISSUES)}
        pick = sorted(tied["gtex_tissue"], key=lambda tissue: order[tissue])[0]
        chosen[target] = {"gtex_tissue": pick, "spearman": best, "n_tied": int(len(tied))}
    payload = {
        "datetime": now_iso(),
        "hpa_file": "data/raw/hpa/rna_tissue_consensus.tsv.zip",
        "transform": "log2(value + 1)",
        "gtex_value": "mean log2(TPM + 1) of GTEx-ref, TPM summed to 1e6 on G",
        "tie_break": "earlier entry in CANDIDATE_TISSUES",
        "ln_proxy": chosen["lymph node"]["gtex_tissue"],
        "bm_proxy": chosen["bone marrow"]["gtex_tissue"],
        "chosen": chosen,
        "meta": meta["used"],
    }
    (out_dir / "proxy_choice.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
