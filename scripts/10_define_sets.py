#!/usr/bin/env python3
"""Stage-3 metadata: eval labels, set sizes, tumour-content edges, power, model hashes.

Does not read POG570 expression values.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import binom, binomtest

from loader import read_pog570_table
from stage3_rules import (
    BETA_CUT,
    COHORT_MARGINAL_N,
    FAKE_ORGAN_WEIGHTS,
    attach_site,
    label_table,
    set_flags,
)
from util import load_paths, sha256_file

POWER_SIM = 10000
POWER_SEED = 20261001


def mcnemar_greater_sf(n_favor: np.ndarray, n_against: np.ndarray) -> np.ndarray:
    """One-sided exact McNemar p for the alternative that the favor count is larger."""
    n_disc = n_favor + n_against
    p = np.ones(len(n_favor), dtype=np.float64)
    ok = n_disc > 0
    # P(X >= n_favor) under Binomial(n_disc, 0.5). n_favor == 0 gives 1.
    p[ok] = binom.sf(n_favor[ok] - 1, n_disc[ok], 0.5)
    return p


def power_one(n: int, p_favor: float, p_against: float, rng: np.random.Generator) -> dict:
    draw = rng.random((POWER_SIM, n))
    favor = draw < p_favor
    against = (draw >= p_favor) & (draw < p_favor + p_against)
    p = mcnemar_greater_sf(favor.sum(axis=1), against.sum(axis=1))
    return {
        "n": int(n),
        "p_favor": p_favor,
        "p_against": p_against,
        "n_sim": POWER_SIM,
        "power": float(np.mean(p < 0.05)),
        "n_significant": int(np.sum(p < 0.05)),
    }


def tumour_content_edges(values: pd.Series) -> dict:
    raw = pd.Series(pd.qcut(values.to_numpy(), 3, duplicates="drop"))
    intervals = list(raw.cat.categories)
    labels = [f"T{i + 1}" for i in range(len(intervals))]
    assigned = raw.map({interval: labels[i] for i, interval in enumerate(intervals)})
    rows = []
    for interval, label in zip(intervals, labels):
        rows.append({
            "tertile": label,
            "left": float(interval.left),
            "right": float(interval.right),
            "closed": interval.closed,
            "n": int((assigned == label).sum()),
        })
    return {"n_bins": len(intervals), "bins": rows, "n": int(values.shape[0])}


def write_references(paths, out_dir: Path) -> None:
    import joblib

    from stage2_common import FIXED_SITE_TISSUES, load_gene_pack, load_tpm, sum_1e6, train_frame
    from stage2_fit import ld_references
    from importlib import import_module
    ev = import_module("09_eval")

    choice = json.loads((paths["results"] / "stage2" / "proxy_choice.json").read_text())
    site_tissues = dict(FIXED_SITE_TISSUES)
    site_tissues["lymph_node"] = [choice["ln_proxy"]]
    site_tissues["bone_marrow"] = [choice["bm_proxy"]]
    genes, _b0, _off, _ind, set_names = load_gene_pack(paths)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    train = train_frame(paths)
    tumor_norm = sum_1e6(tpm[[row[sample] for sample in train["sample"]]])
    mu, _h_dummy, label_names = ld_references(
        tumor_norm, train["learning_label"].to_numpy(), tumor_norm[:1]
    )
    del _h_dummy
    h_by_site = ev.build_h(paths, len(genes), site_tissues, samples, tpm, row)
    v0 = ev.v0_site_means(paths, set_names)
    keep = np.array([not gene.startswith("MT-") for gene in genes])
    payload = {
        "genes": np.array(genes, dtype=object),
        "mu": mu.astype(np.float32),
        "mu_labels": np.array(label_names, dtype=object),
        "keep_mt": keep,
    }
    for site, vector in h_by_site.items():
        payload[f"h__{site}"] = vector.astype(np.float32)
    for site, vector in v0.items():
        payload[f"v0__{site}"] = vector.astype(np.float64)
    np.savez_compressed(out_dir / "references.npz", **payload)
    meta = {
        "mu_labels": label_names,
        "h_sites": sorted(h_by_site),
        "v0_sites": sorted(v0),
        "n_genes": len(genes),
        "n_mu": int(mu.shape[0]),
        "site_tissues": site_tissues,
        "ln_proxy": choice["ln_proxy"],
        "bm_proxy": choice["bm_proxy"],
    }
    (out_dir / "references_meta.json").write_text(json.dumps(meta, indent=2))
    del tpm, tumor_norm, mu


def model_paths(paths) -> list[tuple[str, str]]:
    from stage2_common import SITE_MODEL_ORDER

    rows = []
    model_dir = paths["results"] / "stage2" / "models"
    names = ["BASE_Z", "BASE_K", "SA_Z", "SA_K", "NC_Z"]
    names += [f"SC_{site}_Z" for site in SITE_MODEL_ORDER]
    for name in names:
        rows.append((f"results/stage2/models/{name}.joblib", "model"))
    for path in sorted((paths["config"] / "mappings").glob("*.tsv")):
        rows.append((str(path.relative_to(paths["root"])), "mapping"))
    rows.append(("results/B0_genes.txt", "b0_genes"))
    rows.append(("data/processed/genes/G_symbols.txt", "G"))
    rows.append(("data/processed/genes/gene_sets.npz", "ks_sets"))
    rows.append(("results/stage3/references.npz", "ld_v0_references"))
    return rows


def main() -> None:
    paths = load_paths()
    out = paths["results"] / "stage3"
    out.mkdir(parents=True, exist_ok=True)
    s1 = read_pog570_table("s1")
    labels = label_table(s1)
    if int(labels["PATIENT_ID"].duplicated().sum()) != 0:
        raise SystemExit("duplicate PATIENT_ID")
    # NSCLC codes stay out of rule 1. True small-cell codes stay in rule 1.
    ht = labels["HISTOLOGICAL_TYPE"].astype(str).str.strip().str.upper()
    if not (labels.loc[ht.isin(["NSCLC", "NSCLCPD"]), "rule"] != "1").all():
        raise SystemExit("NSCLC was assigned to rule 1")
    if not (labels.loc[ht == "SCLC", "rule"] == "1").all():
        raise SystemExit("SCLC was not assigned to rule 1")
    label_path = paths["config"] / "pog570_eval_labels.tsv"
    labels.to_csv(label_path, sep="\t", index=False)
    digest = sha256_file(label_path)
    (paths["config"] / "pog570_eval_labels.sha256").write_text(digest + "\n")

    site = attach_site(s1)
    merged = labels.merge(site, on="PATIENT_ID", how="left")
    flags = [set_flags(organ, std) for organ, std in zip(merged["organ"], merged["standard_site"])]
    for key in ("in_eval", "in_native_truth", "at_risk"):
        merged[key] = [row[key] for row in flags]
    organ_counts = merged["organ"].value_counts(dropna=False).rename_axis("organ").reset_index(name="n")
    organ_counts.to_csv(out / "organ_counts.tsv", sep="\t", index=False)
    excluded = merged.loc[merged["organ"] == "exclude"]
    reason = excluded["rule"].value_counts().rename_axis("rule").reset_index(name="n")
    reason.to_csv(out / "exclude_reason_counts.tsv", sep="\t", index=False)
    rule_organ = merged.groupby(["rule", "organ"], dropna=False).size().reset_index(name="n")
    rule_organ.to_csv(out / "rule_organ_counts.tsv", sep="\t", index=False)

    site_rows = []
    for std, sub in merged.groupby("standard_site"):
        site_rows.append({
            "standard_site": std,
            "n": int(len(sub)),
            "n_eval": int(sub["in_eval"].sum()),
            "n_at_risk": int(sub["at_risk"].sum()),
            "n_native_truth": int(sub["in_native_truth"].sum()),
        })
    site_df = pd.DataFrame(site_rows).sort_values("standard_site")
    site_df.to_csv(out / "set_sizes_by_site.tsv", sep="\t", index=False)
    met_rows = []
    for level, sub in merged.groupby("METASTATIC_OR_RECURRENCE", dropna=False):
        met_rows.append({
            "METASTATIC_OR_RECURRENCE": level,
            "n": int(len(sub)),
            "n_eval": int(sub["in_eval"].sum()),
            "n_at_risk": int(sub["at_risk"].sum()),
            "n_native_truth": int(sub["in_native_truth"].sum()),
        })
    met_df = pd.DataFrame(met_rows)
    met_df.to_csv(out / "set_sizes_by_metastatic.tsv", sep="\t", index=False)

    n_eval = int(merged["in_eval"].sum())
    n_risk = int(merged["at_risk"].sum())
    n_native = int(merged["in_native_truth"].sum())
    tc = merged.loc[merged["in_eval"], "TUMOUR_CONTENT_num"]
    if int(tc.isna().sum()) != 0:
        raise SystemExit("TUMOUR_CONTENT missing inside the eval set")
    edges = tumour_content_edges(tc)
    (out / "tc_tertile_edges.json").write_text(json.dumps(edges, indent=2))

    # One Generator, scenario order H1 full, H1 half, H2 full, H2 half.
    rng = np.random.default_rng(POWER_SEED)
    h1_favor, h1_against = 34 / 361, 1 / 361
    h2_favor, h2_against = 29 / 437, 12 / 437
    power = {
        "seed": POWER_SEED,
        "n_sim": POWER_SIM,
        "alpha": 0.05,
        "alternative": "one-sided exact McNemar, favor count greater than 0.5",
        "half_effect": "p_favor is halved and p_against is unchanged",
        "H1_full": power_one(n_risk, h1_favor, h1_against, rng),
        "H1_half": power_one(n_risk, h1_favor / 2, h1_against, rng),
        "H2_full": power_one(n_eval, h2_favor, h2_against, rng),
        "H2_half": power_one(n_eval, h2_favor / 2, h2_against, rng),
    }
    check = binomtest(8, 10, 0.5, alternative="greater").pvalue
    got = float(mcnemar_greater_sf(np.array([8]), np.array([2]))[0])
    if abs(check - got) > 1e-12:
        raise SystemExit(f"mcnemar sf mismatch {check} {got}")
    (out / "power.json").write_text(json.dumps(power, indent=2))

    summary = {
        "n_samples": int(len(merged)),
        "n_eval": n_eval,
        "n_at_risk": n_risk,
        "n_native_truth": n_native,
        "n_exclude": int((merged["organ"] == "exclude").sum()),
        "n_NA": int((merged["organ"] == "NA").sum()),
        "label_sha256": digest,
        "beta_cut": BETA_CUT,
        "fake_organ_weights": FAKE_ORGAN_WEIGHTS,
        "cohort_marginal_sum": int(sum(COHORT_MARGINAL_N.values())),
        "fake_weight_sum": int(sum(FAKE_ORGAN_WEIGHTS.values())),
    }
    (out / "set_summary.json").write_text(json.dumps(summary, indent=2))
    print("sets", summary, flush=True)
    print("tc", edges, flush=True)
    print("power", {key: power[key]["power"] for key in ("H1_full", "H1_half", "H2_full", "H2_half")}, flush=True)

    print("references", flush=True)
    write_references(paths, out)
    hash_rows = []
    for rel, role in model_paths(paths):
        path = paths["root"] / rel
        hash_rows.append({
            "path": rel,
            "role": role,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    pd.DataFrame(hash_rows).to_csv(paths["config"] / "frozen_models_A3.tsv", sep="\t", index=False)
    print("hashes", len(hash_rows), flush=True)


if __name__ == "__main__":
    main()
