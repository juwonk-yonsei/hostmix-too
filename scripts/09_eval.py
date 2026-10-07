#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

"""Section 5-7 evaluation plus the POG570 metadata tables. Requires analysis_plan_A2.yaml."""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr
from sklearn.model_selection import StratifiedKFold

from loader import read_pog570_table
from mapping_rules import MET500_COHORT_TO_ORGAN, PROJECT_TO_ORGAN, SITE_NATIVE, map_site
from metrics import ORGANS, decompose, macro_f1, mask_native, percentile_ci, topk_hit
from stage2_common import (
    FIXED_SITE_TISSUES,
    SIM_RHOS,
    SITE_MODEL_ORDER,
    SITE_TO_MODEL,
    gtex_ref_ids,
    gtex_sim_ids,
    load_gene_pack,
    load_tpm,
    rank_z,
    score_ks,
    simulation_tumor_index,
    site_native_organs,
    sum_1e6,
    test_frame,
    train_frame,
)
from stage2_fit import ld_apply, ld_references, predict_bundle
from util import SEED, load_paths

Z_METHODS = ["BASE", "M1", "NC", "LD", "SA", "SC", "SC+IF20", "SC+IF40", "IF20"]
K_METHODS = ["BASE", "M1", "V0", "NC", "LD", "SA", "SC", "SC+IF20", "SC+IF40", "IF20"]
CANDIDATES = ["SA", "SC", "SC+IF20", "SC+IF40", "IF20"]
SITE_LEVELS = ["liver", "lymph_node", "soft_tissue", "bone_marrow", "lung", "skin", "brain", "other"]
POG_COHORTS = ["LUNG", "PANC", "CNS-PNS", "SECR", "MISC", "BCC", "LYMP", "SARC", "HNSC", "OV"]
POG_WORDS = [
    "neuroendocrine", "small cell", "carcinoid", "net", "salivary",
    "adenoid cystic", "mucoepidermoid", "basal cell", "lymphoma", "glioma", "neuroblastoma",
]


def load_named(model_dir: Path, name: str) -> dict:
    return joblib.load(model_dir / f"{name}.joblib")


def take_X(bundle: dict, z_full: np.ndarray, k_full: np.ndarray, b0_idx: np.ndarray) -> np.ndarray:
    if bundle["representation"] == "Z":
        idx = np.asarray(bundle["feature_index"], dtype=np.int32) if "feature_index" in bundle else b0_idx
        return z_full[:, idx]
    if "feature_index" in bundle:
        return k_full[:, np.asarray(bundle["feature_index"], dtype=np.int32)]
    return k_full


def predict_named(bundle: dict, z_full, k_full, b0_idx):
    organ, pred, raw, _, _ = predict_bundle(bundle, take_X(bundle, z_full, k_full, b0_idx))
    return organ, pred, raw


def v0_site_means(paths, set_names) -> dict[str, np.ndarray]:
    gtex = pd.read_csv(paths["config"] / "split_gtex.tsv", sep="\t", dtype=str)
    ref = gtex.loc[gtex["split"] == "GTEx-ref"]
    ks = pd.read_parquet(paths["data_processed"] / "features/ks_gtex_ref.parquet")
    tissue_mean = {}
    for tissue, sub in ref.groupby("tissue"):
        ids = [sample for sample in sub["sample"] if sample in ks.index]
        if ids:
            tissue_mean[tissue] = ks.loc[ids, set_names].to_numpy(dtype=np.float64).mean(axis=0)
    site_mean = {}
    for site, (tissues, _, _) in SITE_NATIVE.items():
        present = [tissue_mean[tissue] for tissue in tissues if tissue in tissue_mean]
        if present:
            site_mean[site] = np.mean(np.vstack(present), axis=0)
    return site_mean


def build_h(paths, genes_n, site_tissues, samples, tpm, row) -> dict[str, np.ndarray]:
    needed = sorted({tissue for tissues in site_tissues.values() for tissue in tissues})
    ids = gtex_ref_ids(paths, needed)
    out = {}
    for site, tissues in site_tissues.items():
        blocks = []
        for tissue in tissues:
            blocks.append(sum_1e6(tpm[[row[sample] for sample in ids[tissue]]]))
        pooled = np.vstack(blocks)
        out[site] = sum_1e6(pooled.mean(axis=0, keepdims=True))[0]
    return out


def predict_all(tpm, sites, bundles, b0_idx, offsets, indices, mu, h_by_site, keep, v0_means):
    z_full = rank_z(tpm).astype(np.float32)
    k_full = score_ks(tpm, offsets, indices)
    k_v0 = k_full.astype(np.float64).copy()
    for i, site in enumerate(sites):
        host = v0_means.get(site)
        if host is not None:
            k_v0[i] = (k_full[i] - 0.3 * host) / 0.7
    k_v0 = k_v0.astype(np.float32)
    cleaned = np.array(tpm, dtype=np.float32, copy=True)
    beta = np.full(tpm.shape[0], np.nan, dtype=np.float64)
    for site in sorted(set(sites)):
        model_site = SITE_TO_MODEL.get(site)
        idx = np.array([i for i, value in enumerate(sites) if value == site])
        if model_site is None:
            continue
        part, part_beta = ld_apply(tpm[idx], h_by_site[model_site], mu, keep)
        cleaned[idx] = part
        beta[idx] = part_beta
    z_ld = rank_z(cleaned).astype(np.float32)
    k_ld = score_ks(cleaned, offsets, indices)
    natives = [site_native_organs(site) for site in sites]
    out = {}
    raw_nc = {}
    for rep, z_use, k_use, methods in (
        ("Z", z_full, k_full, Z_METHODS),
        ("K", z_full, k_full, K_METHODS),
    ):
        base_organ, base_pred, _ = predict_named(bundles[f"BASE_{rep}"], z_use, k_use, b0_idx)
        out[( "BASE", rep)] = (base_organ, base_pred)
        masked, m1_pred = mask_native(base_organ, [sorted(native) for native in natives])
        out[("M1", rep)] = (masked, m1_pred)
        if rep == "K":
            organ, pred, _ = predict_named(bundles["BASE_K"], z_full, k_v0, b0_idx)
            out[("V0", "K")] = (organ, pred)
        nc_organ, nc_pred, nc_raw = predict_named(bundles[f"NC_{rep}"], z_use, k_use, b0_idx)
        out[("NC", rep)] = (nc_organ, nc_pred)
        raw_nc[rep] = nc_raw
        ld_organ, ld_pred, _ = predict_named(bundles[f"BASE_{rep}"], z_ld, k_ld, b0_idx)
        out[("LD", rep)] = (ld_organ, ld_pred)
        sa_organ, sa_pred, _ = predict_named(bundles[f"SA_{rep}"], z_use, k_use, b0_idx)
        out[("SA", rep)] = (sa_organ, sa_pred)
        if_organ, if_pred, _ = predict_named(bundles[f"IF20_{rep}"], z_use, k_use, b0_idx)
        out[("IF20", rep)] = (if_organ, if_pred)
        for tag, prefix in (("SC", "SC"), ("SC+IF20", "SCIF20"), ("SC+IF40", "SCIF40")):
            organ = np.zeros((len(sites), len(ORGANS)), dtype=np.float64)
            pred = np.empty(len(sites), dtype=object)
            for site in sorted(set(sites)):
                idx = np.array([i for i, value in enumerate(sites) if value == site])
                model_site = SITE_TO_MODEL.get(site)
                name = f"SA_{rep}" if model_site is None else f"{prefix}_{model_site}_{rep}"
                part_organ, part_pred, _ = predict_named(bundles[name], z_use[idx], k_use[idx], b0_idx)
                organ[idx] = part_organ
                pred[idx] = part_pred
            out[(tag, rep)] = (organ, pred)
    return out, beta, natives, raw_nc


def truth_organs(labels: np.ndarray) -> np.ndarray:
    return np.array([PROJECT_TO_ORGAN[label] for label in labels], dtype=object)


def host_flags(truth, pred, natives) -> np.ndarray:
    out = np.zeros(len(truth), dtype=bool)
    for i, native in enumerate(natives):
        if native and truth[i] not in native and pred[i] in native:
            out[i] = True
    return out


def bootstrap_table(truth, pred, proba, natives, draws) -> dict:
    point = decompose(truth, pred, natives)
    point["top1"] = float(np.mean(pred == truth))
    point["top3"] = float(np.mean(topk_hit(proba, ORGANS, truth, 3, pred != "NA")))
    point["macro_f1"] = macro_f1(truth, pred)
    native_truth = np.array([bool(natives[i]) and truth[i] in natives[i] for i in range(len(truth))])
    point["n_native_truth"] = int(native_truth.sum())
    point["native_truth_top1"] = float(np.mean(pred[native_truth] == truth[native_truth])) if native_truth.any() else None
    store = {key: [] for key in ("top1", "top3", "macro_f1", "H", "host_rate")}
    undefined = 0
    for draw in draws:
        sub_t, sub_p = truth[draw], pred[draw]
        store["top1"].append(float(np.mean(sub_p == sub_t)))
        store["top3"].append(float(np.mean(topk_hit(proba[draw], ORGANS, sub_t, 3, sub_p != "NA"))))
        store["macro_f1"].append(macro_f1(sub_t, sub_p))
        stats = decompose(sub_t, sub_p, [natives[i] for i in draw])
        if stats["H"] is None:
            undefined += 1
        else:
            store["H"].append(stats["H"])
        if stats["host_rate"] is not None:
            store["host_rate"].append(stats["host_rate"])
    point["top1_ci"] = percentile_ci(store["top1"])
    point["top3_ci"] = percentile_ci(store["top3"])
    point["macro_f1_ci"] = percentile_ci(store["macro_f1"])
    point["H_ci"] = percentile_ci(store["H"])
    point["host_rate_ci"] = percentile_ci(store["host_rate"])
    point["n_H_undefined"] = undefined
    return point


def mcnemar(base_ok, new_ok):
    n10 = int(np.sum(base_ok & ~new_ok))
    n01 = int(np.sum(~base_ok & new_ok))
    n_disc = n10 + n01
    pvalue = None if n_disc == 0 else float(binomtest(n01, n_disc, 0.5, alternative="two-sided").pvalue)
    return n10, n01, pvalue


def paired_diff(truth, base_pred, new_pred, draws):
    point = float(np.mean(new_pred == truth) - np.mean(base_pred == truth))
    diffs = [float(np.mean(new_pred[d] == truth[d]) - np.mean(base_pred[d] == truth[d])) for d in draws]
    return point, percentile_ci(diffs)


def conformal_threshold(scores: np.ndarray, alpha: float) -> float:
    n_cal = len(scores)
    k = math.ceil((n_cal + 1) * (1.0 - alpha))
    if k > n_cal:
        return math.inf
    ordered = np.sort(scores)
    return float(ordered[k - 1])


def conformal_sets(proba, threshold):
    scores = 1.0 - proba
    return scores <= threshold


def pog_tables(out: Path) -> None:
    s1 = read_pog570_table("s1")
    full = pd.crosstab(s1["ANALYSIS_COHORT"], s1["TUMOUR_TYPE"], dropna=False)
    full.to_csv(out / "pog_cohort_by_tumour_type.tsv", sep="\t")
    hist = pd.crosstab(s1["ANALYSIS_COHORT"], s1["HISTOLOGICAL_TYPE"], dropna=False)
    hist.to_csv(out / "pog_cohort_by_histology.tsv", sep="\t")
    rows = []
    for cohort in POG_COHORTS:
        sub = s1.loc[s1["ANALYSIS_COHORT"] == cohort]
        for column in ("TUMOUR_TYPE", "HISTOLOGICAL_TYPE"):
            counts = sub[column].fillna("<NA>").value_counts()
            for value, count in counts.items():
                rows.append({"cohort": cohort, "field": column, "value": value, "n": int(count)})
    pd.DataFrame(rows).to_csv(out / "pog_listed_cohort_types.tsv", sep="\t", index=False)
    word_rows = []
    text = (s1["TUMOUR_TYPE"].fillna("") + " " + s1["HISTOLOGICAL_TYPE"].fillna("")).str.casefold()
    for word in POG_WORDS:
        hit = text.str.contains(word, regex=False)
        for cohort, sub in s1.groupby("ANALYSIS_COHORT"):
            word_rows.append({
                "word": word,
                "ANALYSIS_COHORT": cohort,
                "n": int(hit.loc[sub.index].sum()),
                "n_cohort": int(len(sub)),
            })
    pd.DataFrame(word_rows).to_csv(out / "pog_keyword_counts.tsv", sep="\t", index=False)


def main() -> None:
    paths = load_paths()
    if not (paths["config"] / "analysis_plan_A2.yaml").exists():
        raise SystemExit("config/analysis_plan_A2.yaml is missing. Section 5 was not started.")
    out = paths["results"] / "stage2"
    choice = json.loads((out / "proxy_choice.json").read_text())
    ln_proxy, bm_proxy = choice["ln_proxy"], choice["bm_proxy"]
    site_tissues = dict(FIXED_SITE_TISSUES)
    site_tissues["lymph_node"] = [ln_proxy]
    site_tissues["bone_marrow"] = [bm_proxy]
    genes, b0_idx, offsets, indices, set_names = load_gene_pack(paths)
    model_dir = out / "models"
    bundles = {}
    for rep in ("Z", "K"):
        for name in ("BASE", "SA", "NC", "IF20"):
            bundles[f"{name}_{rep}"] = load_named(model_dir, f"{name}_{rep}")
        for site in SITE_MODEL_ORDER:
            for prefix in ("SC", "SCIF20", "SCIF40"):
                bundles[f"{prefix}_{site}_{rep}"] = load_named(model_dir, f"{prefix}_{site}_{rep}")
    print("references", flush=True)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    train = train_frame(paths)
    tumor_norm = sum_1e6(tpm[[row[sample] for sample in train["sample"]]])
    mu, _, label_names = ld_references(tumor_norm, train["learning_label"].to_numpy(), tumor_norm[:1])
    # h is rebuilt per site; the dummy host argument above is replaced.
    del _
    h_by_site = build_h(paths, len(genes), site_tissues, samples, tpm, row)
    keep = np.array([not gene.startswith("MT-") for gene in genes])
    (out / "ld_label_order.json").write_text(json.dumps({"n_labels": len(label_names), "labels": label_names}))
    v0_means = v0_site_means(paths, set_names)

    print("MET500", flush=True)
    meta = pd.read_csv(paths["config"] / "met500_eval_samples.tsv", sep="\t", dtype=str)
    meta["organ"] = meta["cohort"].map(MET500_COHORT_TO_ORGAN)
    meta["standard_site"] = [map_site("met500_biopsy_tissue", value) for value in meta["biopsy_tissue"]]
    meta["tc_num"] = pd.to_numeric(meta["tc"], errors="coerce")
    eval_df = meta.loc[meta["organ"].notna() & (meta["organ"] != "unmapped")].reset_index(drop=True)
    met_samples, met_tpm = load_tpm(paths["data_processed"] / "met500/tpm_G.h5")
    met_row = {sample: i for i, sample in enumerate(met_samples)}
    block = sum_1e6(met_tpm[[met_row[sample] for sample in eval_df["Sample_id"]]])
    preds, beta, natives, raw_nc = predict_all(
        block, eval_df["standard_site"].tolist(), bundles, b0_idx, offsets, indices, mu, h_by_site, keep, v0_means
    )
    truth = eval_df["organ"].to_numpy(dtype=object)
    draws = np.random.default_rng(SEED).integers(0, len(truth), size=(2000, len(truth)))
    rows = []
    stored = {}
    for (method, rep), (proba, pred) in preds.items():
        stats = bootstrap_table(truth, pred, proba, natives, draws)
        stats.update({"method": method, "representation": rep})
        rows.append(stats)
        stored[(method, rep)] = (pred, proba, host_flags(truth, pred, natives))
        print("met", method, rep, stats["top1"], flush=True)
    overall = pd.DataFrame(rows)
    overall.to_csv(out / "met500_overall.tsv", sep="\t", index=False)
    nc_rows = []
    for rep, raw in raw_nc.items():
        normal = pd.Series(raw).astype(str)
        normal = normal[normal.str.startswith("N_")]
        for name, count in normal.value_counts().items():
            nc_rows.append({"representation": rep, "argmax_class": name, "n": int(count)})
        nc_rows.append({"representation": rep, "argmax_class": "<any normal>", "n": int(normal.shape[0])})
    pd.DataFrame(nc_rows).to_csv(out / "nc_normal_argmax.tsv", sep="\t", index=False)

    base_z_pred, _, base_z_host = stored[("BASE", "Z")]
    paired_rows = []
    for (method, rep), (pred, _, flags) in stored.items():
        if method == "BASE":
            continue
        for base_name, base_pred, base_host in (
            ("BASE-Z", base_z_pred, base_z_host),
            ("BASE-K", stored[("BASE", "K")][0], stored[("BASE", "K")][2]),
        ):
            if base_name == "BASE-K" and rep != "K":
                continue
            diff, ci = paired_diff(truth, base_pred, pred, draws)
            n10, n01, pvalue = mcnemar(base_pred == truth, pred == truth)
            paired_rows.append({
                "method": method,
                "representation": rep,
                "baseline": base_name,
                "top1_diff": diff,
                "top1_diff_ci_low": ci[0],
                "top1_diff_ci_high": ci[1],
                "n_base_only": n10,
                "n_method_only": n01,
                "mcnemar_p": pvalue,
                "n_host_fixed": int(np.sum(base_host & ~flags)),
                "n_host_new": int(np.sum(~base_host & flags)),
            })
    pd.DataFrame(paired_rows).to_csv(out / "met500_paired.tsv", sep="\t", index=False)

    site_rows = []
    site_group = np.array([site if site in SITE_LEVELS[:-1] else "other" for site in eval_df["standard_site"]])
    for method in ("BASE", "LD", "NC", "SA", "SC", "SC+IF20"):
        for rep in ("Z", "K"):
            pred = stored[(method, rep)][0]
            for level in SITE_LEVELS:
                idx = np.where(site_group == level)[0]
                stats = decompose(truth[idx], pred[idx], [natives[i] for i in idx])
                site_rows.append({
                    "method": method, "representation": rep, "site": level,
                    "n": stats["n"], "top1": float(np.mean(pred[idx] == truth[idx])) if len(idx) else None,
                    "n_error": stats["n_error"], "n_host": stats["n_host"], "host_rate": stats["host_rate"],
                })
    pd.DataFrame(site_rows).to_csv(out / "met500_by_site.tsv", sep="\t", index=False)

    core = ~((eval_df["cohort"] == "HNSC") & (eval_df["tissue"] == "parotid"))
    core_idx = np.where(core.to_numpy())[0]
    core_rows = []
    for (method, rep), (pred, proba, _) in stored.items():
        stats = decompose(truth[core_idx], pred[core_idx], [natives[i] for i in core_idx])
        core_rows.append({
            "method": method, "representation": rep, "n": int(len(core_idx)),
            "n_dropped_parotid": int((~core).sum()),
            "top1": float(np.mean(pred[core_idx] == truth[core_idx])),
            "host_rate": stats["host_rate"], "n_at_risk": stats["n_at_risk"], "n_host": stats["n_host"],
        })
    pd.DataFrame(core_rows).to_csv(out / "met500_core.tsv", sep="\t", index=False)

    beta_df = eval_df[["Sample_id", "standard_site", "tc", "organ"]].copy()
    beta_df["beta"] = beta
    beta_df.to_csv(out / "ld_beta_samples.tsv", sep="\t", index=False)
    beta_rows = []
    finite = np.isfinite(beta)
    for site, idx in beta_df.groupby("standard_site").groups.items():
        values = beta[np.asarray(list(idx))]
        values = values[np.isfinite(values)]
        beta_rows.append({
            "site": site, "n": int(values.size),
            "q25": float(np.quantile(values, 0.25)) if values.size else None,
            "median": float(np.median(values)) if values.size else None,
            "q75": float(np.quantile(values, 0.75)) if values.size else None,
        })
    pd.DataFrame(beta_rows).to_csv(out / "ld_beta_by_site.tsv", sep="\t", index=False)
    tc = eval_df["tc_num"].to_numpy()
    corr_rows = []
    mask = finite & np.isfinite(tc)
    if mask.sum() >= 3:
        result = spearmanr(beta[mask], tc[mask])
        corr_rows.append({"site": "all", "n": int(mask.sum()), "spearman": float(result.statistic), "p": float(result.pvalue)})
    for site in sorted(set(eval_df["standard_site"])):
        part = mask & (eval_df["standard_site"].to_numpy() == site)
        if part.sum() < 3:
            corr_rows.append({"site": site, "n": int(part.sum()), "spearman": None, "p": None})
            continue
        result = spearmanr(beta[part], tc[part])
        corr_rows.append({"site": site, "n": int(part.sum()), "spearman": float(result.statistic), "p": float(result.pvalue)})
    pd.DataFrame(corr_rows).to_csv(out / "ld_beta_spearman.tsv", sep="\t", index=False)
    tertile = pd.Series(np.nan, index=np.arange(len(beta)), dtype=object)
    # pd.qcut(..., 3, labels=["T1","T2","T3"], duplicates="drop") cannot split
    # the mass of exact zeros. Equal-count tertiles use a stable argsort of
    # the finite values and np.array_split into 3. Tied values that cross a
    # cut stay in evaluation-table order.
    if finite.sum() >= 3:
        finite_idx = np.where(finite)[0]
        order = np.argsort(beta[finite_idx], kind="mergesort")
        for i, part in enumerate(np.array_split(order, 3)):
            tertile.iloc[finite_idx[part]] = f"T{i + 1}"
        print("beta rank tertile sizes", [int((tertile == name).sum()) for name in ("T1", "T2", "T3")], flush=True)
    tertile_rows = []
    base_pred = stored[("BASE", "Z")][0]
    for level, idx in tertile.groupby(tertile).groups.items():
        if level != level:
            continue
        idx = np.asarray(list(idx))
        stats = decompose(truth[idx], base_pred[idx], [natives[i] for i in idx])
        tertile_rows.append({
            "tertile": level, "n": int(len(idx)),
            "beta_min": float(np.min(beta[idx])), "beta_max": float(np.max(beta[idx])),
            "BASE_Z_top1": float(np.mean(base_pred[idx] == truth[idx])),
            "host_rate": stats["host_rate"], "n_host": stats["n_host"], "n_at_risk": stats["n_at_risk"],
        })
    pd.DataFrame(tertile_rows).to_csv(out / "ld_beta_tertile.tsv", sep="\t", index=False)

    rank_rows = []
    for method in CANDIDATES:
        for rep in ("Z", "K"):
            stats = overall.loc[(overall["method"] == method) & (overall["representation"] == rep)].iloc[0]
            rank_rows.append({
                "method_id": f"{method}-{rep}",
                "method": method,
                "representation": rep,
                "top1": float(stats["top1"]),
                "host_rate": float(stats["host_rate"]) if pd.notna(stats["host_rate"]) else 1.0,
            })
    rank = pd.DataFrame(rank_rows).sort_values(["top1", "host_rate", "method_id"], ascending=[False, True, True])
    rank.to_csv(out / "candidate_rank.tsv", sep="\t", index=False)
    top2 = rank.head(2)["method_id"].tolist()
    (out / "top2.json").write_text(json.dumps({"top2": top2}, indent=2))
    shift_rows = []
    cell_rows = []
    for method_id in top2:
        method, rep = method_id.rsplit("-", 1)
        pred = stored[(method, rep)][0]
        both = pd.DataFrame({"truth": truth, "base": base_z_pred, "new": pred})
        for base_ok in (True, False):
            for new_ok in (True, False):
                n = int(np.sum(((base_z_pred == truth) == base_ok) & ((pred == truth) == new_ok)))
                shift_rows.append({"method_id": method_id, "base_correct": base_ok, "new_correct": new_ok, "n": n})
        changed = both.loc[both["base"] != both["new"]]
        counts = changed.value_counts().reset_index(name="n").head(10)
        counts.insert(0, "method_id", method_id)
        cell_rows.append(counts)
    pd.DataFrame(shift_rows).to_csv(out / "error_shift_2x2.tsv", sep="\t", index=False)
    pd.concat(cell_rows, ignore_index=True).to_csv(out / "error_shift_cells.tsv", sep="\t", index=False)

    print("conformal", top2, flush=True)
    conf_methods = [("BASE", "Z"), ("BASE", "K")]
    for method_id in top2:
        method, rep = method_id.rsplit("-", 1)
        conf_methods.append((method, rep))
    group = np.array([
        "liver" if site == "liver" else "lymph_node" if site == "lymph_node" else "other"
        for site in eval_df["standard_site"]
    ])
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    conf_rows = []
    for method, rep in conf_methods:
        proba = stored[(method, rep)][1]
        for variant in ("global", "mondrian"):
            for alpha in (0.1, 0.2):
                covered = np.zeros(len(truth), dtype=bool)
                sizes = np.zeros(len(truth), dtype=int)
                single = np.zeros(len(truth), dtype=object)
                for train_idx, test_idx in splitter.split(np.zeros(len(truth)), group):
                    if variant == "global":
                        scores = 1.0 - proba[train_idx, [ORGANS.index(label) for label in truth[train_idx]]]
                        threshold = conformal_threshold(scores, alpha)
                        selected = conformal_sets(proba[test_idx], threshold)
                    else:
                        selected = np.zeros((len(test_idx), len(ORGANS)), dtype=bool)
                        for level in ("liver", "lymph_node", "other"):
                            cal = train_idx[group[train_idx] == level]
                            te = np.where(group[test_idx] == level)[0]
                            if len(te) == 0:
                                continue
                            if len(cal) == 0:
                                selected[te] = True
                                continue
                            scores = 1.0 - proba[cal, [ORGANS.index(label) for label in truth[cal]]]
                            selected[te] = conformal_sets(proba[test_idx[te]], conformal_threshold(scores, alpha))
                    local_sizes = selected.sum(axis=1)
                    sizes[test_idx] = local_sizes
                    for row_i, global_i in enumerate(test_idx):
                        covered[global_i] = bool(selected[row_i, ORGANS.index(truth[global_i])])
                        if local_sizes[row_i] == 1:
                            single[global_i] = ORGANS[int(np.flatnonzero(selected[row_i])[0])]
                def emit(level, mask):
                    n = int(mask.sum())
                    sing = (sizes[mask] == 1) if n else np.array([], dtype=bool)
                    acc = None
                    if n and sing.any():
                        acc = float(np.mean(single[mask][sing] == truth[mask][sing]))
                    conf_rows.append({
                        "method": method, "representation": rep, "variant": variant, "alpha": alpha,
                        "level": level, "n": n,
                        "coverage": float(np.mean(covered[mask])) if n else None,
                        "mean_set_size": float(np.mean(sizes[mask])) if n else None,
                        "singleton_rate": float(np.mean(sizes[mask] == 1)) if n else None,
                        "singleton_accuracy": acc,
                    })
                emit("all", np.ones(len(truth), dtype=bool))
                for level in ("liver", "lymph_node", "other"):
                    emit(level, group == level)
    pd.DataFrame(conf_rows).to_csv(out / "conformal.tsv", sep="\t", index=False)

    print("simulation", flush=True)
    test = test_frame(paths)
    take = simulation_tumor_index(test)
    tumor_ids = test["sample"].to_numpy()[take]
    tumor_organ = truth_organs(test["learning_label"].to_numpy()[take])
    tumor_block = sum_1e6(tpm[[row[sample] for sample in tumor_ids]])
    sim_tissues = ["Liver", "Lung", "Adipose - Subcutaneous", ln_proxy, bm_proxy]
    sim_site = {
        "Liver": "liver", "Lung": "lung", "Adipose - Subcutaneous": "soft_tissue",
        ln_proxy: "lymph_node", bm_proxy: "bone_marrow",
    }
    sim_ids = gtex_sim_ids(paths, sim_tissues)
    host_bank = {}
    for tissue, ids in sim_ids.items():
        mat = sum_1e6(tpm[[row[sample] for sample in ids]])
        for i, sample in enumerate(ids):
            host_bank[sample] = mat[i]
    del tpm
    rng = np.random.default_rng(SEED)
    draws_host = {}
    for tissue in sim_tissues:
        ids = sim_ids[tissue]
        for rho in SIM_RHOS:
            picked = rng.integers(0, len(ids), size=len(tumor_ids))
            draws_host[(tissue, rho)] = [ids[int(i)] for i in picked]
    sim_rows = []
    for tissue in sim_tissues:
        site = sim_site[tissue]
        native = site_native_organs(site)
        for rho in SIM_RHOS:
            hosts = np.stack([host_bank[sample] for sample in draws_host[(tissue, rho)]])
            rho_v = np.float32(rho)
            mixed = rho_v * tumor_block + (np.float32(1.0) - rho_v) * hosts
            print("sim", tissue, rho, flush=True)
            pred_map, _, _, _ = predict_all(
                mixed, [site] * len(tumor_ids), bundles, b0_idx, offsets, indices, mu, h_by_site, keep, v0_means
            )
            at_risk = np.array([bool(native) and organ not in native for organ in tumor_organ])
            native_truth = np.array([bool(native) and organ in native for organ in tumor_organ])
            for (method, rep), (_, pred) in pred_map.items():
                pulled = at_risk & np.array([pred[i] in native for i in range(len(pred))])
                sim_rows.append({
                    "tissue": tissue, "site": site, "rho": rho, "method": method, "representation": rep,
                    "n": int(len(pred)),
                    "top1": float(np.mean(pred == tumor_organ)),
                    "n_at_risk": int(at_risk.sum()),
                    "host_pull_rate": float(pulled.sum() / at_risk.sum()) if at_risk.any() else None,
                    "n_native_truth": int(native_truth.sum()),
                    "native_truth_top1": float(np.mean(pred[native_truth] == tumor_organ[native_truth])) if native_truth.any() else None,
                })
    pd.DataFrame(sim_rows).to_csv(out / "simulation.tsv", sep="\t", index=False)
    print("POG metadata", flush=True)
    pog_tables(out)
    print("eval done", flush=True)


if __name__ == "__main__":
    main()
