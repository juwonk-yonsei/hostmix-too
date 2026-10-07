"""Expected contents of every file in manuscript/bmc/figure_source/ rebuilt from results (A17 §4.5).

The builders do not import the figure scripts. figure_source files were written with
float_format "%.10g", so numeric cells are compared after the same rounding; strings exactly.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

import numcheck_A17 as NC

ROOT = Path(__file__).resolve().parents[1]
FS = ROOT / "manuscript" / "bmc" / "figure_source"
DERIVED = "results/stage10/derived"

POOL = ["Liver", "Lung", "Brain - Cortex", "Adrenal Gland", "Skin - Not Sun Exposed (Suprapubic)",
        "Muscle - Skeletal", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Spleen", "Whole Blood"]
LABEL = {
    "BASE-Z": "Baseline", "SA-Z": "HostMix-TOO", "SC-Z": "Site-specific mixtures", "NC-Z": "Normal classes",
    "LD-Z": "Linear deconvolution", "M1-Z": "Native-organ masking", "IF20-Z": "Gene removal 20%",
    "SC+IF20-Z": "Site-specific + removal 20%", "SC+IF40-Z": "Site-specific + removal 40%", "BASE-K": "Gene sets",
    "SA-K": "Gene sets + mixtures", "V0-K": "Gene sets + host correction", "SA-G": "Gated model",
    "SA-pool22": "22-tissue model", "SA-pool3": "3-tissue model", "BASE-MLP": "Perceptron",
    "SA-MLP": "Perceptron + mixtures", "PURE-Zs": "Standardized, pure (C = 0.03)",
    "PURE-Zs-w": "Standardized, pure (C = 0.15)", "MIX-Z0": "Unstandardized, mixtures (C = 0.1)",
    "SCOPE": "SCOPE", "CUP-AI-Dx": "CUP-AI-Dx",
}
COHORT = {
    "blca_iatlas_imvigor210_2017": "IMvigor210", "brca_iatlas_anders_2022": "Anders", "mel_dfci_2019": "DFCI",
    "paad_iatlas_prince_2022": "PRINCE", "GSE50760": "GSE50760", "prad_su2c_2019": "SU2C/PCF",
    "GSE209998": "AURORA US", "GSE41258": "GSE41258", "GSE14018": "GSE14018", "GSE71729": "GSE71729",
    "GSE74685": "GSE74685", "prad_fhcrc": "FHCRC", "MET500": "MET500", "POG570": "POG570",
}
AUX = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019", "paad_iatlas_prince_2022",
       "GSE50760", "prad_su2c_2019"]
MICROARRAY = ["GSE41258", "GSE14018", "GSE71729", "GSE74685", "prad_fhcrc"]
METHODS4 = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
COHORTS3 = ["MET500", "POG570", "aux_rnaseq"]


def label(method: str) -> str:
    if method in LABEL:
        return LABEL[method]
    if method.startswith("SA-LOHO-"):
        return f"Leave-one-host-out ({method[len('SA-LOHO-'):]})"
    return method


READ: set[str] = set()


def read(path: str) -> pd.DataFrame:
    READ.add(path)
    return pd.read_csv(ROOT / path, sep="\t")


def one(frame: pd.DataFrame, **keys) -> pd.Series:
    mask = np.ones(len(frame), dtype=bool)
    for k, v in keys.items():
        mask &= frame[k].eq(v).to_numpy()
    hit = frame.loc[mask]
    if len(hit) != 1:
        raise KeyError(f"{keys}: {len(hit)} rows")
    return hit.iloc[0]


def fig2_means() -> pd.DataFrame:
    sim = read("results/stage5/sim_ext/full.tsv")
    outside = sorted(set(sim["tissue"]) - set(POOL))
    parts = []
    for method, tissues, where in (("BASE-Z", POOL, "in"), ("SA-Z", POOL, "in"), ("BASE-Z", outside, "out"),
                                   ("SA-Z", outside, "out"), ("SA-pool22", outside, "out")):
        part = sim.loc[sim["method"].eq(method) & sim["tissue"].isin(tissues)]
        by_tissue = part.groupby(["rho", "tissue"], as_index=False)["host_pull_rate"].mean()
        mean = by_tissue.groupby("rho", as_index=False).agg(host_pull_rate=("host_pull_rate", "mean"),
                                                            n_tissues=("tissue", "nunique"))
        parts.append(mean.assign(method=method, pool=where)[["rho", "host_pull_rate", "method", "n_tissues", "pool"]])
    return pd.concat(parts, ignore_index=True)


def fig2_loho() -> pd.DataFrame:
    sim = read("results/stage5/sim_ext/full.tsv")
    proxy = {"Spleen": "Spleen (lymph node proxy)", "Whole Blood": "Whole Blood (bone marrow proxy)"}
    rows = []
    for tissue in POOL:
        hits = {"baseline": one(sim, method="BASE-Z", tissue=tissue, rho=0.6),
                "leave_one_host_out": one(sim, method=f"SA-LOHO-{tissue}", tissue=tissue, rho=0.6),
                "hostmix": one(sim, method="SA-Z", tissue=tissue, rho=0.6)}
        n = {int(h["n_at_risk"]) for h in hits.values()}
        if len(n) != 1:
            raise KeyError(f"{tissue}: at-risk n differs")
        rows.append({"tissue": tissue, "display": proxy.get(tissue, tissue), "n_at_risk": n.pop(),
                     **{k: float(h["host_rate"]) for k, h in hits.items()}})
    return pd.DataFrame(rows).sort_values("baseline", ascending=False, kind="mergesort").reset_index(drop=True)


def fig3_differences() -> pd.DataFrame:
    pog = read("results/stage4/confirm/primary.tsv")
    aux = read("results/stage7/confirm/tables/hypothesis_primary.tsv")
    overall = read("results/stage4/confirm/secondary_overall.tsv")
    external = read("results/stage8/external/tables/metrics.tsv")
    sets = read("results/stage9/set_metrics.tsv")

    def arm(kind: str, method: str) -> float:
        if kind == "pog_host":
            return float(one(overall, method=method)["host_rate"])
        if kind == "pog_top1":
            return float(one(overall, method=method)["top1"])
        if kind == "pog_native":
            return float(one(overall, method=method)["native_truth_top1"])
        if kind == "aux_host":
            return float(one(sets, analysis_cohort="aux_rnaseq", set="risk", method=method)["host_rate"])
        if kind == "aux_top1":
            return float(one(external, subset="common_label", analysis_cohort="aux_rnaseq", method=method)["top1"])
        return float(one(sets, analysis_cohort="aux_rnaseq", set="native_truth", method=method)["top1"])

    rows = []
    for panel, h, cohort, kind in (("a", "H1", "POG570", "pog_host"), ("a", "AH1", "Auxiliary RNA-seq", "aux_host"),
                                   ("b", "H2", "POG570", "pog_top1"), ("b", "AH2", "Auxiliary RNA-seq", "aux_top1"),
                                   ("c", "H3", "POG570", "pog_native"), ("c", "AH3", "Auxiliary RNA-seq", "aux_native")):
        if h.startswith("A"):
            hit = one(aux, hypothesis=h)
            diff = float(hit["diff_sa_minus_comparator"])
        else:
            hit = one(pog, hypothesis=h)
            diff = float(hit["diff"])
        rows.append({"panel": panel, "hypothesis": h, "cohort": cohort, "n": int(hit["n"]), "diff": diff,
                     "ci_low": float(hit["ci_low"]), "ci_high": float(hit["ci_high"]),
                     "onesided_low": float(hit["onesided_low"]) if panel == "c" else np.nan,
                     "baseline": arm(kind, "BASE-Z"), "hostmix": arm(kind, "SA-Z")})
    return pd.DataFrame(rows)


def fig5_ablation() -> pd.DataFrame:
    rna = read("results/stage9/ablation/metrics.tsv")
    micro = read("results/stage10/ablation_microarray/metrics.tsv")
    restriction = read("results/stage10/ablation_microarray/restriction.tsv")
    names = {"MET500": "MET500", "POG570": "POG570", "aux_rnaseq": "Auxiliary RNA-seq", "microarray": "Microarray"}
    rows = []
    for method in ["BASE-Z", "PURE-Zs", "PURE-Zs-w", "MIX-Z0", "SA-Z"]:
        for cohort in COHORTS3:
            hit = one(rna, cohort=cohort, method=method)
            rows.append({"panel": "a", "cohort": cohort, "method": method, "metric": "host_rate",
                         "value": float(hit["host_rate"]), "n": int(hit["n_at_risk"]), "source": "stage9 ablation"})
            rows.append({"panel": "b", "cohort": cohort, "method": method, "metric": "top1",
                         "value": float(hit["top1"]), "n": int(hit["n"]), "source": "stage9 ablation"})
        for panel, set_name, metric in (("a", "at_risk", "host_rate"), ("b", "evaluation", "top1")):
            hit = one(micro, model=method, set=set_name, cohort="all", metric=metric)
            rows.append({"panel": panel, "cohort": "microarray", "method": method, "metric": metric,
                         "value": float(hit["value"]), "n": int(hit["n"]), "source": "stage10 ablation"})
        for cohort in ("MET500", "POG570"):
            part = restriction.loc[restriction["model"].eq(method) & restriction["cohort"].eq(cohort)]
            rows.append({"panel": "c", "cohort": cohort, "method": method, "metric": "random_removal_mean_drop",
                         "value": float(part["random_mean"].mean()), "n": len(part), "source": "stage10 restriction"})
    frame = pd.DataFrame(rows)
    frame["display"] = frame["cohort"].map(names)
    frame["label"] = frame["method"].map(label)
    return frame


def fig6_published() -> pd.DataFrame:
    raw = read("results/stage8/external/tables/metrics.tsv")
    part = raw.loc[raw["subset"].eq("common_label") & raw["method"].isin(METHODS4)].copy()
    part["display"] = part["method"].map(label)
    return part.reset_index(drop=True)


def microarray_rates() -> pd.DataFrame:
    return read(f"{DERIVED}/microarray_rates.tsv")


def fig4a_cohorts() -> pd.DataFrame:
    rows = []
    external = read("results/stage8/external/tables/metrics.tsv")
    for cohort, group in (("MET500", "Development"), ("POG570", "Confirmation 1")):
        for method in ("BASE-Z", "SA-Z"):
            hit = one(external, subset="common_label", analysis_cohort=cohort, method=method)
            rows.append({"group": group, "cohort": cohort, "method": method, "host_rate": float(hit["host_rate"]),
                         "n_at_risk": int(hit["n_at_risk"]), "source": "stage8 metrics"})
    sub = read("results/stage9/subcohort_metrics.tsv")
    for cohort in AUX:
        for method in ("BASE-Z", "SA-Z"):
            hit = one(sub, cohort=cohort, method=method)
            rows.append({"group": "Confirmation 2", "cohort": COHORT[cohort], "method": method,
                         "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                         "source": "subcohort_metrics"})
    metrics = read("results/stage7/confirm/tables/metrics.tsv")
    for method in ("BASE-Z", "SA-Z"):
        hit = one(metrics, cohort="GSE209998", standard_site="all", method=method)
        rows.append({"group": "Descriptive RNA-seq", "cohort": COHORT["GSE209998"], "method": method,
                     "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                     "source": "stage7 metrics"})
    micro = microarray_rates()
    for cohort in MICROARRAY:
        for method in ("BASE-Z", "SA-Z"):
            hit = one(micro, cohort=cohort, method=method)
            rows.append({"group": "Descriptive microarray", "cohort": COHORT[cohort], "method": method,
                         "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                         "source": "post hoc aggregate of confirmatory predictions"})
    frame = pd.DataFrame(rows)
    frame["group_display"] = frame["group"].map({"Development": "Development", "Confirmation 1": "Confirmation 1",
                                                 "Confirmation 2": "Confirmation 2",
                                                 "Descriptive RNA-seq": "Other RNA-seq",
                                                 "Descriptive microarray": "Microarray"})
    return frame


def fig4b(kind: str) -> pd.DataFrame:
    path = {"site": "results/stage4/confirm/secondary_by_site.tsv", "tc": "results/stage4/confirm/secondary_tc.tsv"}[kind]
    frame = read(path)
    return frame.loc[frame["method"].isin(["BASE-Z", "SA-Z"])].reset_index(drop=True)


def copy(path: str) -> pd.DataFrame:
    return read(path)


def fig4d_removed() -> pd.DataFrame:
    cases = read("results/stage8/diagnostics/error_shift_cases.tsv")
    counts = cases.groupby("analysis").agg(n=("SA_correct", "size"), correct=("SA_correct", "sum")).reset_index()
    counts["other"] = counts["n"] - counts["correct"]
    return counts


def s1_heatmap() -> pd.DataFrame:
    rows = []
    for _, r in read("results/stage2/met500_overall.tsv").iterrows():
        rows.append({"cohort": "MET500", "method": f"{r['method']}-{r['representation']}", "top1": float(r["top1"]),
                     "host_rate": float(r["host_rate"]), "status": "exploratory"})
    for _, r in read("results/stage4/confirm/secondary_overall.tsv").iterrows():
        rows.append({"cohort": "POG570", "method": r["method"], "top1": float(r["top1"]),
                     "host_rate": float(r["host_rate"]), "status": "pre-specified descriptive"})
    for cohort, path in (("MET500", "results/stage5/posthoc/met500_overall.tsv"),
                         ("POG570", "results/stage5/posthoc/pog_overall.tsv")):
        overall = read(path)
        for method in ("SA-pool22", "SA-MLP"):
            hit = one(overall, method=method)
            rows.append({"cohort": cohort, "method": method, "top1": float(hit["top1"]),
                         "host_rate": float(hit["host_rate"]), "status": "exploratory"})
    for cohort, path in (("MET500", "results/stage5/posthoc/met500_rules.tsv"),
                         ("POG570", "results/stage5/posthoc/pog_rules.tsv")):
        hit = one(read(path), rule="G-beta(0.02)")
        rows.append({"cohort": cohort, "method": "SA-G", "top1": float(hit["top1"]),
                     "host_rate": float(hit["host_rate"]), "status": "exploratory"})
    micro = microarray_rates()
    for method in ("BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z", "BASE-K", "SA-K",
                   "V0-K", "SA-MLP"):
        hit = one(micro, cohort="all", method=method)
        rows.append({"cohort": "Microarray", "method": method, "top1": float(hit["top1"]),
                     "host_rate": float(hit["host_rate"]), "status": "post hoc aggregate"})
    aux = read(f"{DERIVED}/aux_six_rates.tsv")
    for method in ("BASE-Z", "SA-Z", "SC-Z", "LD-Z", "NC-Z", "M1-Z", "SA-G", "SA-pool22", "BASE-K", "SA-K",
                   "V0-K", "SA-MLP"):
        hit = one(aux, method=method)
        rows.append({"cohort": "Auxiliary RNA-seq", "method": method, "top1": float(hit["top1"]),
                     "host_rate": float(hit["host_rate"]), "status": "post hoc aggregate"})
    for row in rows:
        if row["cohort"] == "POG570" and row["method"] in ("BASE-Z", "SA-Z") \
                and row["status"] == "pre-specified descriptive":
            row["status"] = NC.status_of("results/stage4/confirm/secondary_overall.tsv", f"method == '{row['method']}'")
    return pd.DataFrame(rows)


def s1_n() -> pd.DataFrame:
    """Set sizes printed under the Figure S1 columns (A18)."""
    rows = []
    for cohort, path, key, n_col in (
        ("MET500", "results/stage2/met500_overall.tsv", {}, "n"),
        ("POG570", "results/stage4/confirm/secondary_overall.tsv", {}, "n"),
        ("Auxiliary RNA-seq", f"{DERIVED}/aux_six_rates.tsv", {}, "n_evaluation"),
        ("Microarray", f"{DERIVED}/microarray_rates.tsv", {"cohort": "all"}, "n_evaluation"),
    ):
        frame = read(path)
        for column, value in key.items():
            frame = frame.loc[frame[column].eq(value)]
        for panel, column in (("a", n_col), ("b", "n_at_risk")):
            values = set(frame[column].astype(int))
            if len(values) != 1:
                raise NC.SourceError(f"{path}: {column} differs between rows")
            rows.append({"panel": panel, "cohort": cohort, "set": "evaluation" if panel == "a" else "at risk",
                         "n": values.pop(), "source_file": path, "source_column": column})
    return pd.DataFrame(rows)


S3_LAYER = f"{DERIVED}/s3a_microarray_layer.tsv"
S3_ABSORB = "results/stage8/diagnostics/esophagus_absorption.tsv"


def s3_fraction() -> pd.DataFrame:
    absorb = read(S3_ABSORB)
    layer = read(S3_LAYER).rename(columns={"esophagus_fraction": "esophagus_absorption"})
    rows = []
    for analysis, cohort in (("MET500", "MET500"), ("POG570", "POG570"), ("aux_rnaseq", "all"), ("microarray", "all")):
        frame = layer if analysis == "microarray" else absorb
        for method in ("BASE-Z", "SA-Z", "LD-Z"):
            hit = frame.loc[frame["analysis"].eq(analysis) & frame["cohort"].eq(cohort) & frame["method"].eq(method)]
            rows.append({"analysis": analysis, "cohort": cohort, "method": method,
                         "n_truth_not_esophagus": int(hit["n_truth_not_esophagus"].iloc[0]) if len(hit) else np.nan,
                         "esophagus_fraction": float(hit["esophagus_absorption"].iloc[0]) if len(hit) else np.nan,
                         "shown": "value" if len(hit) else "n/c",
                         "source": S3_LAYER if analysis == "microarray" else S3_ABSORB})
    return pd.DataFrame(rows)


def s3_gse41258() -> pd.DataFrame:
    colon = read("results/stage8/diagnostics/gse41258_colon_confusion.tsv")
    return colon.loc[colon["method"].eq("SA-Z")].sort_values("n", ascending=False, kind="mergesort").reset_index(drop=True)


def s3_platform() -> pd.DataFrame:
    mask = read("results/stage9/mask_metrics.tsv")
    return mask.loc[mask["rep"].eq("platform") & mask["method"].isin(["BASE-Z", "SA-Z"])].reset_index(drop=True)


def s3_random() -> pd.DataFrame:
    random = read("results/stage9/mask_random_summary.tsv")
    keep = random["method"].isin(["BASE-Z", "SA-Z"]) & random["mask"].ne("unmasked")
    return random.loc[keep].reset_index(drop=True)


def s4_contributions() -> pd.DataFrame:
    raw = read("results/stage9/imvigor_contribution_fix1.tsv")
    pieces = []
    for method, part in raw.groupby("method"):
        ordered = part.sort_values("mean_contribution", ascending=False, kind="mergesort")
        pieces.append(ordered.head(10).assign(side="positive", order=range(1, 11)))
        pieces.append(ordered.tail(5).assign(side="negative", order=range(1, 6)))
    shown = pd.concat(pieces, ignore_index=True)
    shown["display"] = shown["method"].map(label)
    return shown


def s5_reconciliation() -> pd.DataFrame:
    raw = read("results/stage9/esophagus_histology_summary.tsv")
    groups = [AUX, ["GSE209998"], MICROARRAY]
    rows = []
    for method in ("BASE-Z", "SA-Z", "LD-Z"):
        for cohorts in groups:
            for cohort in cohorts:
                hit = raw.loc[raw["method"].eq(method) & raw["cohort"].eq(cohort)]
                r = hit.iloc[0] if len(hit) else None
                rows.append({"cohort": COHORT.get(cohort, cohort), "method": label(method),
                             "n": int(r["n"]) if r is not None else 0,
                             "adenocarcinoma": int(r["n_closer_adenocarcinoma"]) if r is not None else 0,
                             "squamous": int(r["n_closer_squamous"]) if r is not None else 0,
                             "tie": int(r["n_tie"]) if r is not None else 0})
    return pd.DataFrame(rows)


def s6_meta() -> pd.DataFrame:
    raw = read("results/stage7/confirm/tables/meta_analysis.tsv")
    cohorts = raw.loc[raw["status"].eq("cohort")]
    rows = []
    for layer, title in (("RNA-seq", "RNA-seq"), ("\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774", "Microarray")):
        for r in cohorts.loc[cohorts["layer"].eq(layer)].itertuples():
            se = math.sqrt(float(r.v))
            rows.append({"layer": title, "row": "cohort", "cohort": COHORT.get(r.cohort, r.cohort),
                         "n_at_risk": int(r.n), "difference": float(r.d), "ci_low": float(r.d) - 1.96 * se,
                         "ci_high": float(r.d) + 1.96 * se, "interval": "Wald, recorded variance",
                         "I2": np.nan, "I2_percent": np.nan, "k": np.nan})
        hit = raw.loc[raw["status"].eq("\uacc4\uc0b0") & raw["layer"].eq(layer)].iloc[0]
        rows.append({"layer": title, "row": "DerSimonian-Laird", "cohort": "", "n_at_risk": np.nan,
                     "difference": float(hit["estimate"]), "ci_low": float(hit["ci_low"]),
                     "ci_high": float(hit["ci_high"]), "interval": "DerSimonian-Laird 95% interval",
                     "I2": float(hit["I2"]), "I2_percent": float(hit["I2"]) * 100, "k": int(hit["k"])})
    return pd.DataFrame(rows)


BUILDERS = {
    "Fig2_means": (fig2_means, ["results/stage5/sim_ext/full.tsv"]),
    "Fig2_loho": (fig2_loho, ["results/stage5/sim_ext/full.tsv"]),
    "Fig3_differences": (fig3_differences, ["results/stage4/confirm/primary.tsv",
                                            "results/stage7/confirm/tables/hypothesis_primary.tsv",
                                            "results/stage4/confirm/secondary_overall.tsv",
                                            "results/stage8/external/tables/metrics.tsv", "results/stage9/set_metrics.tsv"]),
    "Fig5_ablation": (fig5_ablation, ["results/stage9/ablation/metrics.tsv",
                                      "results/stage10/ablation_microarray/metrics.tsv",
                                      "results/stage10/ablation_microarray/restriction.tsv"]),
    "Fig6_published": (fig6_published, ["results/stage8/external/tables/metrics.tsv"]),
    "Fig4a_cohorts": (fig4a_cohorts, ["results/stage8/external/tables/metrics.tsv", "results/stage9/subcohort_metrics.tsv",
                                      "results/stage7/confirm/tables/metrics.tsv", f"{DERIVED}/microarray_rates.tsv"]),
    "Fig4b_site": (lambda: fig4b("site"), ["results/stage4/confirm/secondary_by_site.tsv"]),
    "Fig4b_tc": (lambda: fig4b("tc"), ["results/stage4/confirm/secondary_tc.tsv"]),
    "Fig4c_aux_top1": (lambda: copy("results/stage9/subcohort_metrics.tsv"), ["results/stage9/subcohort_metrics.tsv"]),
    "Fig4d_removed": (fig4d_removed, ["results/stage8/diagnostics/error_shift_cases.tsv"]),
    "S1_heatmap": (s1_heatmap, ["results/stage2/met500_overall.tsv", "results/stage4/confirm/secondary_overall.tsv",
                                "results/stage5/posthoc/met500_overall.tsv", "results/stage5/posthoc/pog_overall.tsv",
                                "results/stage5/posthoc/met500_rules.tsv", "results/stage5/posthoc/pog_rules.tsv",
                                f"{DERIVED}/microarray_rates.tsv", f"{DERIVED}/aux_six_rates.tsv"]),
    "S1_n": (s1_n, ["results/stage2/met500_overall.tsv", "results/stage4/confirm/secondary_overall.tsv",
                    f"{DERIVED}/aux_six_rates.tsv", f"{DERIVED}/microarray_rates.tsv"]),
    "S2_sets": (lambda: copy("results/stage9/set_metrics.tsv"), ["results/stage9/set_metrics.tsv"]),
    "S3_esophagus_fraction": (s3_fraction, [S3_ABSORB, S3_LAYER]),
    "S3_gse41258": (s3_gse41258, ["results/stage8/diagnostics/gse41258_colon_confusion.tsv"]),
    "S3_platform": (s3_platform, ["results/stage9/mask_metrics.tsv"]),
    "S3_random": (s3_random, ["results/stage9/mask_random_summary.tsv"]),
    "S4_contributions": (s4_contributions, ["results/stage9/imvigor_contribution_fix1.tsv"]),
    "S5_histology": (lambda: copy("results/stage9/esophagus_histology_summary.tsv"),
                     ["results/stage9/esophagus_histology_summary.tsv"]),
    "S5_reconciliation": (s5_reconciliation, ["results/stage9/esophagus_histology_summary.tsv"]),
    "S6_meta": (s6_meta, ["results/stage7/confirm/tables/meta_analysis.tsv"]),
}


def g10(value):
    return float(f"{float(value):.10g}")


def compare(name: str) -> dict:
    shown = pd.read_csv(FS / f"{name}.tsv", sep="\t")
    builder, sources = BUILDERS[name]
    expected = builder()
    out = {"file": f"{name}.tsv", "sources": "; ".join(sources), "rows": len(shown), "columns": len(shown.columns),
           "numeric_cells": 0, "string_cells": 0, "mismatch": 0, "first_mismatch": ""}
    if list(shown.columns) != list(expected.columns) or len(shown) != len(expected):
        out["mismatch"] = 1
        out["first_mismatch"] = f"shape {list(shown.columns)} x {len(shown)} vs {list(expected.columns)} x {len(expected)}"
        return out
    for col in shown.columns:
        for i in range(len(shown)):
            a, b = shown.at[i, col], expected.iloc[i][col]
            a_empty = a is None or (isinstance(a, float) and math.isnan(a)) or a == ""
            b_empty = b is None or (isinstance(b, float) and math.isnan(b)) or b == ""
            if isinstance(a, (int, float, np.integer, np.floating)) and not isinstance(a, bool) and not a_empty:
                out["numeric_cells"] += 1
                ok = not b_empty and g10(a) == g10(b)
            else:
                out["string_cells"] += 1
                ok = (a_empty and b_empty) or str(a) == str(b)
            if not ok:
                out["mismatch"] += 1
                if not out["first_mismatch"]:
                    out["first_mismatch"] = f"row {i} {col}: file {a!r} expected {b!r}"
    return out
