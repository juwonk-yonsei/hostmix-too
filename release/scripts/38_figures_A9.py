"""Finalize figure 5 and 6 and add post-hoc supplement figures. Does not retrain."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

def _release_config():
    cfg = Path(__file__).resolve().parents[1] / "config.yaml"
    out = {}
    for line in cfg.read_text().splitlines():
        if ":" not in line or line.strip().startswith("#"):
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if key in {"project_root", "figure_dir", "disk_mount"}:
            out[key] = value.strip().strip('"').strip("'")
    return out


def project_root():
    return Path(_release_config()["project_root"])


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = project_root()
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
figures = importlib.import_module("34_figures_A8")

FIG = figures.FIG
COLORS = figures.COLORS
DOUBLE = (175 / 25.4, 150 / 25.4)
SINGLE = (85 / 25.4, 90 / 25.4)


def tag(axis, text: str) -> None:
    axis.text(0.01, 0.99, text, transform=axis.transAxes, va="top", ha="left", fontsize=7)


def fig5() -> str:
    metrics = pd.read_csv(ROOT / "results/stage8/external/tables/metrics.tsv", sep="\t")
    metrics = metrics.loc[metrics["subset"].eq("common_label")].copy()
    sets = pd.read_csv(ROOT / "results/stage9/set_metrics.tsv", sep="\t")
    figures.source(metrics, "fig5_external.tsv")
    figures.source(sets, "fig5_sets.tsv")
    cohorts = ["MET500", "POG570", "aux_rnaseq"]
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    set_names = ["native_truth", "risk", "pool_out_risk"]
    fig, axes = plt.subplots(2, 2, figsize=DOUBLE)
    x = np.arange(len(cohorts))
    width = 0.18
    for axis, column in zip(axes[0], ["top1", "host_rate"]):
        for i, method in enumerate(methods):
            vals = [float(metrics.loc[metrics["analysis_cohort"].eq(c) & metrics["method"].eq(method), column].iloc[0]) for c in cohorts]
            axis.bar(x + (i - 1.5) * width, vals, width, label=method, color=COLORS[method])
        axis.set_xticks(x)
        axis.set_xticklabels(["MET500", "POG570", "aux"], fontsize=7)
        axis.set_ylim(0, 1)
        axis.set_ylabel(column)
        tag(axis, "pre-specified descriptive")
        axis.legend(fontsize=6, frameon=False)
    for axis, column in zip(axes[1], ["top1", "host_rate"]):
        xpos = np.arange(len(cohorts) * len(set_names))
        labels = []
        for i, method in enumerate(methods):
            vals = []
            for cohort in cohorts:
                for name in set_names:
                    row = sets.loc[sets["analysis_cohort"].eq(cohort) & sets["set"].eq(name) & sets["method"].eq(method)].iloc[0]
                    vals.append(float(row[column]))
                    if i == 0:
                        labels.append(f"{cohort[:3]}\n{name.split('_')[0]}\nn={int(row['n'])}")
            axis.bar(xpos + (i - 1.5) * width, vals, width, label=method, color=COLORS[method], hatch="//")
        axis.set_xticks(xpos)
        axis.set_xticklabels(labels, fontsize=5)
        axis.set_ylim(0, 1)
        axis.set_ylabel(column)
        tag(axis, "post hoc")
    fig.tight_layout()
    figures.save(fig, "fig5_external")
    return (
        "Figure 5. Top row is a pre-specified descriptive comparison, not a hypothesis test, "
        "from results/stage8/external/tables/metrics.tsv subset common_label. "
        "common_label n equals the full set. Top-1 denominators: MET500 437, POG570 512, aux RNA-seq selected_eval 729. "
        "Host-rate denominators: MET500 361, POG570 378, aux selected_risk 427 "
        "(eval intersect risk is 409; 18 GSE50760 risk-only samples are in the host-rate denominator only). "
        "Bottom row is post hoc (hatch) from results/stage9/set_metrics.tsv. "
        "Sets are native_truth, risk, and pool_out_risk. "
        "native_truth n: MET500 76, POG570 91, aux 315. "
        "pool_out_risk n: MET500 9, POG570 10, aux 71. "
        "Host rate is not structurally zero on native_truth when a site has more than one native organ."
    )


def fig6() -> str:
    cases = pd.read_csv(ROOT / "results/stage8/diagnostics/error_shift_cases.tsv", sep="\t")
    absorb = pd.read_csv(ROOT / "results/stage8/diagnostics/esophagus_absorption.tsv", sep="\t")
    colon = pd.read_csv(ROOT / "results/stage8/diagnostics/gse41258_colon_confusion.tsv", sep="\t")
    mask = pd.read_csv(ROOT / "results/stage9/mask_metrics.tsv", sep="\t")
    random_summary = pd.read_csv(ROOT / "results/stage9/mask_random_summary.tsv", sep="\t")
    figures.source(cases, "fig6_error_shift_cases.tsv")
    figures.source(absorb, "fig6_esophagus_absorption.tsv")
    figures.source(colon, "fig6_colon_confusion.tsv")
    platform = mask.loc[mask["rep"].eq("platform") & mask["method"].isin(["BASE-Z", "SA-Z"])].copy()
    control = random_summary.loc[random_summary["method"].isin(["BASE-Z", "SA-Z"]) & random_summary["mask"].ne("unmasked")].copy()
    figures.source(platform, "fig6_mask_platform.tsv")
    figures.source(control, "fig6_mask_random.tsv")
    fig, axes = plt.subplots(2, 2, figsize=DOUBLE)
    analyses = ["aux_rnaseq", "POG570", "MET500"]
    correct_n, other_n, favor_n = [], [], []
    for analysis in analyses:
        sub = cases.loc[cases["analysis"].eq(analysis)]
        flag = sub["SA_correct"].astype(str).str.lower().isin(["true", "1"])
        favor_n.append(int(len(sub)))
        correct_n.append(int(flag.sum()))
        other_n.append(int((~flag).sum()))
    xpos = np.arange(len(analyses))
    axes[0, 0].bar(xpos - 0.15, correct_n, 0.3, label="SA-Z correct", color=COLORS["SA-Z"], hatch="//")
    axes[0, 0].bar(xpos + 0.15, other_n, 0.3, label="other error", color="0.45", hatch="//")
    axes[0, 0].set_xticks(xpos)
    axes[0, 0].set_xticklabels([f"{name}\nn={n}" for name, n in zip(analyses, favor_n)], fontsize=6)
    axes[0, 0].set_ylabel("patients")
    tag(axes[0, 0], "post hoc")
    axes[0, 0].legend(fontsize=6, frameon=False)
    part = absorb.loc[absorb["cohort"].isin(["all", "POG570", "MET500"]) & absorb["method"].isin(["BASE-Z", "SA-Z", "LD-Z"])]
    axes[0, 1].bar(range(len(part)), part["esophagus_absorption"], color="0.45", hatch="//")
    axes[0, 1].set_xticks(range(len(part)))
    axes[0, 1].set_xticklabels([f"{a}\n{m}" for a, m in zip(part["analysis"], part["method"])], fontsize=5)
    axes[0, 1].set_ylabel("Esophagus fraction")
    tag(axes[0, 1], "post hoc")
    colon_sa = colon.loc[colon["method"].eq("SA-Z")].sort_values("n", ascending=False).head(6)
    axes[1, 0].bar(colon_sa["predicted"], colon_sa["n"], color="0.45", hatch="//")
    axes[1, 0].tick_params(axis="x", rotation=60)
    axes[1, 0].set_ylabel("patients")
    tag(axes[1, 0], "post hoc")
    platforms = ["GPL96", "GPL20769", "GPL15659", "prad_fhcrc_agilent"]
    short = ["GPL96", "GPL20769", "GPL15659", "Agilent"]
    x = np.arange(len(platforms))
    width = 0.18
    series = [("POG570", "BASE-Z"), ("POG570", "SA-Z"), ("MET500", "BASE-Z"), ("MET500", "SA-Z")]
    for i, (cohort, method) in enumerate(series):
        vals = []
        for name in platforms:
            row = platform.loc[platform["cohort"].eq(cohort) & platform["mask"].eq(name) & platform["method"].eq(method)].iloc[0]
            vals.append(float(row["top1_change"]))
        axes[1, 1].bar(x + (i - 1.5) * width, vals, width, label=f"{cohort[:3]} {method}", color=COLORS[method], hatch="//", alpha=0.9 if cohort == "POG570" else 0.45)
    axes[1, 1].axhline(0, color="black", linewidth=0.4)
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(short, fontsize=6)
    axes[1, 1].set_ylabel("top-1 change")
    axes[1, 1].legend(fontsize=5, frameon=False)
    tag(axes[1, 1], "post hoc")
    fig.tight_layout()
    figures.save(fig, "fig6_limits")
    return (
        "Figure 6 is post hoc (hatch). Panel A counts patients where BASE-Z was wrong and SA-Z was right, "
        "split by whether SA-Z matched truth. Source results/stage8/diagnostics/error_shift_cases.tsv. "
        "Panel B is the fraction of non-Esophagus truths called Esophagus. Source esophagus_absorption.tsv. "
        "Panel C is GSE41258 selected_eval colon_rectum, n=183, SA-Z predicted class counts. "
        "Panel D is the change in top-1 after keeping only genes measured on each microarray platform "
        "and rebuilding Z without retraining. Denominators stay 437 (MET500) and 512 (POG570). "
        "Random-drop means and ranges are in results/stage9/mask_random_summary.tsv, not drawn as bars. "
        "Unmasked Z matched stored BASE-Z and SA-Z labels on both cohorts (mismatch 0). "
        "Unmasked K matched stored POG570 labels (mismatch 0). MET500 K baseline is the unmasked rerun."
    )


def fig_s_imvigor() -> str:
    table = pd.read_csv(ROOT / "results/stage9/imvigor_contribution.tsv", sep="\t")
    figures.source(table, "figS1_imvigor_contribution.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(175 / 25.4, 120 / 25.4))
    for axis, method in zip(axes, ["BASE-Z", "SA-Z"]):
        part = table.loc[table["method"].eq(method)].sort_values("rank").head(10).iloc[::-1]
        axis.barh(part["gene"], part["mean_contribution"], color=COLORS[method], hatch="//")
        n = int(part["n_samples"].iloc[0])
        axis.set_xlabel("mean (ESCA coef - BLCA coef) x Z")
        axis.set_title(f"{method} n={n}")
        tag(axis, "post hoc")
    fig.tight_layout()
    figures.save(fig, "figS1_imvigor_contribution")
    return (
        "Figure S1 is post hoc. Mean logistic contribution on IMvigor210 bladder native-truth samples "
        "predicted Esophagus. BASE-Z n=145. SA-Z n=65. Native-truth denominator 194. "
        "Source results/stage9/imvigor_contribution.tsv. "
        "Length is GENCODE v23 exon-union kb, first positive row per symbol, in imvigor_length.tsv. "
        "Z comparison with TCGA-test BLCA n=81 is in imvigor_gene_Z.tsv. "
        "TCGA-test Z is rank-normal on all G genes, then B0 columns."
    )


def fig_s_esca() -> str:
    summary = pd.read_csv(ROOT / "results/stage9/esophagus_histology_summary.tsv", sep="\t")
    counts = pd.read_csv(ROOT / "results/stage9/esca_histology_counts.tsv", sep="\t")
    figures.source(summary, "figS2_esca_summary.tsv")
    figures.source(counts, "figS2_esca_counts.tsv")
    fig, axes = plt.subplots(1, 3, figsize=DOUBLE, sharey=True)
    for axis, method in zip(axes, ["BASE-Z", "SA-Z", "LD-Z"]):
        part = summary.loc[summary["method"].eq(method)].sort_values("n", ascending=False)
        y = np.arange(len(part))
        axis.barh(y, part["n_closer_adenocarcinoma"], color="#0072B2", hatch="//", label="closer to adenocarcinoma")
        axis.barh(y, part["n_closer_squamous"], left=part["n_closer_adenocarcinoma"], color="#D55E00", hatch="//", label="closer to squamous")
        axis.set_yticks(y)
        axis.set_yticklabels([f"{c} n={n}" for c, n in zip(part["cohort"], part["n"])], fontsize=5)
        axis.set_title(method)
        tag(axis, "post hoc")
    axes[0].legend(fontsize=5, frameon=False, loc="lower right")
    fig.tight_layout()
    figures.save(fig, "figS2_esca_histology")
    n_ad = int(counts.loc[counts["histology"].eq("adenocarcinoma"), "n"].iloc[0])
    n_sq = int(counts.loc[counts["histology"].eq("squamous"), "n"].iloc[0])
    n_ot = int(counts.loc[counts["histology"].eq("other"), "n"].iloc[0])
    return (
        f"Figure S2 is post hoc. TCGA-train ESCA n=145 from config/split_tcga.tsv. "
        f"GDC primary_diagnosis bins: adenocarcinoma {n_ad}, squamous {n_sq}, other {n_ot}. "
        "Adenosquamous and mixed strings are other. Centers use only pure adenocarcinoma and pure squamous. "
        "Each selected_eval sample predicted Esophagus is assigned by the larger Spearman correlation "
        "of its stored cohort Z to the two training-center mean Z vectors. Ties would be counted as ties. "
        "Source results/stage9/esophagus_histology_summary.tsv. "
        "Microarray Z was ranked on platform genes. TCGA centers were rank-normal on all G genes, then B0 columns."
    )


def main() -> None:
    notes = [
        figures.fig1(),
        figures.fig2(),
        figures.fig3(),
        figures.fig4(figures.rates_from_predictions()),
        fig5(),
        fig6(),
        fig_s_imvigor(),
        fig_s_esca(),
    ]
    header = (
        "Font is Liberation Sans 7-9 pt. Arial and Helvetica are not installed on this machine. "
        "Double-column width is 175 mm. PDF is vector. PNG is 300 dpi. "
        "Solid bars are preregistered or pre-specified descriptive, as labeled. Hatch marks post hoc.\n\n"
    )
    (FIG / "legends_en.md").write_text(header + "\n\n".join(notes) + "\n")
    print("figures A9", FIG, flush=True)


if __name__ == "__main__":
    main()
