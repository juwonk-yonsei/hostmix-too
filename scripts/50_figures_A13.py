"""Redraw manuscript figures for stage A13 from stored tables.

HOSTMIX is the only display name for the mixture model. Change it here,
then replace the same word in the manuscript, to rename the method.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

HOSTMIX = "HostMix-TOO"
BASE = "Baseline"

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.6,
})

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript/bmc/figures"
SRC = ROOT / "manuscript/bmc/figure_source"
MM = 1 / 25.4
FULL = 170 * MM
HALF = 85 * MM
BLUE = "#0072B2"
ORANGE = "#D55E00"
PINK = "#CC79A7"
GREEN = "#009E73"
GRAYS = ["#4D4D4D", "#737373", "#999999", "#B3B3B3"]

POOL = [
    "Liver", "Lung", "Brain - Cortex", "Adrenal Gland",
    "Skin - Not Sun Exposed (Suprapubic)", "Muscle - Skeletal",
    "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Spleen", "Whole Blood",
]
COHORT = {
    "blca_iatlas_imvigor210_2017": "IMvigor210",
    "brca_iatlas_anders_2022": "Anders",
    "mel_dfci_2019": "DFCI",
    "paad_iatlas_prince_2022": "PRINCE",
    "GSE50760": "GSE50760",
    "prad_su2c_2019": "SU2C/PCF",
    "GSE209998": "AURORA US",
    "GSE41258": "GSE41258",
    "GSE14018": "GSE14018",
    "GSE71729": "GSE71729",
    "GSE74685": "GSE74685",
    "prad_fhcrc": "FHCRC",
    "MET500": "MET500",
    "POG570": "POG570",
}
METHOD = {
    "BASE-Z": BASE, "SA-Z": HOSTMIX, "SC-Z": "Site-specific mixtures",
    "NC-Z": "Normal classes", "LD-Z": "Linear deconvolution", "M1-Z": "Native-organ masking",
    "IF20-Z": "Gene removal 20%", "SC+IF20-Z": "Site-specific + removal 20%",
    "SC+IF40-Z": "Site-specific + removal 40%",
    "BASE-K": "Gene sets", "SA-K": "Gene sets + mixtures", "V0-K": "Gene sets + host correction",
    "SA-G": "Gated model", "SA-pool22": "22-tissue model", "SA-pool3": "3-tissue model",
    "BASE-MLP": "Perceptron", "SA-MLP": "Perceptron + mixtures",
    "PURE-Zs": "Standardized, pure (C = 0.03)", "PURE-Zs-w": "Standardized, pure (C = 0.15)",
    "MIX-Z0": "Unstandardized, mixtures (C = 0.1)",
    "SCOPE": "SCOPE", "CUP-AI-Dx": "CUP-AI-Dx",
}
COLOR = {BASE: BLUE, HOSTMIX: ORANGE, "SCOPE": PINK, "CUP-AI-Dx": GREEN}


def label(method: str) -> str:
    if method in METHOD:
        return METHOD[method]
    if method.startswith("SA-LOHO-"):
        return f"Leave-one-host-out ({method[len('SA-LOHO-'):]})"
    return method


def letter(ax, text: str) -> None:
    ax.set_title(" ", fontsize=10, pad=8, loc="left")
    ax.text(-0.02, 1.08, text, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", clip_on=False)


def save(fig, stem: str, tight: bool = False) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    options = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(FIG / f"{stem}.pdf", **options)
    fig.savefig(FIG / f"{stem}.png", dpi=300, **options)
    plt.close(fig)


def source(frame: pd.DataFrame, name: str) -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False, float_format="%.10g")


def box(ax, xy, w, h, text, face, edge="0.2"):
    patch = FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                           facecolor=face, edgecolor=edge, linewidth=0.6, mutation_aspect=0.6)
    ax.add_patch(patch)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=6.5, wrap=True)


def fig1() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 78 * MM))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")
    ax = axes[0]
    letter(ax, "a")
    box(ax, (0.4, 6.2), 4.2, 2.6, "Liver biopsy\ncolorectal cells\n+ hepatocytes", "#F0F0F0")
    box(ax, (5.4, 6.2), 4.2, 2.6, "Primary-tumor\nclassifier", BLUE, "white")
    ax.annotate("", xy=(5.3, 7.5), xytext=(4.7, 7.5), arrowprops=dict(arrowstyle="->", lw=0.7))
    box(ax, (5.4, 3.2), 4.2, 2.2, 'Call: "Liver"\n(host attraction)', "#F4C7A8")
    ax.annotate("", xy=(7.5, 5.5), xytext=(7.5, 6.1), arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.text(0.4, 1.5, "At risk: true organ is outside\nthe site's native organs.", fontsize=6.5, va="top")
    ax.text(0.4, 0.4, "Host-attraction rate: fraction of\nat-risk biopsies called native.", fontsize=6.5, va="top")

    ax = axes[1]
    letter(ax, "b")
    box(ax, (0.3, 7.4), 4.4, 1.8, "TCGA-train\n7486 tumors", BLUE, "white")
    box(ax, (5.2, 7.4), 4.4, 1.8, "GTEx-ref\n10 tissues, 1103", "#E6E6E6")
    box(ax, (1.6, 4.6), 6.8, 2.0, "x = ρ t + (1−ρ) h\nρ ~ U(0.15, 1)\n1 pure + 4 mixtures", ORANGE, "white")
    ax.annotate("", xy=(4.8, 6.7), xytext=(2.5, 7.3), arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.annotate("", xy=(5.2, 6.7), xytext=(7.4, 7.3), arrowprops=dict(arrowstyle="->", lw=0.6))
    box(ax, (1.6, 1.6), 6.8, 2.2, "Rank-normal scores\n5000 genes\n32 labels → 26 organs", "#FFF2E6")
    ax.annotate("", xy=(5.0, 3.9), xytext=(5.0, 4.5), arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.text(0.3, 0.7, f"{BASE}: same features, no mixtures", fontsize=6.5)

    ax = axes[2]
    letter(ax, "c")
    stages = [
        (7.6, "Development\nsimulation, MET500", "#D6E6F5", "Development"),
        (5.8, "Lock", "#EEEEEE", "Lock"),
        (4.0, "Confirmation 1\nPOG570; H1–H3", "#F6D6C6", "Preregistered"),
        (2.2, "Confirmation 2\nsix RNA-seq; AH, AS", "#F6D6C6", "Preregistered"),
        (0.5, "Published tools\nAblation, errors", "#E5E5E5", "Descriptive / post hoc"),
    ]
    for y, text, face, _role in stages:
        box(ax, (0.6, y), 8.6, 1.5, text, face)
    source(pd.DataFrame([{"panel": "c", "stage": t, "role": r} for _y, t, _f, r in stages]), "Fig1_stages.tsv")
    fig.savefig(FIG / "Figure1.pdf")
    fig.savefig(FIG / "Figure1.png", dpi=300)
    plt.close(fig)


def tissue_mean(sim: pd.DataFrame, method: str, tissues: list[str]) -> pd.DataFrame:
    part = sim.loc[sim["method"].eq(method) & sim["tissue"].isin(tissues)]
    grouped = part.groupby(["rho", "tissue"], as_index=False)["host_pull_rate"].mean()
    out = grouped.groupby("rho", as_index=False)["host_pull_rate"].mean()
    out["method"] = method
    out["n_tissues"] = grouped.groupby("rho")["tissue"].nunique().to_numpy()
    return out


def display_tissue(tissue: str) -> str:
    proxy = {"Spleen": "Spleen (lymph node proxy)", "Whole Blood": "Whole Blood (bone marrow proxy)"}
    return proxy.get(tissue, tissue)


def fig2() -> None:
    sim = pd.read_csv(ROOT / "results/stage5/sim_ext/full.tsv", sep="\t")
    out_tissues = sorted(set(sim["tissue"]) - set(POOL))
    rows = []
    for method, tissues, pool in (
        ("BASE-Z", POOL, "in"), ("SA-Z", POOL, "in"),
        ("BASE-Z", out_tissues, "out"), ("SA-Z", out_tissues, "out"),
        ("SA-pool22", out_tissues, "out"),
    ):
        part = tissue_mean(sim, method, tissues)
        part["pool"] = pool
        rows.append(part)
    loho_rows = []
    for tissue in POOL:
        method = f"SA-LOHO-{tissue}"
        hit = sim.loc[sim["method"].eq(method) & sim["tissue"].eq(tissue) & sim["rho"].eq(0.6)]
        if len(hit) != 1:
            raise SystemExit(f"LOHO rows for {tissue}: {len(hit)}")
        base = sim.loc[sim["method"].eq("BASE-Z") & sim["tissue"].eq(tissue) & sim["rho"].eq(0.6)].iloc[0]
        mix = sim.loc[sim["method"].eq("SA-Z") & sim["tissue"].eq(tissue) & sim["rho"].eq(0.6)].iloc[0]
        loho_rows.append({
            "tissue": tissue, "display": display_tissue(tissue),
            "baseline": float(base["host_rate"]), "leave_one_host_out": float(hit["host_rate"].iloc[0]),
            "hostmix": float(mix["host_rate"]),
        })
    table = pd.concat(rows, ignore_index=True)
    loho = pd.DataFrame(loho_rows).sort_values("baseline", ascending=False)
    source(table, "Fig2_means.tsv")
    source(loho, "Fig2_loho.tsv")
    fig = plt.figure(figsize=(FULL, 150 * MM))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.45, wspace=0.28)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    specs = [
        (axes[0], "in", [("BASE-Z", BLUE, "o", BASE), ("SA-Z", ORANGE, "s", HOSTMIX)]),
        (axes[1], "out", [("BASE-Z", BLUE, "o", BASE), ("SA-Z", ORANGE, "s", HOSTMIX),
                          ("SA-pool22", "#666666", "D", "22-tissue model")]),
    ]
    for ax, pool, series in specs:
        for method, color, marker, name in series:
            part = table.loc[table["pool"].eq(pool) & table["method"].eq(method)].sort_values("rho")
            ax.plot(part["rho"], part["host_pull_rate"], color=color, marker=marker, ms=4, lw=1.1, label=name)
        ax.set_xlim(1.05, 0.15)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Tumor RNA fraction ρ")
        ax.legend(frameon=False, loc="upper left")
    axes[0].set_title("In-pool tissues", fontsize=8)
    axes[1].set_title("Out-of-pool tissues", fontsize=8)
    axes[0].set_ylabel("Host-attraction rate")
    letter(axes[0], "a")
    letter(axes[1], "b")
    ax = fig.add_subplot(grid[1, :])
    letter(ax, "c")
    y = np.arange(len(loho))
    ax.scatter(loho["baseline"], y, facecolors=BLUE, edgecolors=BLUE, marker="o", s=22, label=BASE, zorder=3)
    ax.scatter(loho["leave_one_host_out"], y, facecolors="none", edgecolors="black", marker="o", s=22, label="Leave-one-host-out", zorder=3)
    ax.scatter(loho["hostmix"], y, color=ORANGE, marker="s", s=18, label=HOSTMIX, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(loho["display"], fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Host-attraction rate at ρ = 0.6")
    ax.invert_yaxis()
    ax.legend(frameon=False, ncol=3, loc="upper right", bbox_to_anchor=(0.98, 0.42))
    fig.subplots_adjust(left=0.36, right=0.98, top=0.96, bottom=0.08, hspace=0.45, wspace=0.28)
    save(fig, "Figure2")


def fig3() -> None:
    rows = [
        ("a", "POG570", "H1", -0.2037037037037037, -0.24345238095238098, -0.164021164021164, None, "0.228 → 0.024 (n = 378)"),
        ("a", "Auxiliary", "AH1", -0.14754098360655737, -0.18266978922716628, -0.11241217798594849, None, "0.204 → 0.056 (n = 427)"),
        ("b", "POG570", "H2", 0.07421875, 0.0390625, 0.107421875, None, "0.682 → 0.756 (n = 512)"),
        ("b", "Auxiliary", "AH2", 0.17146776406035658, 0.1454046639231824, 0.20027434842249658, None, "0.527 → 0.698 (n = 729)"),
        ("c", "POG570", "H3", -0.06593406593406592, -0.13186813186813184, 0.0, -0.1208791208791209, "0.868 → 0.802 (n = 91)"),
        ("c", "Auxiliary", "AH3", 0.2571428571428571, 0.2095238095238095, 0.3047619047619048, 0.21587301587301588, "0.473 → 0.730 (n = 315)"),
    ]
    frame = pd.DataFrame(rows, columns=["panel", "cohort", "hypothesis", "diff", "low", "high", "bound", "side"])
    source(frame, "Fig3_differences.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 60 * MM), sharey=True)
    titles = {"a": "Host-attraction rate", "b": "Top-1 accuracy", "c": "Native-truth top-1"}
    names = {"POG570": "POG570", "Auxiliary": "Auxiliary RNA-seq"}
    for ax, panel in zip(axes, "abc"):
        part = frame.loc[frame["panel"].eq(panel)].reset_index(drop=True)
        y = np.arange(len(part))[::-1]
        ax.axvline(0, color="0.5", lw=0.6)
        ax.errorbar(part["diff"], y, xerr=[part["diff"] - part["low"], part["high"] - part["diff"]],
                    fmt="o", color="black", ms=4, lw=1, capsize=2)
        if panel == "c":
            ax.axvline(-0.10, color="0.35", ls="--", lw=0.7)
            ax.text(-0.10, 1.08, "margin", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6)
            ax.scatter(part["bound"], y, facecolors="none", edgecolors="black", marker="D", s=18, zorder=3, label="One-sided 95% lower bound")
            ax.set_xticks(sorted(set(list(ax.get_xticks()) + [-0.10])))
            ax.legend(frameon=False, fontsize=5.5, loc="lower left")
        for yi, side in zip(y, part["side"]):
            ax.text(1.04, yi, side, transform=ax.get_yaxis_transform(), va="center", fontsize=6, clip_on=False)
        ax.set_yticks(y)
        ax.set_yticklabels([names[name] for name in part["cohort"]])
        ax.set_title("")
        ax.text(0.14, 1.12, titles[panel], transform=ax.transAxes, fontsize=8, va="bottom", ha="left")
        ax.set_xlabel("Difference (HostMix-TOO − baseline)")
        letter(ax, panel)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.55, right=0.86)
    save(fig, "Figure3", tight=True)


def fig4() -> None:
    metrics = pd.read_csv(ROOT / "results/stage9/ablation/metrics.tsv", sep="\t")
    shares = pd.read_csv(ROOT / "results/stage9/ablation/shares.tsv", sep="\t")
    order = ["BASE-Z", "PURE-Zs", "PURE-Zs-w", "MIX-Z0", "SA-Z"]
    cohorts = ["MET500", "POG570", "aux_rnaseq"]
    names = {"MET500": "MET500", "POG570": "POG570", "aux_rnaseq": "Auxiliary"}
    rows = []
    for cohort in cohorts:
        for method in order:
            hit = metrics.loc[metrics["cohort"].eq(cohort) & metrics["method"].eq(method)].iloc[0]
            share = shares.loc[shares["cohort"].eq(cohort) & shares["share"].eq(method), "value"]
            rows.append({
                "cohort": names[cohort], "method": label(method), "internal": method,
                "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                "share": float(share.iloc[0]) if len(share) else np.nan,
            })
    frame = pd.DataFrame(rows)
    source(frame, "Fig4_ablation.tsv")
    fig, ax = plt.subplots(figsize=(FULL, 85 * MM))
    x = np.arange(len(cohorts))
    width = 0.15
    hatches = ["", "///", "...", "xx", ""]
    colors = [BLUE, GRAYS[0], GRAYS[1], GRAYS[2], ORANGE]
    for i, method in enumerate(order):
        part = frame.loc[frame["internal"].eq(method)]
        xpos = x + (i - 2) * width
        ax.bar(xpos, part["host_rate"], width=width * 0.92, color=colors[i],
               edgecolor="0.15", linewidth=0.4, hatch=hatches[i], label=label(method))
        if method in {"PURE-Zs", "PURE-Zs-w", "MIX-Z0"}:
            for xp, value, share in zip(xpos, part["host_rate"], part["share"]):
                ax.text(xp, value + 0.008, f"{float(share):.3f}", ha="center", va="bottom", fontsize=6, color="black")
    ax.set_xticks(x)
    ax.set_xticklabels([
        f"{names[c]}\n(n = {int(frame.loc[frame.cohort.eq(names[c]), 'n_at_risk'].iloc[0])})"
        for c in cohorts
    ])
    ax.set_ylabel("Host-attraction rate")
    ax.set_ylim(0, 0.40)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    save(fig, "Figure4", tight=True)


def fig5() -> None:
    raw = pd.read_csv(ROOT / "results/stage8/external/tables/metrics.tsv", sep="\t")
    part = raw.loc[raw["subset"].eq("common_label") & raw["method"].isin(["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"])]
    source(part.assign(display=part["method"].map(label)), "Fig5_published.tsv")
    cohorts = ["MET500", "POG570", "aux_rnaseq"]
    shown = {"aux_rnaseq": "Auxiliary"}
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    colors = [BLUE, ORANGE, PINK, GREEN]
    markers = ["o", "s", "D", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 78 * MM))
    x = np.arange(len(cohorts))
    width = 0.18
    for ax, column, title in zip(axes, ["top1", "host_rate"], ["Top-1 accuracy", "Host-attraction rate"]):
        for i, method in enumerate(methods):
            vals = []
            for cohort in cohorts:
                hit = part.loc[part["analysis_cohort"].eq(cohort) & part["method"].eq(method)].iloc[0]
                vals.append(float(hit[column]))
            ax.bar(x + (i - 1.5) * width, vals, width=width * 0.9, color=colors[i],
                   edgecolor="0.15", lw=0.4, label=label(method), hatch=["", "", "//", ".."][i])
        labels = []
        for cohort in cohorts:
            hit = part.loc[part["analysis_cohort"].eq(cohort) & part["method"].eq("BASE-Z")].iloc[0]
            labels.append(f"{shown.get(cohort, cohort)}\nn = {int(hit['n'])}; {int(hit['n_at_risk'])}")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, 1)
        ax.set_ylabel(title)
        ax.legend(frameon=False, fontsize=6)
    letter(axes[0], "a")
    letter(axes[1], "b")
    fig.tight_layout()
    save(fig, "Figure5")


def as_bool(frame: pd.DataFrame, column: str) -> pd.Series:
    if frame[column].dtype == bool:
        return frame[column]
    return frame[column].astype(str).str.lower().isin(["true", "1"])


def dedupe_patients(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or frame["patient_id"].is_unique:
        return frame.reset_index(drop=True)
    part = frame.copy()
    library = part["library"].fillna("")
    part["_library_rank"] = np.where(library.eq("polyA"), 0, 1)
    part = part.sort_values(["_library_rank", "sample_id", "cohort"], kind="mergesort")
    return part.drop_duplicates("patient_id", keep="first").drop(columns="_library_rank").reset_index(drop=True)


def pulled_rate(frame: pd.DataFrame, method: str) -> float:
    if frame.empty:
        return float("nan")
    pred = frame[f"{method}__pred"].to_numpy(dtype=object)
    truth = frame["organ"].to_numpy(dtype=object)
    natives = [set(filter(None, str(text).split("|"))) for text in frame["native_organs"]]
    pulled = [
        bool(natives[i]) and truth[i] not in natives[i] and pred[i] in natives[i]
        for i in range(len(frame))
    ]
    return float(np.mean(pulled))


def microarray_sets() -> tuple[pd.DataFrame, pd.DataFrame]:
    """One patient per set in the microarray descriptive layer. Post hoc aggregate."""
    columns = ["cohort", "patient_id", "sample_id", "library", "organ", "native_organs",
               "layer", "layer_test", "selected_eval", "selected_risk"]
    columns += [f"{method}__pred" for method in (
        "BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z",
        "BASE-K", "SA-K", "V0-K", "SA-MLP",
    )]
    pred = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet", columns=columns)
    array = pred.loc[pred["layer"].eq("마이크로어레이") & as_bool(pred, "layer_test")].copy()
    return (
        dedupe_patients(array.loc[as_bool(array, "selected_eval")]),
        dedupe_patients(array.loc[as_bool(array, "selected_risk")]),
    )


def cohort_rates() -> pd.DataFrame:
    rows = []
    external = pd.read_csv(ROOT / "results/stage8/external/tables/metrics.tsv", sep="\t")
    for cohort, group in (("MET500", "Development"), ("POG570", "Confirmation 1")):
        for method in ("BASE-Z", "SA-Z"):
            hit = external.loc[external["subset"].eq("common_label") & external["analysis_cohort"].eq(cohort) & external["method"].eq(method)].iloc[0]
            rows.append({"group": group, "cohort": cohort, "method": method,
                         "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                         "source": "stage8 metrics"})
    sub = pd.read_csv(ROOT / "results/stage9/subcohort_metrics.tsv", sep="\t")
    order = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
             "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019"]
    for cohort in order:
        for method in ("BASE-Z", "SA-Z"):
            hit = sub.loc[sub["cohort"].eq(cohort) & sub["method"].eq(method)].iloc[0]
            rows.append({"group": "Confirmation 2", "cohort": COHORT[cohort], "method": method,
                         "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                         "source": "subcohort_metrics"})
    metrics = pd.read_csv(ROOT / "results/stage7/confirm/tables/metrics.tsv", sep="\t")
    for method in ("BASE-Z", "SA-Z"):
        hit = metrics.loc[metrics["cohort"].eq("GSE209998") & metrics["standard_site"].eq("all") & metrics["method"].eq(method)].iloc[0]
        rows.append({"group": "Descriptive RNA-seq", "cohort": COHORT["GSE209998"], "method": method,
                     "host_rate": float(hit["host_rate"]), "n_at_risk": int(hit["n_at_risk"]),
                     "source": "stage7 metrics"})
    _, risk = microarray_sets()
    expected = {"GSE41258": 58, "GSE14018": 36, "GSE71729": 56, "GSE74685": 49, "prad_fhcrc": 8}
    for cohort, n_expected in expected.items():
        part = risk.loc[risk["cohort"].eq(cohort)]
        if len(part) != n_expected:
            raise SystemExit(f"{cohort} at-risk n {len(part)} != Table 1 {n_expected}")
        for method in ("BASE-Z", "SA-Z"):
            rows.append({"group": "Descriptive microarray", "cohort": COHORT[cohort], "method": method,
                         "host_rate": pulled_rate(part, method), "n_at_risk": int(len(part)),
                         "source": "post hoc aggregate of confirmatory predictions"})
    return pd.DataFrame(rows)


def fig6() -> None:
    rates = cohort_rates()
    source(rates, "Fig6a_cohorts.tsv")
    site = pd.read_csv(ROOT / "results/stage4/confirm/secondary_by_site.tsv", sep="\t")
    tc = pd.read_csv(ROOT / "results/stage4/confirm/secondary_tc.tsv", sep="\t")
    site = site.loc[site["method"].isin(["BASE-Z", "SA-Z"])]
    tc = tc.loc[tc["method"].isin(["BASE-Z", "SA-Z"])]
    source(site, "Fig6b_site.tsv")
    source(tc, "Fig6b_tc.tsv")
    sub = pd.read_csv(ROOT / "results/stage9/subcohort_metrics.tsv", sep="\t")
    source(sub, "Fig6c_aux_top1.tsv")
    shift = pd.read_csv(ROOT / "results/figures/source/fig6_error_shift_cases.tsv", sep="\t")
    # No sample identifiers are copied onward.
    counts = shift.groupby("analysis").agg(n=("SA_correct", "size"), correct=("SA_correct", "sum")).reset_index()
    counts["other"] = counts["n"] - counts["correct"]
    source(counts, "Fig6d_removed.tsv")

    fig = plt.figure(figsize=(FULL, 150 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    letter(ax, "a")
    cohorts = list(dict.fromkeys(rates["cohort"]))
    x = np.arange(len(cohorts))
    for method, color, marker in (( "BASE-Z", BLUE, "o"), ("SA-Z", ORANGE, "s")):
        part = rates.loc[rates["method"].eq(method)].set_index("cohort").loc[cohorts]
        ax.plot(x, part["host_rate"], color=color, marker=marker, ms=4, lw=1, label=label(method))
    ax.set_xticks(x)
    labels = []
    base = rates.loc[rates["method"].eq("BASE-Z")].set_index("cohort")
    for cohort in cohorts:
        labels.append(f"{cohort}\n({int(base.loc[cohort, 'n_at_risk'])})")
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Host-attraction rate")
    ax.legend(frameon=False, ncol=2)
    # group separators
    bounds = []
    groups = rates.loc[rates["method"].eq("BASE-Z"), "group"].tolist()
    for i in range(1, len(groups)):
        if groups[i] != groups[i - 1]:
            ax.axvline(i - 0.5, color="0.8", lw=0.6)

    ax = fig.add_subplot(gs[1, 0])
    letter(ax, "b")
    site_order = ["liver", "lung", "lymph_node", "soft_tissue", "remainder"]
    site_names = ["Liver", "Lung", "Lymph node", "Soft tissue", "Other"]
    xs = np.arange(len(site_order))
    for method, color, offset in (("BASE-Z", BLUE, -0.15), ("SA-Z", ORANGE, 0.15)):
        vals = [float(site.loc[site["method"].eq(method) & site["site_group"].eq(g), "host_rate"].iloc[0]) for g in site_order]
        ax.bar(xs + offset, vals, width=0.28, color=color, label=label(method))
    ax.set_xticks(xs)
    ax.set_xticklabels(site_names, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("Host-attraction rate")
    ax2 = ax.twinx()
    tc_order = ["T1", "T2", "T3"]
    # draw tumor-content as open markers on a second small axis would crowd.
    # Use a companion axis created below instead of twinx.
    ax2.remove()
    ax.legend(frameon=False, fontsize=6)

    ax = fig.add_subplot(gs[1, 1])
    letter(ax, "c")
    order = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
             "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019"]
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    colors = [BLUE, ORANGE, PINK, GREEN]
    xs = np.arange(len(order))
    width = 0.18
    for i, method in enumerate(methods):
        vals = [float(sub.loc[sub["cohort"].eq(c) & sub["method"].eq(method), "top1"].iloc[0]) for c in order]
        ax.bar(xs + (i - 1.5) * width, vals, width=width * 0.9, color=colors[i], label=label(method))
    ax.set_xticks(xs)
    ax.set_xticklabels([COHORT[c] for c in order], rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("Top-1 accuracy")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=5.5, ncol=2)

    # panel d replaces the unused twin. Rebuild grid: put d by splitting was wrong.
    # Add a new figure row by using a 3-row layout instead. Recreate below if d missing.
    fig.savefig(FIG / "Figure6_partial.pdf")
    plt.close(fig)
    fig6_full(rates, site, tc, sub, counts)


def fig6_full(rates, site, tc, sub, counts) -> None:
    fig = plt.figure(figsize=(FULL, 175 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.45, 1, 0.9], hspace=0.55, wspace=0.38)
    ax = fig.add_subplot(gs[0, :])
    letter(ax, "a")
    groups = ["Development", "Confirmation 1", "Confirmation 2", "Descriptive RNA-seq", "Descriptive microarray"]
    y_cohort = {}
    yticks = []
    ylabels = []
    cursor = 0
    for group in groups:
        yticks.append(cursor)
        ylabels.append(group)
        cursor += 1
        cohorts = list(dict.fromkeys(rates.loc[rates["group"].eq(group), "cohort"]))
        for cohort in cohorts:
            y_cohort[cohort] = cursor
            yticks.append(cursor)
            n = int(rates.loc[rates["cohort"].eq(cohort) & rates["method"].eq("BASE-Z"), "n_at_risk"].iloc[0])
            ylabels.append(f"{cohort} ({n})")
            cursor += 1
        cursor += 0.6
    for cohort, y in y_cohort.items():
        base = float(rates.loc[rates["cohort"].eq(cohort) & rates["method"].eq("BASE-Z"), "host_rate"].iloc[0])
        mix = float(rates.loc[rates["cohort"].eq(cohort) & rates["method"].eq("SA-Z"), "host_rate"].iloc[0])
        ax.plot([base, mix], [y, y], color="0.45", lw=0.8, zorder=1)
        ax.scatter([base], [y], color=BLUE, marker="o", s=16, zorder=2, label=label("BASE-Z") if cohort == "MET500" else None)
        ax.scatter([mix], [y], color=ORANGE, marker="s", s=16, zorder=2, label=label("SA-Z") if cohort == "MET500" else None)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    for tick, text in zip(ax.get_yticklabels(), ylabels):
        if text in groups:
            tick.set_fontweight("bold")
    ax.set_xlabel("Host-attraction rate")
    ax.set_ylim(cursor - 0.4, -0.8)
    ax.legend(frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[1, 0])
    letter(ax, "b")
    site_order = ["liver", "lung", "lymph_node", "soft_tissue", "remainder"]
    site_names = ["Liver", "Lung", "Lymph\nnode", "Soft\ntissue", "Other"]
    xs = np.arange(len(site_order))
    for method, color, offset in (("BASE-Z", BLUE, -0.15), ("SA-Z", ORANGE, 0.15)):
        vals = [float(site.loc[site["method"].eq(method) & site["site_group"].eq(g), "host_rate"].iloc[0]) for g in site_order]
        ax.bar(xs + offset, vals, width=0.28, color=color, label=label(method))
        for x, value in zip(xs + offset, vals):
            if value == 0:
                ax.text(x, 0.01, "0", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(xs)
    ax.set_xticklabels(site_names, fontsize=6)
    ax.set_ylabel("Host-attraction rate")
    ax.set_title("POG570: biopsy site", fontsize=7)
    ax.legend(frameon=False, fontsize=6)

    ax = fig.add_subplot(gs[1, 1])
    letter(ax, "b")
    # second half of b: tumor content. Use a prime so the pair reads as one panel's two cells.
    ax.texts[0].set_text("")
    tc_order = ["T1", "T2", "T3"]
    xs = np.arange(3)
    for method, color, offset in (("BASE-Z", BLUE, -0.15), ("SA-Z", ORANGE, 0.15)):
        vals = [float(tc.loc[tc["method"].eq(method) & tc["tc_bin"].eq(g), "host_rate"].iloc[0]) for g in tc_order]
        ax.bar(xs + offset, vals, width=0.28, color=color, label=label(method))
        for x, value in zip(xs + offset, vals):
            if value == 0:
                ax.text(x, 0.01, "0", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(xs)
    ax.set_xticklabels(["Low", "Mid", "High"], fontsize=7)
    ax.set_title("POG570: tumor-content tertile", fontsize=7)
    ax.set_ylabel("Host-attraction rate")

    ax = fig.add_subplot(gs[2, 0])
    letter(ax, "c")
    order = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
             "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019"]
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    colors = [BLUE, ORANGE, PINK, GREEN]
    xs = np.arange(len(order))
    width = 0.18
    for i, method in enumerate(methods):
        vals = [float(sub.loc[sub["cohort"].eq(c) & sub["method"].eq(method), "top1"].iloc[0]) for c in order]
        ax.bar(xs + (i - 1.5) * width, vals, width=width * 0.9, color=colors[i], label=label(method))
    ax.set_xticks(xs)
    ax.set_xticklabels([COHORT[c] for c in order], rotation=35, ha="right", fontsize=6)
    ax.set_ylabel("Top-1 accuracy")
    ax.set_ylim(0, 1.05)
    for text in ax.texts:
        if text.get_text() == "c":
            text.set_position((-0.18, 1.02))
    ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.55))

    ax = fig.add_subplot(gs[2, 1])
    letter(ax, "d")
    show = {"POG570": "POG570", "MET500": "MET500", "aux_rnaseq": "Auxiliary"}
    # analysis column values
    print("shift analyses", counts["analysis"].tolist(), counts[["n", "correct"]].to_string())
    names = []
    correct = []
    other = []
    for key, title in (("MET500", "MET500"), ("POG570", "POG570"), ("aux_rnaseq", "Auxiliary RNA-seq")):
        hit = counts.loc[counts["analysis"].astype(str).str.contains(key, case=False)]
        if hit.empty:
            continue
        row = hit.iloc[0]
        names.append(f"{title}\n(n = {int(row['n'])})")
        correct.append(int(row["correct"]))
        other.append(int(row["other"]))
    xs = np.arange(len(names))
    ax.bar(xs, correct, color=ORANGE, label="Correct")
    ax.bar(xs, other, bottom=correct, color="#BDBDBD", hatch="//", label="Other error")
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("Biopsies")
    ax.legend(frameon=False, fontsize=6)
    save(fig, "Figure6", tight=True)
    partial = FIG / "Figure6_partial.pdf"
    if partial.exists():
        partial.unlink()


def aux_six() -> pd.DataFrame:
    """Six confirmatory RNA-seq cohorts, one sample per patient within each set."""
    pred_cols = [f"{method}__pred" for method in (
        "BASE-Z", "SA-Z", "SC-Z", "LD-Z", "NC-Z", "M1-Z", "SA-G", "SA-pool22",
        "BASE-K", "SA-K", "V0-K", "SA-MLP",
    )]
    frame = pd.read_parquet(
        ROOT / "results/stage7/confirm/predictions.parquet",
        columns=["cohort", "layer", "organ", "native_organs", "selected_eval", "selected_risk", *pred_cols],
    )
    keep = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
            "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019"]
    part = frame.loc[frame["cohort"].isin(keep)]
    rows = []
    for method in [c[:-6] for c in pred_cols]:
        ev = part.loc[part["selected_eval"]]
        rows_top = float((ev[f"{method}__pred"] == ev["organ"]).mean())
        risk = part.loc[part["selected_risk"]]
        natives = [set(filter(None, str(text).split("|"))) for text in risk["native_organs"]]
        pred = risk[f"{method}__pred"].to_numpy()
        truth = risk["organ"].to_numpy()
        flags = [bool(natives[i]) and truth[i] not in natives[i] and pred[i] in natives[i] for i in range(len(risk))]
        rows.append({"cohort": "Auxiliary RNA-seq", "method": method, "top1": rows_top, "host_rate": float(np.mean(flags))})
    return pd.DataFrame(rows)


def locked_heatmap() -> pd.DataFrame:
    rows = []
    met = pd.read_csv(ROOT / "results/stage2/met500_overall.tsv", sep="\t")
    for _, row in met.iterrows():
        rows.append({"cohort": "MET500", "method": f"{row['method']}-{row['representation']}",
                     "top1": float(row["top1"]), "host_rate": float(row["host_rate"]),
                     "status": "exploratory"})
    pog = pd.read_csv(ROOT / "results/stage4/confirm/secondary_overall.tsv", sep="\t")
    for _, row in pog.iterrows():
        rows.append({"cohort": "POG570", "method": row["method"], "top1": float(row["top1"]),
                     "host_rate": float(row["host_rate"]), "status": "pre-specified descriptive"})
    for cohort, path in (
        ("MET500", ROOT / "results/stage5/posthoc/met500_overall.tsv"),
        ("POG570", ROOT / "results/stage5/posthoc/pog_overall.tsv"),
    ):
        overall = pd.read_csv(path, sep="\t")
        for method in ("SA-pool22", "SA-MLP"):
            hit = overall.loc[overall["method"].eq(method)].iloc[0]
            rows.append({"cohort": cohort, "method": method, "top1": float(hit["top1"]),
                         "host_rate": float(hit["host_rate"]), "status": "exploratory"})
    for cohort, path in (
        ("MET500", ROOT / "results/stage5/posthoc/met500_rules.tsv"),
        ("POG570", ROOT / "results/stage5/posthoc/pog_rules.tsv"),
    ):
        rules = pd.read_csv(path, sep="\t")
        hit = rules.loc[rules["rule"].eq("G-beta(0.02)")].iloc[0]
        rows.append({"cohort": cohort, "method": "SA-G", "top1": float(hit["top1"]),
                     "host_rate": float(hit["host_rate"]), "status": "exploratory"})
    evaluation, risk = microarray_sets()
    if len(evaluation) != 536 or len(risk) != 207:
        raise SystemExit(f"microarray sets {len(evaluation)} {len(risk)}")
    for method in ("BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z",
                   "BASE-K", "SA-K", "V0-K", "SA-MLP"):
        pred = evaluation[f"{method}__pred"].to_numpy(dtype=object)
        truth = evaluation["organ"].to_numpy(dtype=object)
        rows.append({"cohort": "Microarray", "method": method,
                     "top1": float(np.mean(pred == truth)), "host_rate": pulled_rate(risk, method),
                     "status": "post hoc aggregate"})
    return pd.concat([pd.DataFrame(rows), aux_six()], ignore_index=True)


def fig_s1() -> None:
    frame = locked_heatmap()
    methods = [
        "BASE-Z", "SA-Z", "SC-Z", "NC-Z", "LD-Z", "M1-Z", "IF20-Z", "SC+IF20-Z", "SC+IF40-Z",
        "BASE-K", "SA-K", "V0-K", "SA-G", "SA-pool22", "SA-MLP",
    ]
    cohorts = ["MET500", "POG570", "Auxiliary RNA-seq", "Microarray"]
    source(frame, "S1_heatmap.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 120 * MM))
    for ax, column, title in zip(axes, ["top1", "host_rate"], ["Top-1 accuracy", "Host-attraction rate"]):
        grid = np.full((len(methods), len(cohorts)), np.nan)
        for i, method in enumerate(methods):
            for j, cohort in enumerate(cohorts):
                hit = frame.loc[frame["method"].eq(method) & frame["cohort"].eq(cohort), column]
                if len(hit):
                    grid[i, j] = float(hit.iloc[0])
        im = ax.imshow(grid, cmap="cividis", vmin=0, vmax=1, aspect="auto")
        for i in range(len(methods)):
            for j in range(len(cohorts)):
                if np.isnan(grid[i, j]):
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="///", edgecolor="0.5", lw=0))
                    ax.text(j, i, "–", ha="center", va="center", fontsize=6, color="0.3")
                else:
                    rgba = im.cmap(im.norm(grid[i, j]))
                    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
                    ax.text(j, i, f"{grid[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                            color="black" if luminance > 0.55 else "white")
        ax.set_xticks(range(len(cohorts)))
        ax.set_xticklabels(["MET500", "POG570", "Auxiliary\nRNA-seq", "Microarray"], fontsize=6)
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels([label(m) for m in methods], fontsize=6)
        ax.set_title(title, fontsize=8)
    letter(axes[0], "a")
    letter(axes[1], "b")
    fig.tight_layout()
    save(fig, "FigureS1")


def fig_s2() -> None:
    raw = pd.read_csv(ROOT / "results/figures/source/fig5_sets.tsv", sep="\t")
    source(raw, "S2_sets.tsv")
    sets = ["native_truth", "risk", "pool_out_risk"]
    set_names = ["Native truth", "At risk", "Pool-out at risk"]
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    colors = [BLUE, ORANGE, PINK, GREEN]
    cohorts = ["MET500", "POG570", "aux_rnaseq"]
    cohort_names = {"MET500": "MET500", "POG570": "POG570", "aux_rnaseq": "Auxiliary"}
    fig, axes = plt.subplots(2, 3, figsize=(FULL, 110 * MM), sharey="row")
    for col, (set_name, title) in enumerate(zip(sets, set_names)):
        for row, column, ylab in ((0, "top1", "Top-1 accuracy"), (1, "host_rate", "Host-attraction rate")):
            ax = axes[row, col]
            xs = np.arange(len(cohorts))
            width = 0.18
            for i, method in enumerate(methods):
                vals = []
                for cohort in cohorts:
                    hit = raw.loc[raw["analysis_cohort"].eq(cohort) & raw["set"].eq(set_name) & raw["method"].eq(method)]
                    vals.append(float(hit[column].iloc[0]) if len(hit) and pd.notna(hit[column].iloc[0]) else np.nan)
                bars = ax.bar(xs + (i - 1.5) * width, vals, width=width * 0.9, color=colors[i], label=label(method) if col == 0 and row == 0 else None)
                for bar, value in zip(bars, vals):
                    if value == 0:
                        ax.text(bar.get_x() + bar.get_width() / 2, 0.02, "0", ha="center", va="bottom", fontsize=5)
            ax.set_xticks(xs)
            ax.set_xticklabels([cohort_names[c] for c in cohorts], fontsize=6)
            if row == 0:
                ax.set_title(title, fontsize=8)
            if col == 0:
                ax.set_ylabel(ylab)
            ax.set_ylim(0, 1.05)
    axes[0, 0].legend(frameon=False, fontsize=5.5, ncol=2)
    letter(axes[0, 0], "a")
    letter(axes[0, 1], "b")
    letter(axes[0, 2], "c")
    fig.tight_layout()
    save(fig, "FigureS2")


def fig_s3() -> None:
    absorb = pd.read_csv(ROOT / "results/figures/source/fig6_esophagus_absorption.tsv", sep="\t")
    colon = pd.read_csv(ROOT / "results/figures/source/fig6_colon_confusion.tsv", sep="\t")
    platform = pd.read_csv(ROOT / "results/figures/source/fig6_mask_platform.tsv", sep="\t")
    random = pd.read_csv(ROOT / "results/figures/source/fig6_mask_random.tsv", sep="\t")
    source(absorb, "S3_esophagus_fraction.tsv")
    source(colon, "S3_gse41258.tsv")
    source(platform, "S3_platform.tsv")
    source(random, "S3_random.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 85 * MM))
    # panel a: groups in the absorption file
    ax = axes[0]
    letter(ax, "a")
    groups = list(dict.fromkeys(absorb["analysis"]))
    methods = ["BASE-Z", "SA-Z", "LD-Z"]
    colors = [BLUE, ORANGE, "#666666"]
    xs = np.arange(len(groups))
    for i, method in enumerate(methods):
        vals = []
        for group in groups:
            hit = absorb.loc[absorb["analysis"].eq(group) & absorb["method"].eq(method) & absorb["cohort"].isin(["all", group, "POG570"])]
            if hit.empty:
                hit = absorb.loc[absorb["analysis"].eq(group) & absorb["method"].eq(method)]
            vals.append(float(hit["esophagus_absorption"].iloc[0]) if len(hit) else np.nan)
        ax.bar(xs + (i - 1) * 0.25, vals, width=0.24, color=colors[i], label=label(method))
    ax.set_xticks(xs)
    ax.set_xticklabels(["Auxiliary RNA-seq", "POG570", "MET500", "Microarray"], fontsize=6, rotation=25, ha="right")
    ax.set_ylabel("Esophagus call fraction")
    ax.legend(frameon=False, fontsize=6)
    ax = axes[1]
    letter(ax, "b")
    part = colon.loc[colon["method"].eq("SA-Z")].sort_values("n", ascending=False)
    ax.barh(part["predicted"][::-1], part["n"][::-1], color=ORANGE)
    ax.set_xlabel("GSE41258 native-truth calls")
    ax = axes[2]
    letter(ax, "c")
    platforms = ["GPL96", "GPL20769", "GPL15659", "prad_fhcrc_agilent"]
    platform_names = ["GPL96", "GPL20769", "GPL15659", "FHCRC Agilent"]
    series = [
        ("POG570", "BASE-Z", BLUE, "POG570 baseline"),
        ("POG570", "SA-Z", ORANGE, "POG570 HostMix-TOO"),
        ("MET500", "BASE-Z", "#56B4E9", "MET500 baseline"),
        ("MET500", "SA-Z", "#E69F00", "MET500 HostMix-TOO"),
    ]
    xs = np.arange(len(platforms))
    width = 0.18
    for i, (cohort, method, color, name) in enumerate(series):
        vals = []
        for mask in platforms:
            hit = platform.loc[platform["cohort"].eq(cohort) & platform["mask"].eq(mask) & platform["method"].eq(method)]
            vals.append(float(hit["top1_change"].iloc[0]))
        ax.bar(xs + (i - 1.5) * width, vals, width=width * 0.9, color=color, label=name)
        for j, mask in enumerate(platforms):
            rnd = random.loc[random["cohort"].eq(cohort) & random["mask"].eq(mask) & random["method"].eq(method)].iloc[0]
            origin = float(rnd["top1_mean"]) - float(rnd["top1_change_mean"])
            low = float(rnd["top1_min"]) - origin
            high = float(rnd["top1_max"]) - origin
            mean = float(rnd["top1_change_mean"])
            x = xs[j] + (i - 1.5) * width
            ax.plot([x - width * 0.35, x + width * 0.35], [mean, mean], color="black", lw=0.8, zorder=3)
            ax.vlines(x, low, high, color="0.15", lw=0.6, zorder=3)
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(["GPL96", "GPL20769", "GPL15659", "FHCRC Agilent"], fontsize=5.5, rotation=30, ha="right")
    ax.set_ylabel("Change in top-1 accuracy")
    ax.legend(frameon=False, fontsize=5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.42))
    fig.tight_layout()
    save(fig, "FigureS3", tight=True)


def fig_s4() -> None:
    raw = pd.read_csv(ROOT / "results/stage9/imvigor_contribution_fix1.tsv", sep="\t")
    pieces = []
    for method, part in raw.groupby("method"):
        ordered = part.sort_values("mean_contribution", ascending=False)
        pieces.append(ordered.head(10).assign(side="positive"))
        pieces.append(ordered.tail(5).assign(side="negative"))
    shown = pd.concat(pieces, ignore_index=True)
    source(shown, "S4_contributions.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 130 * MM))
    for ax, method, panel, color in (
        (axes[0], "BASE-Z", "a", BLUE),
        (axes[1], "SA-Z", "b", ORANGE),
    ):
        part = shown.loc[shown["method"].eq(method)].sort_values("mean_contribution", ascending=True)
        positive = part.loc[part["side"].eq("positive")]
        negative = part.loc[part["side"].eq("negative")]
        labels = list(negative["gene"]) + [""] + list(positive["gene"])
        values = list(negative["mean_contribution"]) + [np.nan] + list(positive["mean_contribution"])
        y = np.arange(len(labels))
        ax.barh(y, values, color=color)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6)
        n = int(part["n_samples"].iloc[0])
        ax.set_xlabel("Mean contribution to the esophagus − bladder logit")
        ax.set_title(f"{label(method)} (n = {n})", fontsize=8)
        letter(ax, panel)
    fig.tight_layout()
    save(fig, "FigureS4")


def fig_s5() -> None:
    raw = pd.read_csv(ROOT / "results/stage9/esophagus_histology_summary.tsv", sep="\t")
    order = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
             "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019", "GSE209998", "GSE41258",
             "GSE14018", "GSE71729", "GSE74685", "prad_fhcrc"]
    methods = ["BASE-Z", "SA-Z", "LD-Z"]
    source(raw, "S5_histology.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 145 * MM), sharey=True)
    y = np.arange(len(order))
    for ax, method in zip(axes, methods):
        adeno, squamous, ns = [], [], []
        for cohort in order:
            hit = raw.loc[raw["method"].eq(method) & raw["cohort"].eq(cohort)]
            if hit.empty:
                adeno.append(0)
                squamous.append(0)
                ns.append(0)
            else:
                adeno.append(int(hit["n_closer_adenocarcinoma"].iloc[0]))
                squamous.append(int(hit["n_closer_squamous"].iloc[0]))
                ns.append(int(hit["n"].iloc[0]))
        ax.barh(y, adeno, color="#E69F00", label="Adenocarcinoma")
        ax.barh(y, squamous, left=adeno, color="#56B4E9", label="Squamous")
        for yi, n in zip(y, ns):
            ax.text(max(adeno[yi] + squamous[yi], 0) + 2, yi, f"n = {n}", va="center", fontsize=5.5, clip_on=False)
        ax.set_yticks(y)
        ax.set_yticklabels([COHORT.get(c, c) for c in order], fontsize=6)
        ax.set_title("Linear\ndeconvolution" if method == "LD-Z" else label(method), fontsize=8)
        ax.set_xlim(0, 340)
        ax.set_xlabel("Esophagus calls" if method == "BASE-Z" else "")
        ax.invert_yaxis()
    axes[0].legend(frameon=False, fontsize=6, loc="lower center", bbox_to_anchor=(1.55, 1.02), ncol=2)
    letter(axes[0], "a")
    letter(axes[1], "b")
    letter(axes[2], "c")
    fig.tight_layout()
    save(fig, "FigureS5", tight=True)
    # reconciliation table
    rows = []
    for method in methods:
        for cohort in order:
            hit = raw.loc[raw["method"].eq(method) & raw["cohort"].eq(cohort)]
            if hit.empty:
                continue
            rows.append({
                "cohort": COHORT.get(cohort, cohort), "method": label(method),
                "n": int(hit["n"].iloc[0]),
                "adenocarcinoma": int(hit["n_closer_adenocarcinoma"].iloc[0]),
                "squamous": int(hit["n_closer_squamous"].iloc[0]),
            })
    source(pd.DataFrame(rows), "S5_reconciliation.tsv")


def fig_s6() -> None:
    raw = pd.read_csv(ROOT / "results/stage7/confirm/tables/meta_analysis.tsv", sep="\t")
    cohorts = raw.loc[raw["status"].eq("cohort")].copy()
    source(cohorts.drop(columns=[c for c in cohorts.columns if c == "note"], errors="ignore"), "S6_meta.tsv")
    se = np.sqrt(cohorts["v"].astype(float))
    cohorts["low"] = cohorts["d"] - 1.96 * se
    cohorts["high"] = cohorts["d"] + 1.96 * se
    fig, ax = plt.subplots(figsize=(FULL, 110 * MM))
    blocks = [("RNA-seq", "RNA-seq"), ("마이크로어레이", "Microarray")]
    y = 0
    yticks = []
    ylabels = []
    for layer, title in blocks:
        yticks.append(y)
        ylabels.append(title)
        y += 1
        part = cohorts.loc[cohorts["layer"].eq(layer)]
        for row in part.itertuples():
            name = COHORT.get(row.cohort, row.cohort)
            if row.cohort == "GSE209998":
                name += " (descriptive)"
            ax.errorbar(row.d, y, xerr=[[row.d - row.low], [row.high - row.d]],
                    fmt="o", color="black", ms=4, lw=1, capsize=2)
            yticks.append(y)
            n_risk = int(row.n) if hasattr(row, "n") else ""
            ylabels.append(f"{name} ({n_risk})")
            y += 1
        overall = raw.loc[raw["status"].eq("계산") & raw["layer"].eq(layer)]
        if len(overall) != 1:
            raise SystemExit(f"meta summary rows for {layer}: {len(overall)}")
        estimate = float(overall["estimate"].iloc[0])
        low = float(overall["ci_low"].iloc[0]) if "ci_low" in overall.columns else estimate
        high = float(overall["ci_high"].iloc[0]) if "ci_high" in overall.columns else estimate
        i2 = overall["I2"].iloc[0] if "I2" in overall.columns else overall.filter(like="I").iloc[0, 0]
        ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]],
                    fmt="D", color="black", ms=5, lw=1, capsize=2)
        ax.text(1.02, y, f"{estimate:.3f} [{low:.3f} to {high:.3f}]; I² = {float(i2):.2f}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=6, clip_on=False)
        yticks.append(y)
        ylabels.append("DerSimonian–Laird")
        y += 1.4
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    for tick, text in zip(ax.get_yticklabels(), ylabels):
        if text in ("RNA-seq", "Microarray"):
            tick.set_fontweight("bold")
    ax.set_ylim(y - 0.2, -0.8)
    ax.set_xlabel("Host-attraction difference")
    fig.tight_layout()
    save(fig, "FigureS6")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig_s1()
    fig_s2()
    fig_s3()
    fig_s4()
    fig_s5()
    fig_s6()
    print("figures written", FIG)


if __name__ == "__main__":
    raise SystemExit("57_figures_A15.py draws the figures; this module only supplies its loaders.")
