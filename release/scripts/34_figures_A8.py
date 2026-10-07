"""Draft manuscript figures. Axis text names the quantity, the set, and n."""
from __future__ import annotations

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

# Arial and Helvetica are not installed. Liberation Sans is the fallback.
plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
DOUBLE_W = 175 / 25.4

sys.path.insert(0, str(project_root() / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
diagnostics = importlib.import_module("32_diagnostics_A8")

ROOT = project_root()
_figure = _release_config()["figure_dir"]
FIG = Path(_figure) if Path(_figure).is_absolute() else project_root() / _figure
SRC = FIG / "source"
COLORS = {
    "BASE-Z": "#0072B2", "SA-Z": "#D55E00", "SA-pool22": "#009E73", "LD-Z": "#666666",
    "SCOPE": "#CC79A7", "CUP-AI-Dx": "#56B4E9", "BASE-K": "#E69F00", "SA-K": "#F0E442",
    "SA-G": "#000000", "SA-MLP": "#882255", "SC-Z": "#44AA99", "NC-Z": "#88CCEE",
    "M1-Z": "#CC6677", "V0-K": "#DDCC77",
}
POOL_OUT = {"kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach", "colon_rectum", "prostate", "esophagus", "head_neck", "cervix"}


def save(fig, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.pdf")
    fig.savefig(FIG / f"{stem}.png", dpi=300)
    plt.close(fig)


def source(frame: pd.DataFrame, name: str) -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False, float_format="%.16g")


def fig1() -> str:
    raw = pd.read_csv(ROOT / "results/stage5/sim_ext/full.tsv", sep="\t")
    raw["pool"] = np.where(raw["site"].isin(POOL_OUT), "pool_out", "in_pool_site")
    methods = ["BASE-Z", "SA-Z", "SA-pool22"]
    rows = []
    for keys, sub in raw.groupby(["rho", "method", "pool"]):
        if keys[1] not in methods and not str(keys[1]).startswith("SA-LOHO"):
            continue
        rows.append({"rho": keys[0], "method": keys[1], "pool": keys[2], "host_pull_rate": float(sub["host_pull_rate"].mean()), "n_rows": int(len(sub))})
    table = pd.DataFrame(rows)
    source(table, "fig1_host_pull.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_W, 90 / 25.4), sharey=True)
    for method in methods:
        part = table.loc[table["method"].eq(method)].groupby("rho")["host_pull_rate"].mean()
        axes[0].plot(part.index, part.values, marker="o", color=COLORS[method], label=method)
    axes[0].set_title("mean across sites")
    for pool, style in (("in_pool_site", "-"), ("pool_out", "--")):
        part = table.loc[table["pool"].eq(pool) & table["method"].eq("SA-Z")].groupby("rho")["host_pull_rate"].mean()
        axes[1].plot(part.index, part.values, style, color=COLORS["SA-Z"], label=f"SA-Z {pool}")
        part = table.loc[table["pool"].eq(pool) & table["method"].eq("BASE-Z")].groupby("rho")["host_pull_rate"].mean()
        axes[1].plot(part.index, part.values, style, color=COLORS["BASE-Z"], label=f"BASE-Z {pool}")
    loho = table.loc[table["method"].str.startswith("SA-LOHO")].groupby("rho")["host_pull_rate"].mean()
    axes[1].plot(loho.index, loho.values, color="0.3", marker="s", label="SA-LOHO mean")
    for axis in axes:
        axis.set_xlabel("rho")
        axis.set_ylabel("host_pull_rate")
        axis.legend(fontsize=7)
    fig.tight_layout()
    save(fig, "fig1_simulation")
    return "fig1 uses results/stage5/sim_ext/full.tsv. Each point is the unweighted mean of host_pull_rate across tissue-site rows at that rho. in_pool_site is every site not in the pool-out list. SA-LOHO is the unweighted mean of SA-LOHO-* methods."


def fig2() -> str:
    pog = pd.read_csv(ROOT / "results/stage4/confirm/primary.tsv", sep="\t")
    aux = pd.read_csv(ROOT / "results/stage7/confirm/tables/hypothesis_primary.tsv", sep="\t")
    sec = pd.read_csv(ROOT / "results/stage7/confirm/tables/hypothesis_secondary.tsv", sep="\t")
    rows = []
    for _, row in pog.iterrows():
        rows.append({"label": f"POG {row['hypothesis']}", "diff": row["diff"], "low": row["ci_low"], "high": row["ci_high"], "n": row["n"], "kind": "preregistered"})
    for _, row in aux.iterrows():
        rows.append({"label": f"aux {row['hypothesis']}", "diff": row["diff_sa_minus_comparator"], "low": row["ci_low"], "high": row["ci_high"], "n": row["n"], "kind": "preregistered"})
    for _, row in sec.iterrows():
        rows.append({"label": f"aux {row['hypothesis']}", "diff": row["diff"], "low": row["ci_low"], "high": row["ci_high"], "n": row["n"], "kind": "preregistered"})
    table = pd.DataFrame(rows)
    source(table, "fig2_forest.tsv")
    fig, axis = plt.subplots(figsize=(DOUBLE_W, 110 / 25.4))
    y = np.arange(len(table))[::-1]
    axis.errorbar(table["diff"], y, xerr=[table["diff"] - table["low"], table["high"] - table["diff"]], fmt="o", color=COLORS["SA-Z"], capsize=3)
    axis.axvline(0, color="0.5", lw=0.8)
    axis.set_yticks(y)
    axis.set_yticklabels([f"{label} n={n}" for label, n in zip(table["label"], table["n"])])
    axis.set_xlabel("stored difference")
    fig.tight_layout()
    save(fig, "fig2_prereg_forest")
    return "fig2 uses stage4 confirm primary.tsv (column diff) and stage7 hypothesis_primary.tsv (column diff_sa_minus_comparator) and hypothesis_secondary.tsv (column diff). The drawn interval is the stored interval. AS contrasts are the stored secondary contrasts, not SA-Z minus BASE-Z."


def fig3() -> str:
    meta = pd.read_csv(ROOT / "results/stage7/confirm/tables/meta_analysis.tsv", sep="\t")
    source(meta, "fig3_meta.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_W, 90 / 25.4), sharex=False)
    for axis, layer in zip(axes, ["RNA-seq", "마이크로어레이"]):
        part = meta.loc[meta["layer"].eq(layer)]
        cohorts = part.loc[part["status"].eq("cohort")]
        pooled = part.loc[part["status"].eq("계산")]
        y = np.arange(len(cohorts))[::-1]
        axis.errorbar(cohorts["d"], y, xerr=[cohorts["d"] - (cohorts["d"] - 1.96 * np.sqrt(cohorts["v"])), (cohorts["d"] + 1.96 * np.sqrt(cohorts["v"])) - cohorts["d"]], fmt="o", color="0.3")
        if len(pooled):
            axis.axvline(float(pooled["estimate"].iloc[0]), color=COLORS["SA-Z"], label=f"DL {float(pooled['estimate'].iloc[0]):.3g}")
            axis.axvspan(float(pooled["ci_low"].iloc[0]), float(pooled["ci_high"].iloc[0]), color=COLORS["SA-Z"], alpha=0.15)
        axis.set_yticks(y)
        axis.set_yticklabels([f"{c} n={int(n)}" for c, n in zip(cohorts["cohort"], cohorts["n"])])
        axis.axvline(0, color="0.6", lw=0.6)
        axis.set_title({"RNA-seq": "RNA-seq", "마이크로어레이": "microarray"}.get(layer, layer))
        axis.legend(fontsize=7)
    fig.tight_layout()
    save(fig, "fig3_cohort_forest")
    return "fig3 uses results/stage7/confirm/tables/meta_analysis.tsv. Cohort intervals are d plus or minus 1.96 times the square root of the stored variance v. The vertical line is the stored DerSimonian-Laird estimate. n is the stored cohort n."


def rates_from_predictions() -> pd.DataFrame:
    frame = diagnostics.confirm()
    rows = []
    subsets = {
        "aux_rnaseq": frame.loc[frame["include_confirm"] & frame["selected_eval"] & frame["cohort"].isin(diagnostics.RNA)],
        "microarray": frame.loc[frame["layer"].eq("마이크로어레이") & frame["selected_eval"]],
    }
    methods = sorted({column[: -len("__pred")] for column in frame.columns if column.endswith("__pred")})
    for name, sub in subsets.items():
        risk = sub.loc[sub["selected_risk"]]
        for method in methods:
            pred = sub[f"{method}__pred"]
            host = [diagnostics.host(t, p, n) for t, p, n in zip(risk["organ"], risk[f"{method}__pred"], diagnostics.natives(risk["native_organs"]))]
            rows.append({
                "cohort": name, "method": method, "n": int(len(sub)), "top1": float((pred == sub["organ"]).mean()),
                "n_at_risk": int(len(risk)), "host_rate": float(np.mean(host)) if len(risk) else np.nan,
            })
    for cohort, path, truth_col, pred_prefix in (
        ("MET500", ROOT / "results/stage5/posthoc/met500_overall.tsv", None, None),
        ("POG570", ROOT / "results/stage5/posthoc/pog_overall.tsv", None, None),
    ):
        overall = pd.read_csv(path, sep="\t")
        for _, row in overall.iterrows():
            rows.append({"cohort": cohort, "method": row["method"], "n": row["n"], "top1": row["top1"], "n_at_risk": row["n_at_risk"], "host_rate": row["host_rate"]})
    return pd.DataFrame(rows)


def fig4(table: pd.DataFrame) -> str:
    source(table, "fig4_methods.tsv")
    methods = [m for m in ["BASE-Z", "SA-Z", "SA-G", "SA-pool22", "LD-Z", "SC-Z", "NC-Z", "M1-Z", "BASE-K", "SA-K", "V0-K", "SA-MLP"] if m in set(table["method"])]
    cohorts = ["MET500", "POG570", "aux_rnaseq", "microarray"]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_W, 100 / 25.4))
    for axis, column in zip(axes, ["top1", "host_rate"]):
        grid = table.pivot(index="method", columns="cohort", values=column).reindex(index=methods, columns=cohorts)
        image = axis.imshow(grid.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap="viridis")
        axis.set_xticks(range(len(cohorts)))
        axis.set_xticklabels(cohorts, rotation=30, ha="right")
        axis.set_yticks(range(len(methods)))
        axis.set_yticklabels(methods)
        axis.set_title(column)
        fig.colorbar(image, ax=axis, fraction=0.046)
    fig.tight_layout()
    save(fig, "fig4_method_heatmap")
    return "fig4 MET500 and POG570 are the stored overall tables. aux_rnaseq is include_confirm and selected_eval in the six RNA-seq cohorts. microarray is selected_eval in the microarray layer. host_rate uses the at-risk rows of the same selection. MET500 and POG570 rows that are not in the aux prediction file stay blank."


def fig5() -> str:
    path = ROOT / "results/stage8/external/tables/metrics.tsv"
    if not path.exists():
        return "fig5 not drawn; comparison table absent"
    table = pd.read_csv(path, sep="\t")
    table = table.loc[table["subset"].eq("common_label")]
    source(table, "fig5_external.tsv")
    cohorts = ["MET500", "POG570", "aux_rnaseq"]
    methods = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    x = np.arange(len(cohorts))
    width = 0.18
    for axis, column in zip(axes, ["top1", "host_rate"]):
        for i, method in enumerate(methods):
            vals = [float(table.loc[table["analysis_cohort"].eq(c) & table["method"].eq(method), column].iloc[0]) for c in cohorts]
            axis.bar(x + (i - 1.5) * width, vals, width, label=method, color=COLORS[method])
        axis.set_xticks(x)
        axis.set_xticklabels(cohorts)
        axis.set_ylim(0, 1)
        axis.set_title(column + " common_label")
        axis.legend(fontsize=7)
    fig.tight_layout()
    save(fig, "fig5_external")
    bits = []
    for cohort in cohorts:
        row = table.loc[table["analysis_cohort"].eq(cohort) & table["method"].eq("BASE-Z")].iloc[0]
        bits.append(f"{cohort} top1_n={int(row['n'])} host_n={int(row['n_at_risk'])}")
    return (
        "fig5 uses results/stage8/external/tables/metrics.tsv subset common_label. "
        "Pre-specified description, not a hypothesis test. "
        + "; ".join(bits)
        + ". common_label equals the full set in these three analysis cohorts. "
        "aux_rnaseq top-1 denominator is selected_eval 729. host_rate denominator is selected_risk 427, "
        "including 18 GSE50760 risk-only samples scored in a supplement run."
    )


def fig6() -> str:
    cases = pd.read_csv(ROOT / "results/stage8/diagnostics/error_shift_cases.tsv", sep="\t")
    absorb = pd.read_csv(ROOT / "results/stage8/diagnostics/esophagus_absorption.tsv", sep="\t")
    colon = pd.read_csv(ROOT / "results/stage8/diagnostics/gse41258_colon_confusion.tsv", sep="\t")
    source(cases, "fig6_error_shift_cases.tsv")
    source(absorb, "fig6_esophagus_absorption.tsv")
    source(colon, "fig6_colon_confusion.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    analyses = ["aux_rnaseq", "POG570", "MET500"]
    correct_n, other_n, favor_n = [], [], []
    for analysis in analyses:
        sub = cases.loc[cases["analysis"].eq(analysis)]
        flag = sub["SA_correct"].astype(str).str.lower().isin(["true", "1"])
        favor_n.append(int(len(sub)))
        correct_n.append(int(flag.sum()))
        other_n.append(int((~flag).sum()))
    xpos = np.arange(len(analyses))
    axes[0].bar(xpos - 0.15, correct_n, 0.3, label="SA-Z correct", color=COLORS["SA-Z"], hatch="//")
    axes[0].bar(xpos + 0.15, other_n, 0.3, label="other error", color="0.45", hatch="//")
    axes[0].set_xticks(xpos)
    axes[0].set_xticklabels([f"{name}\nn={n}" for name, n in zip(analyses, favor_n)], fontsize=7)
    axes[0].set_title("AH1-style favor")
    axes[0].legend(fontsize=6)
    part = absorb.loc[absorb["cohort"].isin(["all", "POG570", "MET500"]) & absorb["method"].isin(["BASE-Z", "SA-Z", "LD-Z"])]
    labels = [f"{a}\n{m}" for a, m in zip(part["analysis"], part["method"])]
    axes[1].bar(range(len(part)), part["esophagus_absorption"], color="0.45", hatch="//")
    axes[1].set_xticks(range(len(part)))
    axes[1].set_xticklabels(labels, fontsize=6, rotation=90)
    axes[1].set_title("Esophagus absorption")
    colon_sa = colon.loc[colon["method"].eq("SA-Z")].sort_values("n", ascending=False).head(6)
    axes[2].bar(colon_sa["predicted"], colon_sa["n"], color="0.45", hatch="//")
    axes[2].set_title(f"GSE41258 colon SA-Z n={int(colon_sa['n_slice'].iloc[0]) if len(colon_sa) else 0}")
    axes[2].tick_params(axis="x", rotation=60)
    fig.tight_layout()
    save(fig, "fig6_limits")
    return "fig6 is post hoc (hatch). Panel 1 is AH1-style favor counts: SA-Z correct versus other error. Class names are in the source TSV. Panel 2 is the fraction of non-Esophagus truths predicted Esophagus. Panel 3 is GSE41258 selected_eval colon_rectum, n=183."


def main() -> None:
    notes = [fig1(), fig2(), fig3(), fig4(rates_from_predictions()), fig5(), fig6()]
    (FIG / "legends.md").write_text("\n\n".join(notes) + "\n")
    print("figures", FIG, flush=True)


if __name__ == "__main__":
    main()
