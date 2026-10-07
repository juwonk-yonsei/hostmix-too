"""Printed numbers of Figures 1-6 and S1-S6 (A17 section 4.5, A18 section 2.2).

Each figure has a list of expected printed items built from the result-based frames of
fig_expected_A17 (not from figure_source), from facts and constants, and from the tick
positions fixed in the drawing script DRAW. An item is a template with "{}" for each
number and one link id per number:
  G:<id>   a drawn value (frame cell or axis tick), resolved from CELLS
  F:<id>   a fact of manuscript/checks/facts.tsv
  C:<id>   a constant of manuscript/checks/constants.tsv

Value labels (bar and cell values, n, arrow and interval texts, I2, zero marks) carry the
element key that the drawing script wrote to manuscript/checks/figure_labels.tsv with the
extent of the text. The checker takes the PDF words inside that extent, requires their text to
equal the recorded text, and links their numbers to the item with the same key. Items without
a key (box text, titles with fixed words, legend constants) are matched to a printed line by
their fixed text; axis ticks are matched to a remaining printed token equal to the expected
tick text. Tick positions are read from the line of DRAW that sets them (code search inside
the drawing function), not from fixed line numbers.
"""
from __future__ import annotations

import html
import re
import subprocess
from dataclasses import dataclass, field
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd

import fig_expected_A17 as FE
import numcheck_A17 as NC

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "manuscript" / "bmc" / "figures"
DRAW = "scripts/84_figures_A19.py"
LABELS = ROOT / "manuscript" / "checks" / "figure_labels.tsv"
FIGURES = [f"Figure{i}" for i in range(1, 7)] + [f"FigureS{i}" for i in range(1, 7)]
SHORT = {**{f"Figure{i}": f"Fig{i}" for i in range(1, 7)}, **{f"FigureS{i}": f"FigS{i}" for i in range(1, 7)}}
FUNCTION = {**{f"Figure{i}": f"fig{i}" for i in range(1, 7)}, **{f"FigureS{i}": f"fig_s{i}" for i in range(1, 7)}}


def code_line(fig: str, pattern: str, nth: int = 0) -> int:
    """Line number in DRAW of the nth line inside the drawing function of fig that contains pattern."""
    lines = (ROOT / DRAW).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(f"def {FUNCTION[fig]}("))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("def ")), len(lines))
    hits = [i + 1 for i in range(start, end) if pattern in lines[i]]
    if len(hits) <= nth:
        raise ValueError(f"{DRAW}: {pattern!r} found {len(hits)} times in {FUNCTION[fig]}")
    return hits[nth]


def fmt(value: float, digits: int = 3) -> str:
    number = Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP)
    if number.is_zero():
        number = abs(number)
    return f"{number:.{digits}f}".replace("-", "\u2212")


def tick_text(value: float, digits: int | None) -> str:
    if digits is None:
        return str(int(value)) if float(value).is_integer() else repr(float(value))
    return "0" if abs(value) < 1e-9 else fmt(value, digits)


@dataclass
class Item:
    template: str
    links: list[str]
    whole: bool = True
    texts: list[str] = field(default_factory=list)
    key: str = ""


class Expect:
    """Items and drawn cells of one figure."""

    def __init__(self, fig: str) -> None:
        self.fig = fig
        self.items: list[Item] = []
        self.cells: dict[str, dict] = {}

    def cell(self, frame: str, row: str, col: str, value, sources: list[str], meaning: str,
             digits: int | None = 3, selector: str = "") -> str:
        ref = f"{self.fig}:{frame}[{row}].{col}"
        if digits is None:
            number = Decimal(int(value))
            display = str(int(value))
        else:
            number = Decimal(repr(float(value)))
            display = fmt(value, digits)
        src = sources[0]
        self.cells[ref] = {"value": number, "display": display, "raw": repr(value), "cells": set(),
                           "boundary": digits is not None and NC.boundary(number, digits),
                           "row": {"meaning": meaning, "source_file": "; ".join(sources),
                                   "selector": f"{frame}: {row}" + (f" ({selector})" if selector else ""),
                                   "column": col, "status": NC.status_of(src, selector), "cohort": "",
                                   "method": "", "decimals": "" if digits is None else str(digits)}}
        return "G:" + ref

    def ticks(self, axis: str, values, digits: int | None, pattern: str, nth: int = 0) -> None:
        line = code_line(self.fig, pattern, nth)
        for v in values:
            ref = f"{self.fig}:tick:{axis}[{tick_text(v, digits)}]"
            shown = tick_text(v, digits)
            self.cells[ref] = {"value": Decimal(NC.norm_token(shown)), "display": shown, "raw": repr(float(v)),
                               "cells": set(), "boundary": False,
                               "row": {"meaning": f"axis tick, {axis}", "source_file": DRAW,
                                       "selector": f"line {line}", "column": "tick", "status": "design",
                                       "cohort": "", "method": "", "decimals": ""}}
            self.token("G:" + ref)

    def line(self, template: str, *links: str, key: str = "") -> None:
        self.items.append(Item(template, list(links), True, key=self.full(key)))

    def token(self, link: str, key: str = "") -> None:
        self.items.append(Item("{}", [link], False, key=self.full(key)))

    def label(self, key: str, *links: str) -> None:
        """A value label: its numbers, in printed order, are the links; matched by key and position."""
        self.items.append(Item("", list(links), False, key=self.full(key)))

    def full(self, key: str) -> str:
        return f"{SHORT[self.fig]}/{key}" if key else ""


# ---------------------------------------------------------------- figure by figure

def fig1(e: Expect) -> None:
    e.line("{} primary tumors,", "F:tcga_train_n")
    e.line("{} labels", "F:n_labels")
    e.line("{} normal profiles,", "F:gtexref_pool_n")
    e.line("{} tissues", "C:host_pool_n")
    e.line("\u03c1 ~ U({}, {}); per tumor {} pure", "C:rho_mix_low", "C:rho_mix_high", "C:pure_per_tumor")
    e.line("+ {} mixtures ({} profiles)", "C:mixes_per_tumor", "F:sa_train_rows")
    e.line("Rank-normal scores, {} genes;", "F:sa_n_features")
    e.line("Multinomial logistic regression (C = {}):", "F:sa_c")
    e.line("{} labels \u2192 {} organs", "F:n_labels", "F:n_organs")
    e.line("{} RNA-seq cohorts", "F:aux_cohorts_n")
    e.line("unstandardized, C = {}", "F:base_c")


def fig2(e: Expect) -> None:
    src = ["results/stage5/sim_ext/full.tsv"]
    e.line("{} host-pool tissues", "C:host_pool_n")
    e.line("{} tissues outside the host pool", "F:sim_outside_n")
    e.line("Host-attraction rate at \u03c1 = {}", "C:rho_display")
    for panel in ("a", "b"):
        e.ticks(f"{panel} x", [1.0, 0.8, 0.6, 0.4, 0.2], 1, "FixedLocator([1.0, 0.8, 0.6, 0.4, 0.2])")
        e.ticks(f"{panel} y", np.linspace(0, 1, 6), 1, "ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))")
    e.ticks("c x", np.linspace(0, 1, 6), 1, "ax.xaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))")
    for r in FE.fig2_loho().itertuples():
        e.label(f"c/{r.tissue}/n_at_risk", e.cell("Fig2_loho", r.tissue, "n_at_risk", r.n_at_risk, src,
                                                  f"simulation at rho 0.6, {r.display}, at-risk n", None,
                                                  f"tissue == '{r.tissue}'"))


def fig3(e: Expect) -> None:
    src = FE.BUILDERS["Fig3_differences"][1]
    arm_src = {"a": ("results/stage4/confirm/secondary_overall.tsv", "results/stage9/set_metrics.tsv"),
               "b": ("results/stage4/confirm/secondary_overall.tsv", "results/stage8/external/tables/metrics.tsv"),
               "c": ("results/stage4/confirm/secondary_overall.tsv", "results/stage9/set_metrics.tsv")}
    for r in FE.fig3_differences().itertuples():
        pog = not r.hypothesis.startswith("A")
        source = arm_src[r.panel][0 if pog else 1]
        links = [e.cell("Fig3_differences", r.hypothesis, col, getattr(r, col), [source],
                        f"{r.hypothesis} arm, {name}", 3, "method == 'BASE-Z'" if col == "baseline" else
                        "method == 'SA-Z'") for col, name in (("baseline", "baseline"), ("hostmix", "HostMix-TOO"))]
        e.label(f"{r.panel}/{r.hypothesis}/arms", *links)
        hyp_src = "results/stage4/confirm/primary.tsv" if pog else "results/stage7/confirm/tables/hypothesis_primary.tsv"
        e.label(f"{r.panel}/{r.hypothesis}/n", e.cell("Fig3_differences", r.hypothesis, "n", r.n, [hyp_src],
                                                      f"{r.hypothesis}, n", None))
    e.ticks("x", [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3], 2, "FixedLocator([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3])")


def fig5(e: Expect) -> None:
    frame = FE.fig5_ablation()
    src = {"stage9 ablation": "results/stage9/ablation/metrics.tsv",
           "stage10 ablation": "results/stage10/ablation_microarray/metrics.tsv",
           "stage10 restriction": "results/stage10/ablation_microarray/restriction.tsv"}
    e.line("Standardized, pure (C = {})", "C:c_pure03")
    e.line("Standardized, pure (C = {})", "C:c_pure15")
    e.line("Unstandardized, mixtures (C = {})", "C:c_mix")
    for panel in ("a", "b"):
        part = frame.loc[frame["panel"].eq(panel)]
        for cohort in ["MET500", "POG570", "aux_rnaseq", "microarray"]:
            hit = part.loc[part["cohort"].eq(cohort)].iloc[0]
            e.label(f"{panel}/{cohort}/n", e.cell("Fig5_ablation", f"{panel}/{cohort}", "n", hit["n"],
                                                  [src[hit["source"]]], f"Figure 5{panel}, {hit['display']}, n", None))
    for r in frame.loc[frame["panel"].eq("c")].itertuples():
        e.label(f"c/{r.cohort}/{r.method}",
                e.cell("Fig5_ablation", f"c/{r.cohort}/{r.method}", "value", r.value, [src[r.source]],
                       f"{r.display}, {r.label}, mean top-1 loss under random gene removal", 3))
    tick = "ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, ymax, nticks)))"
    e.ticks("a y", np.linspace(0, 0.4, 5), 1, tick)
    e.ticks("b y", np.linspace(0, 1.0, 6), 1, tick)
    e.ticks("c y", np.linspace(0, 0.1, 6), 2, tick)


def fig6(e: Expect) -> None:
    frame = FE.fig6_published()
    src = ["results/stage8/external/tables/metrics.tsv"]
    for panel, col in (("a", "n"), ("b", "n_at_risk")):
        for cohort in FE.COHORTS3:
            hit = frame.loc[frame["analysis_cohort"].eq(cohort)].iloc[0]
            e.label(f"{panel}/{cohort}/n", e.cell("Fig6_published", cohort, col, hit[col], src, f"{cohort}, {col}",
                                                  None, f"analysis_cohort == '{cohort}'"))
    for r in frame.itertuples():
        e.label(f"b/{r.analysis_cohort}/{r.method}",
                e.cell("Fig6_published", f"{r.analysis_cohort}/{r.method}", "host_rate", r.host_rate, src,
                       f"{r.analysis_cohort}, {r.display}, host-attraction rate", 3,
                       f"analysis_cohort == '{r.analysis_cohort}' and method == '{r.method}'"))
    tick = "ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, top, 5 if top == 0.4 else 6)))"
    e.ticks("a y", np.linspace(0, 1.0, 6), 1, tick)
    e.ticks("b y", np.linspace(0, 0.4, 5), 1, tick)


def fig4(e: Expect) -> None:
    rates = FE.fig4a_cohorts()
    src_a = {"stage8 metrics": "results/stage8/external/tables/metrics.tsv",
             "subcohort_metrics": "results/stage9/subcohort_metrics.tsv",
             "stage7 metrics": "results/stage7/confirm/tables/metrics.tsv",
             "post hoc aggregate of confirmatory predictions": f"{FE.DERIVED}/microarray_rates.tsv"}
    for r in rates.loc[rates["method"].eq("BASE-Z")].itertuples():
        e.label(f"a/{r.cohort}/n_at_risk", e.cell("Fig4a_cohorts", r.cohort, "n_at_risk", r.n_at_risk,
                                                  [src_a[r.source]], f"{r.cohort}, at-risk n", None))
    e.ticks("a x", [0, 0.2, 0.4, 0.6, 0.8], 1, "FixedLocator([0, 0.2, 0.4, 0.6, 0.8])")
    for kind, key, order in (("site", "site_group", ["liver", "lung", "lymph_node", "soft_tissue", "remainder"]),
                             ("tc", "tc_bin", ["T1", "T2", "T3"])):
        frame = FE.fig4b(kind)
        src = [FE.BUILDERS[f"Fig4b_{kind}"][1][0]]
        for value in order:
            hit = frame.loc[frame[key].eq(value)].iloc[0]
            e.label(f"b/{kind}/{value}/n", e.cell(f"Fig4b_{kind}", value, "n_at_risk", hit["n_at_risk"], src,
                                                  f"POG570 {value}, at-risk n", None))
            for r in frame.loc[frame[key].eq(value)].itertuples():
                if float(r.host_rate) == 0:
                    e.label(f"b/{kind}/{value}/{r.method}/zero",
                            e.cell(f"Fig4b_{kind}", f"{value}/{r.method}", "host_rate_zero_mark", 0, src,
                                   f"POG570 {value}, {r.method}, host-attraction rate 0 (zero mark)", None,
                                   f"method == '{r.method}'"))
        e.ticks(f"b {kind} y", np.linspace(0, 0.5, 6), 1, "FixedLocator(np.linspace(0, 0.5, 6))")
    sub = FE.copy("results/stage9/subcohort_metrics.tsv")
    for cohort in FE.AUX:
        hit = sub.loc[sub["cohort"].eq(cohort)].iloc[0]
        e.label(f"c/{cohort}/n", e.cell("Fig4c_aux_top1", cohort, "n", hit["n"],
                                        ["results/stage9/subcohort_metrics.tsv"], f"{FE.COHORT[cohort]}, evaluation n",
                                        None))
    e.ticks("c y", np.linspace(0, 1, 6), 1, "ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))")
    counts = FE.fig4d_removed()
    src = ["results/stage8/diagnostics/error_shift_cases.tsv"]
    for analysis in ("MET500", "POG570", "aux_rnaseq"):
        hit = counts.loc[counts["analysis"].eq(analysis)].iloc[0]
        for col in ("correct", "other"):
            e.label(f"d/{analysis}/{col}", e.cell("Fig4d_removed", analysis, col, hit[col], src,
                                                  f"{analysis}, removed errors, {col}", None))
        e.label(f"d/{analysis}/n", e.cell("Fig4d_removed", analysis, "n", hit["n"], src,
                                          f"{analysis}, removed errors, n", None))
    e.ticks("d y", [0, 20, 40, 60, 80], None, "FixedLocator([0, 20, 40, 60, 80])")


S1_METHODS = ["BASE-Z", "SA-Z", "SC-Z", "NC-Z", "LD-Z", "M1-Z", "IF20-Z", "SC+IF20-Z", "SC+IF40-Z",
              "BASE-K", "SA-K", "V0-K", "SA-G", "SA-pool22", "SA-MLP"]
S1_COHORTS = ["MET500", "POG570", "Auxiliary RNA-seq", "Microarray"]
S1_SOURCE = {"MET500": "results/stage2/met500_overall.tsv", "POG570": "results/stage4/confirm/secondary_overall.tsv",
             "Microarray": f"{FE.DERIVED}/microarray_rates.tsv", "Auxiliary RNA-seq": f"{FE.DERIVED}/aux_six_rates.tsv"}


def fig_s1(e: Expect) -> None:
    frame = FE.s1_heatmap()
    e.line("Gene removal {}%", "C:if_q20")
    e.line("Site-specific + removal {}%", "C:if_q20")
    e.line("Site-specific + removal {}%", "C:if_q40")
    for panel, col in (("a", "top1"), ("b", "host_rate")):
        for method in S1_METHODS:
            for cohort in S1_COHORTS:
                hit = frame.loc[frame["method"].eq(method) & frame["cohort"].eq(cohort)]
                if not len(hit) or pd.isna(hit[col].iloc[0]):
                    continue
                source = S1_SOURCE[cohort]
                if method in {"SA-pool22", "SA-MLP"} and cohort in {"MET500", "POG570"}:
                    source = f"results/stage5/posthoc/{'met500' if cohort == 'MET500' else 'pog'}_overall.tsv"
                if method == "SA-G" and cohort in {"MET500", "POG570"}:
                    source = f"results/stage5/posthoc/{'met500' if cohort == 'MET500' else 'pog'}_rules.tsv"
                e.label(f"{panel}/{cohort}/{method}",
                        e.cell("S1_heatmap", f"{cohort}/{method}", col, hit[col].iloc[0], [source],
                               f"{cohort}, {FE.label(method)}, {col}", 3, f"method == '{method}'"))
    sizes = FE.s1_n()
    for r in sizes.itertuples():
        e.label(f"{r.panel}/{r.cohort}/n", e.cell("S1_n", f"{r.panel}/{r.cohort}", "n", r.n, [r.source_file],
                                                  f"{r.cohort}, {r.set} n", None, f"column {r.source_column}"))
    e.ticks("a colorbar", [0, 0.5, 1.0], 1, "bar.set_ticks(")
    e.ticks("b colorbar", [0, 0.1, 0.2], 1, "bar.set_ticks(")


def fig_s2(e: Expect) -> None:
    sets = FE.copy("results/stage9/set_metrics.tsv")
    src = ["results/stage9/set_metrics.tsv"]
    for col, set_name in enumerate(("native_truth", "risk", "pool_out_risk")):
        for row, (metric, upper) in enumerate((("top1", 1.0), ("host_rate", 0.4))):
            for cohort in FE.COHORTS3:
                part = sets.loc[sets["analysis_cohort"].eq(cohort) & sets["set"].eq(set_name)
                                & sets["method"].isin(FE.METHODS4)]
                panel = "abc"[col]
                e.label(f"{panel}/{set_name}/{metric}/{cohort}/n",
                        e.cell("S2_sets", f"{set_name}/{metric}/{cohort}", "n", part["n"].iloc[0], src,
                               f"{cohort}, {set_name} set, n", None))
                for r in part.itertuples():
                    if float(getattr(r, metric)) == 0:
                        e.label(f"{panel}/{set_name}/{metric}/{cohort}/{r.method}/zero",
                                e.cell("S2_sets", f"{set_name}/{cohort}/{r.method}", f"{metric}_zero_mark", 0, src,
                                       f"{cohort}, {set_name}, {r.method}, {metric} 0 (zero mark)", None))
            shown = (row == 0 and col == 0) or (row == 1 and col in (0, 1))
            if shown:
                e.ticks(f"{'abc'[col]} {metric} y", np.linspace(0, upper, 6 if upper == 1.0 else 5), 1,
                        "FixedLocator(np.linspace(0, upper, 6 if upper == 1.0 else 5))")


def fig_s3(e: Expect) -> None:
    drawn = FE.s3_fraction()
    for analysis, cohort in (("MET500", "MET500"), ("POG570", "POG570"), ("aux_rnaseq", "all"), ("microarray", "all")):
        part = drawn.loc[drawn["analysis"].eq(analysis) & drawn["cohort"].eq(cohort)]
        n = part["n_truth_not_esophagus"].dropna()
        e.label(f"a/{analysis}/n", e.cell("S3_esophagus_fraction", analysis, "n_truth_not_esophagus", n.iloc[0],
                                          [part["source"].iloc[0]],
                                          f"{analysis}, tumors whose truth is not esophagus", None))
    e.ticks("a y", np.linspace(0, 0.4, 5), 1, "FixedLocator(np.linspace(0, 0.4, 5))")
    colon = FE.s3_gse41258()
    csrc = ["results/stage8/diagnostics/gse41258_colon_confusion.tsv"]
    e.label("b/title/n", e.cell("S3_gse41258", "SA-Z", "n_slice", colon["n_slice"].iloc[0], csrc,
                                "GSE41258 native-truth slice, n", None))
    for r in colon.itertuples():
        e.label(f"b/{r.predicted}/n", e.cell("S3_gse41258", r.predicted, "n", r.n, csrc,
                                             f"GSE41258, called {r.predicted}, n", None))
    e.ticks("b x", [0, 25, 50, 75, 100], None, "FixedLocator([0, 25, 50, 75, 100])")
    platform = FE.s3_platform()
    psrc = ["results/stage9/mask_metrics.tsv"]
    for cohort in ("POG570", "MET500"):
        n = platform.loc[platform["cohort"].eq(cohort), "n"].iloc[0]
        for name in ("baseline", "HostMix-TOO"):
            e.label(f"c/{cohort}/{name}/n", e.cell("S3_platform", f"{cohort}/{name}", "n", n, psrc,
                                                   f"{cohort}, evaluation n", None))
    e.line("Mean of {} random removals", "C:random_removals")
    e.line("Minimum to maximum of {}", "C:random_removals")
    e.ticks("c y", [-0.15, -0.10, -0.05, 0.0], 2, "FixedLocator([-0.15, -0.10, -0.05, 0.0])")


def s4_auto_ticks(values: list[float], left_mm: float) -> list[str]:
    """Tick labels of an automatic x axis, reproduced with the axes size of the drawing script DRAW."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with matplotlib.rc_context({"font.family": "Liberation Sans", "font.size": 7, "xtick.labelsize": 7,
                                "ytick.labelsize": 7, "xtick.major.pad": 2, "pdf.fonttype": 42}):
        w, h = 170.0, 118.0
        fig = plt.figure(figsize=(w / 25.4, h / 25.4))
        ax = fig.add_axes([left_mm / w, 1 - (11.0 + 90.0) / h, 60.0 / w, 90.0 / h])
        low, high = min(values), max(values)
        pad = 0.06 * (high - low)
        ax.set_xlim(low - pad, high + pad)
        fig.canvas.draw()
        labels = [t.get_text() for t in ax.get_xticklabels()
                  if low - pad <= t.get_position()[0] <= high + pad and t.get_text()]
        plt.close(fig)
    return labels


def fig_s4(e: Expect) -> None:
    shown = FE.s4_contributions()
    src = ["results/stage9/imvigor_contribution_fix1.tsv"]
    line = code_line(e.fig, "ax.set_xlim(low - pad, high + pad)")
    for method, left, panel in (("BASE-Z", 21.0, "a"), ("SA-Z", 106.0, "b")):
        part = shown.loc[shown["method"].eq(method)]
        e.label(f"{panel}/{method}/n", e.cell("S4_contributions", method, "n_samples", part["n_samples"].iloc[0], src,
                                              f"IMvigor210, {FE.label(method)}, samples", None))
        for text in s4_auto_ticks([float(v) for v in part["mean_contribution"]], left):
            ref = f"{e.fig}:tick:{method} x auto[{text}]"
            e.cells[ref] = {"value": Decimal(NC.norm_token(text)), "display": text, "raw": text, "cells": set(),
                            "boundary": False,
                            "row": {"meaning": f"axis tick, {method} x (automatic locator, limits from the data)",
                                    "source_file": DRAW, "selector": f"line {line}", "column": "tick",
                                    "status": "design", "cohort": "", "method": "", "decimals": ""}}
            e.token("G:" + ref)


def fig_s5(e: Expect) -> None:
    table = FE.s5_reconciliation()
    src = ["results/stage9/esophagus_histology_summary.tsv"]
    internal = {"Baseline": ("BASE-Z", "a"), "HostMix-TOO": ("SA-Z", "b"), "Linear deconvolution": ("LD-Z", "c")}
    for r in table.itertuples():
        method, panel = internal[r.method]
        e.label(f"{panel}/{r.cohort}/{method}/n", e.cell("S5_reconciliation", f"{r.cohort}/{r.method}", "n", r.n, src,
                                                         f"{r.cohort}, {r.method}, esophagus calls", None))
    for k in range(3):
        e.ticks(f"{'abc'[k]} x", [0, 100, 200], None, "FixedLocator([0, 100, 200])")


def fig_s6(e: Expect) -> None:
    meta = FE.s6_meta()
    src = ["results/stage7/confirm/tables/meta_analysis.tsv"]
    for r in meta.itertuples():
        row = r.cohort if r.row == "cohort" else f"{r.layer}/pooled"
        links = [e.cell("S6_meta", row, col, getattr(r, col), src, f"{row}, {name}", 3)
                 for col, name in (("difference", "difference"), ("ci_low", "lower 95% limit"),
                                   ("ci_high", "upper 95% limit"))]
        e.label(f"{row}/interval", *links)
        if r.row == "cohort":
            e.label(f"{r.cohort}/n_at_risk", e.cell("S6_meta", row, "n_at_risk", r.n_at_risk, src,
                                                    f"{r.cohort}, at-risk n", None))
        else:
            e.label(f"{row}/I2", e.cell("S6_meta", row, "I2_percent", r.I2_percent, src,
                                        f"{r.layer}, I-squared (percent)", 1))
    e.ticks("x", [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2], 1, "FixedLocator([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2])")


BUILD = {"Figure1": fig1, "Figure2": fig2, "Figure3": fig3, "Figure4": fig4, "Figure5": fig5, "Figure6": fig6,
         "FigureS1": fig_s1, "FigureS2": fig_s2, "FigureS3": fig_s3, "FigureS4": fig_s4, "FigureS5": fig_s5,
         "FigureS6": fig_s6}


def expected(fig: str) -> Expect:
    e = Expect(fig)
    BUILD[fig](e)
    return e


def printed_lines(fig: str) -> list[str]:
    out = subprocess.run(["pdftotext", "-bbox-layout", str(PDF / f"{fig}.pdf"), "-"], capture_output=True,
                         text=True, check=True).stdout
    return [html.unescape(" ".join(re.findall(r"<word[^>]*>(.*?)</word>", m.group(1))))
            for m in re.finditer(r"<line[^>]*>(.*?)</line>", out, re.S)]


WORD_RX = re.compile(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>')


def printed_words(fig: str) -> list[list[tuple]]:
    """Words of each printed line with their extent in mm from the top-left of the page."""
    out = subprocess.run(["pdftotext", "-bbox-layout", str(PDF / f"{fig}.pdf"), "-"], capture_output=True,
                         text=True, check=True).stdout
    scale = 25.4 / 72
    return [[(float(a) * scale, float(b) * scale, float(c) * scale, float(d) * scale, html.unescape(t))
             for a, b, c, d, t in WORD_RX.findall(m.group(1))]
            for m in re.finditer(r"<line[^>]*>(.*?)</line>", out, re.S)]


def labels(fig: str, table: pd.DataFrame | None = None) -> pd.DataFrame:
    table = pd.read_csv(LABELS, sep="\t", dtype={"text": str}, keep_default_na=False) if table is None else table
    return table.loc[table["figure"].eq(fig)].reset_index(drop=True)


def template_rx(template: str) -> re.Pattern:
    parts = [re.escape(p.replace("\u2212", "-").replace(",", "")) for p in template.split("{}")]
    return re.compile("^" + r"(\S+?)".join(parts) + "$")
