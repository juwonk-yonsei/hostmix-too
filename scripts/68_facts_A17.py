"""Facts for the A17 number check (A17 section 4.6), written to manuscript/checks/facts.tsv.

The A14 facts of scripts/40_check_numbers.py keep their ids, files, selectors and columns; their
float display rules become Decimal transforms (int -> none/0, r3 -> none/3, r1pct -> x100/1,
p -> p). New facts come from family generators over named result files, and from single
definitions written next to the sentence that needs them. Every fact names its cohort, method,
metric and analysis set, and a one-line meaning. Status follows numcheck_A17.status_of.
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numcheck_A17 as NC  # noqa: E402

FACTS: list[dict] = []
SEEN: set[str] = set()

P_PRIMARY = "results/stage4/confirm/primary.tsv"
P_OVERALL = "results/stage4/confirm/secondary_overall.tsv"
P_SITE = "results/stage4/confirm/secondary_by_site.tsv"
P_TC = "results/stage4/confirm/secondary_tc.tsv"
P_AHP = "results/stage7/confirm/tables/hypothesis_primary.tsv"
P_AHS = "results/stage7/confirm/tables/hypothesis_secondary.tsv"
P_S7 = "results/stage7/confirm/tables/metrics.tsv"
P_MA = "results/stage7/confirm/tables/microarray_effects.tsv"
P_SENS = "results/stage7/confirm/tables/sensitivity.tsv"
P_META = "results/stage7/confirm/tables/meta_analysis.tsv"
P_EXT = "results/stage8/external/tables/metrics.tsv"
P_ABL = "results/stage9/ablation/metrics.tsv"
P_SHARE = "results/stage9/ablation/shares.tsv"
P_SUB = "results/stage9/subcohort_metrics.tsv"
P_SET = "results/stage9/set_metrics.tsv"
P_S10 = "results/stage10/ablation_microarray/metrics.tsv"
P_S10R = "results/stage10/ablation_microarray/restriction.tsv"
P_S10S = "results/stage10/ablation_microarray/shares.tsv"
P_S10C = "results/stage10/ablation_microarray/category.tsv"
P_SIM = "results/stage5/sim_ext/full.tsv"
P_TT = "results/stage5/tcga_test/overall.tsv"
P_TTP = "results/stage5/tcga_test/paired.tsv"
P_TM = "results/stage5/tcga_met/overall.tsv"
P_DEV = "results/stage2/met500_overall.tsv"
P_PHM = "results/stage5/posthoc/met500_overall.tsv"
P_PHP = "results/stage5/posthoc/pog_overall.tsv"
P_LOC = "results/stage8/diagnostics/leave_one_cohort.tsv"

POOL = ["Liver", "Lung", "Brain - Cortex", "Adrenal Gland", "Skin - Not Sun Exposed (Suprapubic)",
        "Muscle - Skeletal", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Spleen", "Whole Blood"]
OUT = ["Bladder", "Breast - Mammary Tissue", "Cervix - Ectocervix", "Colon - Transverse", "Esophagus - Mucosa",
       "Kidney - Cortex", "Minor Salivary Gland", "Ovary", "Pancreas", "Prostate", "Stomach", "Thyroid"]
LISTS = {"POOL": POOL, "OUT": OUT}

METHOD_KEY = {
    "BASE-Z": "baseline", "SA-Z": "HostMix-TOO", "SCOPE": "SCOPE", "CUP-AI-Dx": "CUP-AI-Dx",
    "SC-Z": "site-specific", "NC-Z": "normal classes", "LD-Z": "deconvolution", "M1-Z": "masking",
    "IF20-Z": "gene removal", "SC+IF20-Z": "gene removal", "SC+IF40-Z": "gene removal",
    "BASE-K": "gene sets", "SA-K": "gene sets", "V0-K": "gene sets",
    "SA-G": "gated", "SA-pool22": "22-tissue", "SA-pool3": "3-tissue",
    "BASE-MLP": "perceptron", "SA-MLP": "perceptron",
    "PURE-Zs": "pure C=0.03", "PURE-Zs-w": "pure C=0.15", "MIX-Z0": "mixture model",
}
TAG = {
    "BASE-Z": "base", "SA-Z": "sa", "SCOPE": "scope", "CUP-AI-Dx": "cup", "SC-Z": "sc", "NC-Z": "nc",
    "LD-Z": "ld", "M1-Z": "m1", "IF20-Z": "if20", "SC+IF20-Z": "scif20", "SC+IF40-Z": "scif40",
    "BASE-K": "basek", "SA-K": "sak", "V0-K": "v0k", "M1-K": "m1k", "NC-K": "nck", "LD-K": "ldk",
    "IF20-K": "if20k", "SC-K": "sck", "SC+IF20-K": "scif20k", "SC+IF40-K": "scif40k",
    "SA-G": "gated", "SA-pool22": "pool22", "SA-pool3": "pool3", "BASE-MLP": "basemlp", "SA-MLP": "samlp",
    "PURE-Zs": "pure03", "PURE-Zs-w": "pure15", "MIX-Z0": "mix", "SA-C0.01": "c001", "SA-C0.1": "c01",
    "SA-m0": "m0", "SA-m1": "mx1", "SA-m8": "m8", "SA-r30": "r30", "SA-r50": "r50",
}
TISSUE_TAG = {t: t.split(" ")[0].lower().replace("-", "") for t in POOL + OUT}
TISSUE_TAG.update({"Adipose - Subcutaneous": "adiposesc", "Adipose - Visceral (Omentum)": "adiposevis",
                   "Skin - Not Sun Exposed (Suprapubic)": "skin", "Brain - Cortex": "brain",
                   "Muscle - Skeletal": "muscle", "Whole Blood": "blood", "Adrenal Gland": "adrenal",
                   "Kidney - Cortex": "kidney", "Colon - Transverse": "colon", "Esophagus - Mucosa": "esophagus",
                   "Breast - Mammary Tissue": "breast", "Cervix - Ectocervix": "cervix",
                   "Minor Salivary Gland": "salivary"})
for t in POOL:
    TAG[f"SA-LOHO-{t}"] = f"loho{TISSUE_TAG[t]}"
    METHOD_KEY[f"SA-LOHO-{t}"] = "leave-one-host-out"
for m in ("SA-C0.01", "SA-C0.1", "SA-m0", "SA-m1", "SA-m8", "SA-r30", "SA-r50"):
    METHOD_KEY[m] = "variants"
COHORT_KEY = {"MET500": "MET500", "POG570": "POG570", "aux_rnaseq": "aux",
              "blca_iatlas_imvigor210_2017": "IMvigor210", "brca_iatlas_anders_2022": "Anders",
              "mel_dfci_2019": "DFCI", "paad_iatlas_prince_2022": "PRINCE", "GSE50760": "GSE50760",
              "prad_su2c_2019": "SU2C", "GSE209998": "AURORA", "GSE41258": "GSE41258", "GSE14018": "GSE14018",
              "GSE71729": "GSE71729", "GSE74685": "GSE74685", "prad_fhcrc": "FHCRC", "all": "microarray"}
CTAG = {"MET500": "met", "POG570": "pog", "aux_rnaseq": "aux", "blca_iatlas_imvigor210_2017": "imv",
        "brca_iatlas_anders_2022": "anders", "mel_dfci_2019": "dfci", "paad_iatlas_prince_2022": "prince",
        "GSE50760": "gse50760", "prad_su2c_2019": "su2c", "GSE209998": "aurora", "GSE41258": "gse41258",
        "GSE14018": "gse14018", "GSE71729": "gse71729", "GSE74685": "gse74685", "prad_fhcrc": "fhcrc",
        "all": "layer"}
METRIC = {
    "top1": ("top-1 accuracy", "evaluation", 3), "host_rate": ("host-attraction rate", "at risk", 3),
    "native_truth_top1": ("top-1 accuracy", "native truth", 3), "n": ("n", "evaluation", 0),
    "n_at_risk": ("n", "at risk", 0), "n_native_truth": ("n", "native truth", 0),
    "top3": ("top-3 accuracy", "evaluation", 3), "macro_f1": ("macro-F1", "evaluation", 3),
}
MTAG = {"top1": "top1", "host_rate": "host", "native_truth_top1": "native", "n": "n", "n_at_risk": "nrisk",
        "n_native_truth": "nnative", "top3": "top3", "macro_f1": "f1"}


def add(fid: str, src: str, sel: str, col: str, transform: str = "none", decimals: int = 3, form: str = "num",
        cohort: str = "", method: str = "", metric: str = "", set_: str = "", meaning: str = "",
        note: str = "", approx: str = "") -> None:
    if fid in SEEN:
        raise SystemExit(f"duplicate fact id {fid}")
    if not fid.isidentifier():
        raise SystemExit(f"fact id is not an identifier: {fid}")
    SEEN.add(fid)
    FACTS.append({"fact_id": fid, "source_file": src, "selector": sel, "column": col, "transform": transform,
                  "decimals": decimals, "format": form, "status": "", "cohort": cohort, "method": method,
                  "metric": metric, "set": set_, "meaning": meaning, "note": note, "approx": approx})


def metric_fact(fid: str, src: str, sel: str, col: str, cohort: str, method: str, where: str, pct: bool = False,
                note: str = "") -> None:
    name, set_name, places = METRIC[col]
    mkey = "" if name == "n" else METHOD_KEY.get(method, method)
    meaning = f"{where}, {mkey + ', ' if mkey else ''}{name}, {set_name} set"
    add(fid, src, sel, col, "x100" if pct else "none", 1 if pct else places,
        cohort=COHORT_KEY.get(cohort, cohort), method=mkey, metric=name,
        set_=set_name, meaning=meaning + (" (percent)" if pct else ""), note=note)


# ---------------------------------------------------------------- A14 facts, kept by id

def a14() -> None:
    old = importlib.import_module("40_check_numbers").build_facts()
    mapping = {"int": ("none", 0, "num"), "r3": ("none", 3, "num"), "r1pct": ("x100", 1, "num"), "p": ("none", 0, "p")}
    hyp_meta = {
        "H": ("POG570", P_PRIMARY), "AH": ("aux", P_AHP), "AS": ("aux", P_AHS),
    }
    hyp_set = {"H1": "at risk", "H2": "evaluation", "H3": "native truth", "AH1": "at risk", "AH2": "evaluation",
               "AH3": "native truth", "AS1": "native truth", "AS2": "pool-out at risk", "AS3": "evaluation"}
    hyp_metric = {"H1": "host-attraction rate", "H2": "top-1 accuracy", "H3": "top-1 accuracy",
                  "AH1": "host-attraction rate", "AH2": "top-1 accuracy", "AH3": "top-1 accuracy",
                  "AS1": "top-1 accuracy", "AS2": "host-attraction rate", "AS3": "top-1 accuracy"}
    hyp_method = {"AS1": "gated", "AS2": "22-tissue", "AS3": "gated"}
    words = {"n": "n", "diff": "difference", "ci_low": "95% CI lower end", "ci_high": "95% CI upper end",
             "low": "95% CI lower end", "high": "95% CI upper end", "n_favor": "discordant pairs in favor",
             "n_against": "discordant pairs against", "p": "one-sided exact McNemar p",
             "bound": "one-sided 95% lower bound"}
    for row in old.rows:
        fid = row["fact_id"]
        transform, places, form = mapping[row["transform"]]
        meta = {"cohort": "", "method": "", "metric": "", "set_": "", "meaning": ""}
        head = fid.split("_")[0]
        if head[:2] in {"AH", "AS"} or head[:1] == "H" and head[1:].isdigit():
            fam = head[:2] if head[:2] in {"AH", "AS"} else "H"
            cohort = hyp_meta[fam][0]
            part = fid[len(head) + 1:]
            what = words.get(part, part)
            meth = "HostMix-TOO" if fam != "AS" else hyp_method[head]
            meta = {"cohort": cohort, "method": meth, "metric": hyp_metric[head], "set_": hyp_set[head],
                    "meaning": f"{head}, {cohort}, {hyp_metric[head]}, {hyp_set[head]} set, {what}"}
        elif fid.startswith("pog_"):
            _, tag, met = fid.split("_")
            method = {"base": "BASE-Z", "sa": "SA-Z"}[tag]
            col = {"top1": "top1", "host": "host_rate", "native": "native_truth_top1", "n": "n",
                   "nrisk": "n_at_risk", "nnative": "n_native_truth"}[met]
            name, set_name, _ = METRIC[col]
            meta = {"cohort": "POG570", "method": "" if name == "n" else METHOD_KEY[method], "metric": name, "set_": set_name,
                    "meaning": f"POG570, {METHOD_KEY[method]}, {name}, {set_name} set"}
        elif fid.startswith(("ext_", "abl_")):
            parts = fid.split("_")
            cohort = "aux_rnaseq" if parts[1] == "aux" else parts[1]
            rest = parts[3:] if parts[1] == "aux" else parts[2:]
            tag, met = rest[0], rest[1]
            inv = {v: k for k, v in TAG.items()}
            method = inv[tag]
            col = {"top1": "top1", "host": "host_rate", "native": "native_truth_top1", "n": "n",
                   "nrisk": "n_at_risk", "nnative": "n_native_truth"}[met]
            name, set_name, _ = METRIC[col]
            where = "pre-specified comparison" if fid.startswith("ext_") else "post hoc ablation"
            meta = {"cohort": COHORT_KEY[cohort], "method": "" if name == "n" else METHOD_KEY[method], "metric": name, "set_": set_name,
                    "meaning": f"{COHORT_KEY[cohort]}, {METHOD_KEY[method]}, {name}, {set_name} set ({where})"}
        elif fid in {"loho_liver", "base_liver"}:
            meth = "leave-one-host-out" if fid == "loho_liver" else "baseline"
            meta = {"cohort": "sim", "method": meth, "metric": "host-attraction rate", "set_": "liver, rho 0.6",
                    "meaning": f"simulation, liver mixtures at rho 0.6, {meth}, host-attraction rate"}
        add(fid, row["source_file"], row["selector"], row["column"], transform, places, form,
            meta["cohort"], meta["method"], meta["metric"], meta["set_"], meta["meaning"], note="A14 fact")


# ---------------------------------------------------------------- families

def hypotheses() -> None:
    for h in ("H1", "H2", "H3"):
        sel = f"hypothesis == '{h}'"
        set_name = {"H1": "at risk", "H2": "evaluation", "H3": "native truth"}[h]
        metric = "host-attraction rate" if h == "H1" else "top-1 accuracy"
        for col, tag, what in (("diff", "diff_pct", "difference"), ("ci_low", "ci_low_pct", "95% CI lower end"),
                               ("ci_high", "ci_high_pct", "95% CI upper end"), ("onesided_low", "bound_pct", "one-sided 95% lower bound")):
            if col == "onesided_low" and h != "H3":
                continue
            add(f"{h}_{tag}", P_PRIMARY, sel, col, "x100", 1, cohort="POG570", method="HostMix-TOO", metric=metric,
                set_=set_name, meaning=f"{h}, POG570, {metric}, {set_name} set, {what} (percentage points)")
    for col, tag in (("top1", "top1"), ("host_rate", "host"), ("native_truth_top1", "native")):
        for method, mt in (("BASE-Z", "base"), ("SA-Z", "sa")):
            metric_fact(f"pog_{mt}_{tag}_pct", P_OVERALL, f"method == '{method}'", col, "POG570", method, "POG570", True)
    for h, src in (("AH1", P_AHP), ("AH2", P_AHP), ("AH3", P_AHP)):
        pass


def overall_pog() -> None:
    table = pd.read_csv(NC.ROOT / P_OVERALL, sep="\t")
    for method in table["method"]:
        if method in ARMS:
            continue
        for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth"):
            metric_fact(f"pog_{TAG[method]}_{MTAG[col]}", P_OVERALL, f"method == '{method}'", col, "POG570", method,
                        "POG570 confirmation")
    site = pd.read_csv(NC.ROOT / P_SITE, sep="\t")
    for (method, group) in site[["method", "site_group"]].itertuples(index=False):
        for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth"):
            metric_fact(f"pogsite_{group}_{TAG[method]}_{MTAG[col]}", P_SITE,
                        f"method == '{method}' and site_group == '{group}'", col, "POG570", method,
                        f"POG570 {group} biopsies")
    tc = pd.read_csv(NC.ROOT / P_TC, sep="\t")
    for (method, b) in tc[["method", "tc_bin"]].itertuples(index=False):
        for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth"):
            metric_fact(f"pogtc_{b}_{TAG[method]}_{MTAG[col]}", P_TC, f"method == '{method}' and tc_bin == '{b}'",
                        col, "POG570", method, f"POG570 tumor-content tertile {b}")


ARMS = {"BASE-Z", "SA-Z"}


def external() -> None:
    table = pd.read_csv(NC.ROOT / P_EXT, sep="\t")
    for (cohort, subset, method) in table[["analysis_cohort", "subset", "method"]].itertuples(index=False):
        if subset != "common_label":
            continue
        sel = f"subset == 'common_label' and analysis_cohort == '{cohort}' and method == '{method}'"
        prefix = f"ext_{cohort}_{TAG[method]}"
        for col in ("native_truth_top1", "n_native_truth"):
            metric_fact(f"{prefix}_{MTAG[col]}", P_EXT, sel, col, cohort, method, f"{COHORT_KEY[cohort]} comparison")
        for col in ("top1", "host_rate", "native_truth_top1"):
            metric_fact(f"{prefix}_{MTAG[col]}_pct", P_EXT, sel, col, cohort, method,
                        f"{COHORT_KEY[cohort]} comparison", True)
        if method == "SCOPE":
            add(f"{prefix}_top1ns", P_EXT, sel, "top1_ns_as_organ", cohort=COHORT_KEY[cohort], method="SCOPE",
                metric="top-1 accuracy", set_="evaluation",
                meaning=f"{COHORT_KEY[cohort]}, SCOPE, top-1 accuracy counting a normal-tissue call as its organ")
            add(f"{prefix}_nns", P_EXT, sel, "n_scope_ns_top1", decimals=0, cohort=COHORT_KEY[cohort], method="SCOPE",
                metric="count", set_="evaluation",
                meaning=f"{COHORT_KEY[cohort]}, SCOPE, number of top-1 calls of a normal-tissue class")


def subcohorts() -> None:
    table = pd.read_csv(NC.ROOT / P_SUB, sep="\t")
    for (cohort, method) in table[["cohort", "method"]].itertuples(index=False):
        sel = f"cohort == '{cohort}' and method == '{method}'"
        for col in ("n", "top1", "n_at_risk", "host_rate"):
            metric_fact(f"sub_{CTAG[cohort]}_{TAG[method]}_{MTAG[col]}", P_SUB, sel, col, cohort, method,
                        f"{COHORT_KEY[cohort]} (auxiliary cohort, post hoc split)")


def stage7_cohorts() -> None:
    table = pd.read_csv(NC.ROOT / P_S7, sep="\t")
    part = table.loc[table["standard_site"].eq("all")]
    for (cohort, method) in part[["cohort", "method"]].itertuples(index=False):
        if cohort == "all":
            continue
        sel = f"cohort == '{cohort}' and standard_site == 'all' and method == '{method}'"
        for col in ("n", "top1", "n_at_risk", "host_rate", "n_native_truth", "native_truth_top1"):
            metric_fact(f"s7_{CTAG[cohort]}_{TAG[method]}_{MTAG[col]}", P_S7, sel, col, cohort, method,
                        f"{COHORT_KEY[cohort]} (auxiliary confirmation run, per cohort)")


def ablation() -> None:
    table = pd.read_csv(NC.ROOT / P_ABL, sep="\t")
    for (cohort, method) in table[["cohort", "method"]].itertuples(index=False):
        sel = f"cohort == '{cohort}' and method == '{method}'"
        for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth"):
            fid = f"abl_{cohort}_{TAG[method]}_{MTAG[col]}"
            if fid in SEEN:
                continue
            metric_fact(fid, P_ABL, sel, col, cohort, method, f"{COHORT_KEY[cohort]} post hoc ablation")
    shares = pd.read_csv(NC.ROOT / P_SHARE, sep="\t")
    for (cohort, share) in shares[["cohort", "share"]].itertuples(index=False):
        add(f"share_{cohort}_{TAG[share]}", P_SHARE, f"cohort == '{cohort}' and share == '{share}'", "value",
            cohort=COHORT_KEY[cohort], method=METHOD_KEY[share], metric="share of host-attraction gap",
            set_="at risk", meaning=f"{COHORT_KEY[cohort]}, {METHOD_KEY[share]}, share of the baseline minus "
                                    f"HostMix-TOO host-attraction gap (post hoc ablation)")


def stage10() -> None:
    table = pd.read_csv(NC.ROOT / P_S10, sep="\t")
    for (model, set_, cohort, metric) in table[["model", "set", "cohort", "metric"]].itertuples(index=False):
        sel = f"model == '{model}' and set == '{set_}' and cohort == '{cohort}' and metric == '{metric}'"
        base = f"ma_{CTAG[cohort]}_{TAG[model]}"
        set_name = set_.replace("_", " ")
        mname = {"top1": "top-1 accuracy", "n_esophagus": "esophagus calls", "host_rate": "host-attraction rate",
                 "native_truth_top1": "top-1 accuracy"}[metric]
        places = 0 if metric == "n_esophagus" else 3
        where = "microarray layer" if cohort == "all" else f"microarray cohort {COHORT_KEY[cohort]}"
        add(f"{base}_{metric}", P_S10, sel, "value", decimals=places, cohort=COHORT_KEY[cohort],
            method=METHOD_KEY[model], metric=mname, set_=set_name,
            meaning=f"{where}, {METHOD_KEY[model]}, {mname}, {set_name} set (post hoc)")
        nid = f"ma_{CTAG[cohort]}_n_{set_}"
        if nid not in SEEN and model == "BASE-Z":
            add(nid, P_S10, sel, "n", decimals=0, cohort=COHORT_KEY[cohort], metric="n", set_=set_name,
                meaning=f"{where}, n of the {set_name} set (one sample per patient across the layer)")
        if metric in {"top1", "host_rate", "native_truth_top1"} and cohort == "all":
            add(f"{base}_{metric}_pct", P_S10, sel, "value", "x100", 1, cohort="microarray",
                method=METHOD_KEY[model], metric=mname, set_=set_name,
                meaning=f"microarray layer, {METHOD_KEY[model]}, {mname}, {set_name} set (percent, post hoc)")
    rest = pd.read_csv(NC.ROOT / P_S10R, sep="\t")
    for model in rest["model"].unique():
        for col in ("restriction_drop", "random_mean"):
            add(f"mg_{TAG[model]}_{col}_mean", P_S10R, f"model == '{model}'", "", f"mean:{col}",
                cohort="POG570 and MET500", method=METHOD_KEY[model], metric=f"mean {col}",
                set_="evaluation", meaning=f"missing-gene simulation, {METHOD_KEY[model]}, mean {col} over the "
                                           f"8 platform-cohort pairs (post hoc)")
        for (platform, cohort) in rest.loc[rest["model"].eq(model), ["platform", "cohort"]].itertuples(index=False):
            sel = f"model == '{model}' and platform == '{platform}' and cohort == '{cohort}'"
            for col in ("restriction_drop", "random_mean", "random_min", "random_max"):
                add(f"mg_{TAG[model]}_{platform.lower().replace(' ', '')}_{cohort.lower()}_{col}", P_S10R, sel, col,
                    cohort=cohort, method=METHOD_KEY[model], metric=col, set_="evaluation",
                    meaning=f"missing-gene simulation, {cohort}, platform {platform}, {METHOD_KEY[model]}, {col}")
    shares = pd.read_csv(NC.ROOT / P_S10S, sep="\t")
    for (contrast, name) in shares[["contrast", "name"]].itertuples(index=False):
        model = {"s_mix": "MIX-Z0", "s_std": "PURE-Zs", "s_stdw": "PURE-Zs-w"}.get(name)
        add(f"ms_{contrast}_{name}", P_S10S, f"contrast == '{contrast}' and name == '{name}'", "value",
            cohort="microarray", method=METHOD_KEY.get(model, name), metric="share of loss", set_="evaluation",
            meaning=f"microarray ablation, {contrast} contrast, share {name} of HostMix-TOO's loss (post hoc)")
    add("ms_accuracy_gap", P_S10S, "contrast == 'accuracy' and name == 's_mix'", "gap", cohort="microarray",
        metric="gap", set_="evaluation", meaning="microarray ablation, baseline minus HostMix-TOO top-1 gap (post hoc)")
    add("ms_random_gap", P_S10C, "contrast == 'random_removal'", "gap", cohort="POG570 and MET500",
        metric="gap", set_="evaluation",
        meaning="missing-gene simulation, HostMix-TOO minus baseline mean loss under random removal (post hoc)")
    cat = pd.read_csv(NC.ROOT / P_S10C, sep="\t")
    for name in ("s2_mix", "s2_std", "s2_stdw"):
        add(f"ms_random_{name}", P_S10C, "contrast == 'random_removal'", name, cohort="POG570 and MET500",
            method=METHOD_KEY[{"s2_mix": "MIX-Z0", "s2_std": "PURE-Zs", "s2_stdw": "PURE-Zs-w"}[name]],
            metric="share of loss", set_="evaluation",
            meaning=f"missing-gene simulation, share {name} of HostMix-TOO's extra loss under random removal (post hoc)")
    del cat


def simulation() -> None:
    sim = pd.read_csv(NC.ROOT / P_SIM, sep="\t")
    for method in ("BASE-Z", "SA-Z", "SA-pool22"):
        for rho in sorted(sim["rho"].unique()):
            for pool, name in (("POOL", "in"), ("OUT", "out")):
                sel = f"method == '{method}' and rho == {rho} and tissue in @{pool}"
                add(f"sim_{name}_{str(rho).replace('.', '')}_{TAG[method]}", P_SIM, sel, "", "mean:host_pull_rate",
                    cohort="sim", method=METHOD_KEY[method], metric="host-attraction rate",
                    set_=f"{'in-pool' if name == 'in' else 'pool-out'} tissues, rho {rho}",
                    meaning=f"simulation, {METHOD_KEY[method]}, mean host-attraction rate over the "
                            f"{'ten host-pool' if name == 'in' else 'twelve pool-out'} tissues at rho {rho}")
    for tissue in POOL + OUT:
        for method in ["BASE-Z", "SA-Z", "SA-pool22"] + ([f"SA-LOHO-{tissue}"] if tissue in POOL else []):
            sel = f"method == '{method}' and tissue == '{tissue}' and rho == 0.6"
            t = TISSUE_TAG[tissue]
            mt = "loho" if method.startswith("SA-LOHO") else TAG[method]
            for col, mtag in (("host_rate", "host"), ("top1", "top1"), ("n_at_risk", "nrisk"), ("n", "n")):
                add(f"simt_{t}_06_{mt}_{mtag}", P_SIM, sel, col, decimals=0 if col.startswith("n") else 3,
                    cohort="sim", method=METHOD_KEY[method], metric=METRIC[col][0], set_=f"{tissue}, rho 0.6",
                    meaning=f"simulation, {tissue} mixtures at rho 0.6, {METHOD_KEY[method]}, {METRIC[col][0]}")


def tcga() -> None:
    for src, cohort, tag in ((P_TT, "TCGA-test", "tt"), (P_TM, "TCGA-met", "tm")):
        table = pd.read_csv(NC.ROOT / src, sep="\t")
        for method in table["method"]:
            for col in ("top1", "n"):
                add(f"{tag}_{TAG[method]}_{MTAG[col]}", src, f"method == '{method}'", col,
                    decimals=0 if col == "n" else 3, cohort=cohort, method=METHOD_KEY[method],
                    metric=METRIC[col][0], set_="all", meaning=f"{cohort}, {METHOD_KEY[method]}, {METRIC[col][0]}")
    for col, what in (("diff", "difference"), ("ci_low", "95% CI lower end"), ("ci_high", "95% CI upper end"),
                      ("n", "n"), ("p_two_sided", "two-sided p")):
        add(f"tt_paired_{col}", P_TTP, "comparison == 'SA-Z - BASE-Z'", col,
            decimals=0 if col == "n" else 3, form="p" if col.startswith("p_") else "num", cohort="TCGA-test",
            method="HostMix-TOO", metric="top-1 accuracy", set_="all",
            meaning=f"TCGA-test, HostMix-TOO minus baseline top-1 accuracy, {what}")


def development() -> None:
    table = pd.read_csv(NC.ROOT / P_DEV, sep="\t")
    for (method, rep) in table[["method", "representation"]].itertuples(index=False):
        mid = f"{method}-{rep}"
        sel = f"method == '{method}' and representation == '{rep}'"
        for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth", "top3",
                    "macro_f1", "n_host"):
            if col == "n_host":
                add(f"dev_{TAG[mid]}_nhost", P_DEV, sel, col, decimals=0, cohort="MET500", method=METHOD_KEY.get(mid, mid),
                    metric="host-attraction errors", set_="at risk",
                    meaning=f"MET500 development, {METHOD_KEY.get(mid, mid)}, number of host-attraction errors")
                continue
            metric_fact(f"dev_{TAG[mid]}_{MTAG[col]}", P_DEV, sel, col, "MET500", mid, "MET500 development")


def variants() -> None:
    for src, cohort, tag in ((P_PHM, "MET500", "phm"), (P_PHP, "POG570", "php")):
        table = pd.read_csv(NC.ROOT / src, sep="\t")
        for method in table["method"]:
            for col in ("top1", "host_rate", "native_truth_top1", "n", "n_at_risk", "n_native_truth"):
                metric_fact(f"{tag}_{TAG[method]}_{MTAG[col]}", src, f"method == '{method}'", col, cohort, method,
                            f"{cohort} exploratory variants")


def leave_one_cohort() -> None:
    table = pd.read_csv(NC.ROOT / P_LOC, sep="\t")
    for (dropped, h) in table[["dropped", "hypothesis"]].itertuples(index=False):
        sel = f"dropped == '{dropped}' and hypothesis == '{h}'"
        d = CTAG.get(dropped, dropped)
        for col, what in (("n", "n"), ("diff_sa_minus_base", "difference"), ("ci_low", "95% CI lower end"),
                          ("ci_high", "95% CI upper end"), ("onesided_p05", "one-sided 95% lower bound")):
            tag = {"n": "n", "diff_sa_minus_base": "diff", "ci_low": "low", "ci_high": "high", "onesided_p05": "bound"}[col]
            add(f"loc_{d}_{h}_{tag}", P_LOC, sel, col, decimals=0 if col == "n" else 3,
                cohort="aux", method="HostMix-TOO", metric={"AH1": "host-attraction rate"}.get(h, "top-1 accuracy"),
                set_={"AH1": "at risk", "AH2": "evaluation", "AH3": "native truth"}[h],
                meaning=f"auxiliary {h} without {dropped}, {what} (leave-one-cohort)")


def microarray_effects() -> None:
    table = pd.read_csv(NC.ROOT / P_MA, sep="\t")
    for h in table["hypothesis"]:
        sel = f"hypothesis == '{h}'"
        for col, what in (("n", "n"), ("diff_sa_minus_comparator", "difference"), ("ci_low", "95% CI lower end"),
                          ("ci_high", "95% CI upper end")):
            tag = {"n": "n", "diff_sa_minus_comparator": "diff", "ci_low": "low", "ci_high": "high"}[col]
            add(f"mae_{h}_{tag}", P_MA, sel, col, decimals=0 if col == "n" else 3,
                cohort="microarray", method="HostMix-TOO",
                metric={"AH1": "host-attraction rate"}.get(h, "top-1 accuracy"),
                set_={"AH1": "at risk", "AH2": "evaluation", "AH3": "native truth"}.get(h, ""),
                meaning=f"microarray layer, {h} contrast, {what} (pre-specified descriptive)")


def build() -> list[dict]:
    a14()
    hypotheses()
    overall_pog()
    external()
    subcohorts()
    stage7_cohorts()
    ablation()
    stage10()
    simulation()
    tcga()
    development()
    variants()
    leave_one_cohort()
    microarray_effects()
    import facts_specific_A17 as SPEC  # noqa: E402
    SPEC.define(add, globals())
    for row in FACTS:
        row["status"] = NC.status_of(row["source_file"], row["selector"]) if not row["transform"].startswith("derived:") \
            else ""
    by_id = {row["fact_id"]: row for row in FACTS}
    for row in FACTS:
        if row["transform"].startswith("derived:"):
            names = [n for n in set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", row["transform"][8:])) if n in by_id]
            statuses = {by_id[n]["status"] for n in names if by_id[n]["status"]}
            if not names:
                row["status"] = "design"
            else:
                row["status"] = "post hoc" if not statuses or len(statuses) > 1 else statuses.pop()
            if not row["source_file"]:
                row["source_file"] = "derived"
    return FACTS


def main() -> None:
    facts = build()
    out = NC.CHECKS / "facts.tsv"
    pd.DataFrame(facts).to_csv(out, sep="\t", index=False)
    print(out, len(facts), "facts")


if __name__ == "__main__":
    main()
