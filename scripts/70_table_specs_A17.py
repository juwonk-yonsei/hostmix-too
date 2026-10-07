"""Write the table specs for Additional file 1 and MODEL_CARD (A17 §4.4).

Each spec maps a row label to source selectors by name, never by value.
Run once; the YAML files are the committed specs read by 69_check_numbers_A17.py.
"""
from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "checks" / "table_specs"

DEV = "results/stage2/met500_overall.tsv"
RANK = "results/stage2/candidate_rank.tsv"
PROXY = "results/stage2/proxy_correlation.tsv"
GTEX = "config/split_gtex.tsv"
VAR_TT = "results/stage5/tcga_test/variant_overall.tsv"
SIM = "results/stage5/sim_ext/full.tsv"
PHM = "results/stage5/posthoc/met500_overall.tsv"
PHP = "results/stage5/posthoc/pog_overall.tsv"
TRAIN_LOG = "results/stage5/variants/train_log.tsv"
OVERALL = "results/stage4/confirm/secondary_overall.tsv"
SITE = "results/stage4/confirm/secondary_by_site.tsv"
TC = "results/stage4/confirm/secondary_tc.tsv"
POWER3 = "results/stage3/power.json"
POWER6 = "results/stage6/power.tsv"
PAIRED = "results/stage2/met500_paired.tsv"
INPUTS = "config/external_inputs_A8.tsv"
PIPE = "results/stage9/ablation/pipeline_table.tsv"

S4_ROWS = {
    "Baseline": ("BASE", "Z"), "Native-organ masking": ("M1", "Z"), "Normal classes": ("NC", "Z"),
    "Linear deconvolution": ("LD", "Z"), "HostMix-TOO": ("SA", "Z"), "Gene removal 20%": ("IF20", "Z"),
    "Site-specific mixtures": ("SC", "Z"), "Site-specific + removal 20%": ("SC+IF20", "Z"),
    "Site-specific + removal 40%": ("SC+IF40", "Z"), "Gene sets": ("BASE", "K"),
    "Gene-set form of Native-organ masking": ("M1", "K"), "Gene sets + host correction": ("V0", "K"),
    "Gene-set form of normal classes": ("NC", "K"), "Gene-set form of Linear deconvolution": ("LD", "K"),
    "Gene sets + mixtures": ("SA", "K"), "Gene-set form of Gene removal 20%": ("IF20", "K"),
    "Gene-set form of Site-specific mixtures": ("SC", "K"),
    "Gene-set form of Site-specific + removal 20%": ("SC+IF20", "K"),
    "Gene-set form of Site-specific + removal 40%": ("SC+IF40", "K"),
}
REMOVAL = {"IF20": "if_q20", "SC+IF20": "if_q20", "SC+IF40": "if_q40"}

S5_ROWS = {
    "Baseline": "BASE-Z", "HostMix-TOO": "SA-Z", "Standardized, pure tumors": "SA-m0",
    "One mixture per tumor": "SA-m1", "Eight mixtures per tumor": "SA-m8",
    "Mixture lower bound 0.3": "SA-r30", "Mixture lower bound 0.5": "SA-r50",
    "Mixtures, C = 0.01": "SA-C0.01", "Mixtures, C = 0.1": "SA-C0.1", "3-tissue model": "SA-pool3",
    "22-tissue model": "SA-pool22", "Perceptron": "BASE-MLP", "Perceptron + mixtures": "SA-MLP",
}

S6_ROWS = {
    "Baseline": "BASE-Z", "Gene sets": "BASE-K", "HostMix-TOO": "SA-Z", "Site-specific mixtures": "SC-Z",
    "Linear deconvolution": "LD-Z", "Normal classes": "NC-Z", "Native-organ masking": "M1-Z",
    "Gene sets + mixtures": "SA-K", "Gene sets + host correction": "V0-K",
}

ABL_MODELS = {
    "Baseline": ("base", "BASE-Z"), "Standardized, pure (C = 0.03)": ("pure03", "PURE-Zs"),
    "Standardized, pure (C = 0.15)": ("pure15", "PURE-Zs-w"),
    "Unstandardized, mixtures (C = 0.1)": ("mix", "MIX-Z0"), "HostMix-TOO": ("sa", "SA-Z"),
    "Unstandardized mixture model": ("mix", "MIX-Z0"),
}
ABL_COHORT = {"MET500": "MET500", "POG570": "POG570", "Auxiliary RNA-seq": "aux_rnaseq"}

METRIC_COLS = ["n", "top1", "n_at_risk", "host_rate", "n_native_truth", "native_truth_top1"]
METRIC_IDS = ["n", "top1", "nrisk", "host", "nnative", "native"]


def num(src: str, sel: str, col: str, dec: int = 3, **extra) -> dict:
    return {"src": src, "sel": sel, "col": col, "dec": dec, **extra}


def metric_columns(start: int, src: str, sel: str, cols=METRIC_COLS, ids=METRIC_IDS) -> dict:
    out = {}
    for j, (col, cid) in enumerate(zip(cols, ids)):
        dec = 0 if col in ("n", "n_at_risk", "n_native_truth") else 3
        out[start + j] = {"id": cid, "t": "{v}", "f": {"v": num(src, sel, col, dec)}}
    return out


def s3() -> list[dict]:
    pool = ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Brain - Cortex", "Liver",
            "Lung", "Muscle - Skeletal", "Skin - Not Sun Exposed (Suprapubic)", "Spleen", "Whole Blood"]
    counts = {
        "id": "S3a", "doc": "af1", "caption": "Table S3.", "which": 0,
        "header": ["Host-pool tissue", "GTEx-ref samples"], "key_cols": [0], "row_id": "{t}",
        "maps": {0: {t: {"t": t} for t in pool}},
        "columns": {1: {"id": "n", "t": "{v}", "f": {"v": num(GTEX, "split == 'GTEx-ref' and tissue == '{t}'", "",
                                                             0, tf="count")}}},
    }
    cands = ["Spleen", "Whole Blood", "Cells - Ebv-Transformed Lymphocytes", "Small Intestine - Terminal Ileum"]
    sel = "hpa_query == '{q}' and gtex_tissue == '{g}'"
    corr = {
        "id": "S3b", "doc": "af1", "caption": "Table S3.", "which": 1,
        "header": ["Query", "Candidate tissue", "Genes compared", "GTEx-ref samples", "Spearman"],
        "key_cols": [0, 1], "row_id": "{q}/{g}",
        "maps": {0: {q: {"q": q} for q in ("lymph node", "bone marrow")}, 1: {g: {"g": g} for g in cands}},
        "columns": {
            2: {"id": "genes", "t": "{v}", "f": {"v": num(PROXY, sel, "n_genes", 0)}},
            3: {"id": "gtex_ref", "t": "{v}", "f": {"v": num(PROXY, sel, "n_gtex_ref", 0)}},
            4: {"id": "spearman", "t": "{v}", "f": {"v": num(PROXY, sel, "spearman", 0, form="raw")}},
        },
    }
    return [counts, corr]


def s4() -> dict:
    rows = {}
    for label, (m, r) in S4_ROWS.items():
        rows[label] = {"m": m, "r": r, "mid": f"{m}-{r}"}
    sel = "method == '{m}' and representation == '{r}'"
    col0_rows = {}
    for label, (m, r) in S4_ROWS.items():
        if m in REMOVAL:
            col0_rows[f"{m}-{r}"] = {"t": label.replace("20%", "{q}%").replace("40%", "{q}%"),
                                     "f": {"q": {"const": REMOVAL[m]}}}
    return {
        "id": "S4", "doc": "af1", "caption": "Table S4.",
        "header": ["Classifier", "Top-1", "Top-1 interval", "Host-attraction rate", "Host-rate interval",
                   "Native-truth top-1", "Candidate rank"],
        "key_cols": [0], "row_id": "{mid}", "maps": {0: rows},
        "columns": {
            0: {"id": "classifier", "t": "{*}", "rows": col0_rows},
            1: {"id": "top1", "t": "{v}", "f": {"v": num(DEV, sel, "top1")}},
            2: {"id": "top1_ci", "t": "{lo} to {hi}",
                "f": {"lo": num(DEV, sel, "top1_ci[0]"), "hi": num(DEV, sel, "top1_ci[1]")}},
            3: {"id": "host", "t": "{v}", "f": {"v": num(DEV, sel, "host_rate")}},
            4: {"id": "host_ci", "t": "{lo} to {hi}",
                "f": {"lo": num(DEV, sel, "host_rate_ci[0]"), "hi": num(DEV, sel, "host_rate_ci[1]")}},
            5: {"id": "native", "t": "{v}", "f": {"v": num(DEV, sel, "native_truth_top1")}},
            6: {"id": "rank", "t": "{*}", "rows": {
                mid: {"t": "{v}", "f": {"v": num(RANK, "method_id == '{mid}'", "__row__", 0,
                                                 meaning="candidate rank (row order of the ranking file)")}}
                for mid in ("SA-Z", "IF20-Z", "SC-Z", "SC+IF20-Z", "SC+IF40-Z", "SA-K", "IF20-K", "SC-K",
                            "SC+IF20-K", "SC+IF40-K")}},
        },
    }


def s5() -> dict:
    sel = "method == '{m}'"
    sim_sel = "method == '{m}' and rho == 0.6 and tissue in @{p}"
    col0 = {
        "SA-r30": {"t": "Mixture lower bound {v}", "f": {"v": {"const": "rho_low_r30"}}},
        "SA-r50": {"t": "Mixture lower bound {v}", "f": {"v": {"const": "rho_low_r50"}}},
        "SA-C0.01": {"t": "Mixtures, C = {v}", "f": {"v": num(TRAIN_LOG, "model == 'SA-C0.01'", "C", 2)}},
        "SA-C0.1": {"t": "Mixtures, C = {v}", "f": {"v": num(TRAIN_LOG, "model == 'SA-C0.1'", "C", 1)}},
        "SA-m8": {"t": "{v_w} mixtures per tumor", "f": {"v_w": {"fact": "var_m8_nmix"}}},
    }
    return {
        "id": "S5", "doc": "af1", "caption": "Table S5.",
        "header": ["Variant", "TCGA-test top-1", "In-pool host rate at rho 0.6", "Pool-out host rate at rho 0.6",
                   "MET500 top-1", "MET500 host rate", "POG570 top-1", "POG570 host rate",
                   "POG570 native-truth top-1"],
        "key_cols": [0], "row_id": "{m}", "maps": {0: {k: {"m": v} for k, v in S5_ROWS.items()}},
        "header_cells": {
            2: {"id": "h_in", "t": "In-pool host rate at rho {v}", "f": {"v": {"const": "rho_display"}}},
            3: {"id": "h_out", "t": "Pool-out host rate at rho {v}", "f": {"v": {"const": "rho_display"}}},
        },
        "columns": {
            0: {"id": "variant", "t": "{*}", "rows": col0},
            1: {"id": "tcga_top1", "t": "{v}", "f": {"v": num(VAR_TT, sel, "top1")}},
            2: {"id": "in_pool", "t": "{v}", "f": {"v": num(SIM, sim_sel.replace("{p}", "POOL"), "", 3,
                                                             tf="mean:host_rate")}},
            3: {"id": "pool_out", "t": "{v}", "f": {"v": num(SIM, sim_sel.replace("{p}", "OUT"), "", 3,
                                                              tf="mean:host_rate")}},
            4: {"id": "met_top1", "t": "{v}", "f": {"v": num(PHM, sel, "top1")}},
            5: {"id": "met_host", "t": "{v}", "f": {"v": num(PHM, sel, "host_rate")}},
            6: {"id": "pog_top1", "t": "{v}", "f": {"v": num(PHP, sel, "top1")}},
            7: {"id": "pog_host", "t": "{v}", "f": {"v": num(PHP, sel, "host_rate")}},
            8: {"id": "pog_native", "t": "{v}", "f": {"v": num(PHP, sel, "native_truth_top1")}},
        },
    }


def s6() -> list[dict]:
    model_map = {k: {"m": v} for k, v in S6_ROWS.items()}
    overall = {
        "id": "S6a", "doc": "af1", "caption": "Table S6.", "which": 0,
        "header": ["Classifier", "n", "Top-1", "n at risk", "Host-attraction rate", "n native truth",
                   "Native-truth top-1"],
        "key_cols": [0], "row_id": "{m}", "maps": {0: model_map},
        "columns": metric_columns(1, OVERALL, "method == '{m}'"),
    }
    groups = ["liver", "lung", "lymph_node", "remainder", "soft_tissue"]
    site = {
        "id": "S6b", "doc": "af1", "caption": "Table S6.", "which": 1,
        "header": ["Classifier", "Site group", "n", "Top-1", "n at risk", "Host-attraction rate", "n native truth",
                   "Native-truth top-1"],
        "key_cols": [0, 1], "row_id": "{m}/{g}", "maps": {0: model_map, 1: {g: {"g": g} for g in groups}},
        "columns": metric_columns(2, SITE, "method == '{m}' and site_group == '{g}'"),
    }
    tc = {
        "id": "S6c", "doc": "af1", "caption": "Table S6.", "which": 2,
        "header": ["Classifier", "Tumor-content tertile", "n", "Top-1", "n at risk", "Host-attraction rate"],
        "key_cols": [0, 1], "row_id": "{m}/{b}", "maps": {0: model_map, 1: {b: {"b": b} for b in ("T1", "T2", "T3")}},
        "columns": metric_columns(2, TC, "method == '{m}' and tc_bin == '{b}'", METRIC_COLS[:4], METRIC_IDS[:4]),
    }
    return [overall, site, tc]


PREREG_ROWS = {
    ("P1", 1): {"H3": {"t": "{*} by more than {m} percentage points", "f": {"m": {"const": "ni_margin_pts", "abs": True}}}},
    ("P1", 2): {
        "H1": {"t": "{*} Binomial(n_disc, {p}) variable is at least the observed favorable count. If n_disc = {z}, p = {o}",
               "f": {"p": {"const": "mcnemar_null_p"}, "z": {"const": "ndisc_zero"}, "o": {"const": "p_no_discordant"}}},
        "H3": {"t": "{*} is greater than {b}. Judged only if H2 is met", "f": {"b": {"const": "ni_bound"}}}},
    ("P2", 1): {"AH3": {"t": "{*} by more than {m} percentage points",
                        "f": {"m": {"const": "ni_margin_pts_a6", "abs": True}}}},
    ("P2", 2): {"AH3": {"t": "{*} ({n} replicates) > {b}. Judged only if AH2 is met",
                        "f": {"n": {"const": "bootstrap_n_a6"}, "b": {"const": "ni_bound_a6"}}}},
}


def prereg_tables() -> list[dict]:
    specs = []
    for sid, caption, which, rows in (("P1", "POG570 confirmatory hypotheses", 0, ["H1", "H2", "H3"]),
                                      ("P2", "Auxiliary confirmatory hypotheses", 0, ["AH1", "AH2", "AH3"]),
                                      ("P3", "Auxiliary confirmatory hypotheses", 1, ["AS1", "AS2", "AS3"])):
        specs.append({
            "id": sid, "doc": "af1", "caption": caption, "which": which, "key_cols": [0], "row_id": "{h}",
            "maps": {0: {h: {"h": h} for h in rows}},
            "columns": {1: {"id": "statement", "t": "{*}", "rows": PREREG_ROWS.get((sid, 1), {})},
                        2: {"id": "test", "t": "{*}", "rows": PREREG_ROWS.get((sid, 2), {})}},
        })
    return specs


def power_tables() -> list[dict]:
    rows0 = {"H1": {"h": "H1", "fav": "n_host_fixed", "aga": "n_host_new", "den": "dev_base_nrisk"},
             "H2": {"h": "H2", "fav": "n_method_only", "aga": "n_base_only", "den": "dev_base_n"}}
    sc = {"MET500 proportions": {"s": "full"}, "Half effect": {"s": "half"}}
    psel = "method == 'SA' and representation == 'Z'"
    pog = {
        "id": "S10a", "doc": "af1", "caption": "10. Design power", "which": 0,
        "header": ["Hypothesis", "Scenario", "n", "p_favor", "p_against", "α", "Simulations", "Power"],
        "key_cols": [0, 1], "row_id": "{h}_{s}", "maps": {0: rows0, 1: sc},
        "columns": {
            1: {"id": "scenario", "t": "{*}", "rows": {f"{h}_half": {"t": "{h_w} effect", "f": {"h_w": {"const": "half_factor"}}}
                                                      for h in ("H1", "H2")}},
            2: {"id": "n", "t": "{v}", "f": {"v": num(POWER3, "all", "{h}_{s}.n", 0)}},
            3: {"id": "p_favor", "t": "{v} ({a}/{b})",
                "f": {"v": num(POWER3, "all", "{h}_{s}.p_favor"), "a": num(PAIRED, psel, "{fav}", 0),
                      "b": {"fact": "{den}"}},
                "rows": {
                    "H1_half": {"f": {"a": {"fact": "pw_h1_half_favor"}}},
                    "H2_half": {"f": {"b": {"fact": "pw_h2_half_den"}}},
                }},
            4: {"id": "p_against", "t": "{v} ({a}/{b})",
                "f": {"v": num(POWER3, "all", "{h}_{s}.p_against"), "a": num(PAIRED, psel, "{aga}", 0),
                      "b": {"fact": "{den}"}}},
            5: {"id": "alpha", "t": "{v}", "f": {"v": num(POWER3, "all", "alpha", 2, form="num")}},
            6: {"id": "n_sim", "t": "{v}", "f": {"v": num(POWER3, "all", "{h}_{s}.n_sim", 0)}},
            7: {"id": "power", "t": "{v}", "f": {"v": num(POWER3, "all", "{h}_{s}.power")}},
        },
    }
    sel = "hypothesis == '{h}' and scenario == {s}"
    alpha_rows = {f"{h}_{k}": {"f": {"v": num(POWER6, sel, "alpha", 2)}} for h in ("AH1", "AH2") for k in "12"}
    aux = {
        "id": "S10b", "doc": "af1", "caption": "10. Design power", "which": 1,
        "header": ["hypothesis", "scenario", "n", "p_favor", "p_against", "alpha", "power"],
        "key_cols": [0, 1], "row_id": "{h}_{s}",
        "maps": {0: {h: {"h": h} for h in ("AH1", "AH2", "AS1", "AS2", "AS3")}, 1: {s: {"s": s} for s in ("1", "2")}},
        "columns": {
            1: {"id": "scenario", "t": "{v}", "f": {"v": num(POWER6, sel, "scenario", 0, meaning="scenario number")}},
            2: {"id": "n", "t": "{v}", "f": {"v": num(POWER6, sel, "n", 0)}},
            3: {"id": "p_favor", "t": "{v}", "f": {"v": num(POWER6, sel, "p_favor")}},
            4: {"id": "p_against", "t": "{v}", "f": {"v": num(POWER6, sel, "p_against")}},
            5: {"id": "alpha", "t": "{v}", "f": {"v": num(POWER6, sel, "alpha")}, "rows": alpha_rows},
            6: {"id": "power", "t": "{v}", "f": {"v": num(POWER6, sel, "power")}},
        },
    }
    return [pog, aux]


def inputs_table() -> dict:
    cohorts = {"MET500": "MET500", "POG570": "POG570", "IMvigor210": "blca_iatlas_imvigor210_2017",
               "Anders": "brca_iatlas_anders_2022", "DFCI melanoma": "mel_dfci_2019",
               "PRINCE": "paad_iatlas_prince_2022", "GSE50760": "GSE50760", "SU2C/PCF": "prad_su2c_2019"}
    return {
        "id": "S11", "doc": "af1", "caption": "11. Published classifiers", "key_cols": [0], "row_id": "{c}",
        "maps": {0: {k: {"c": v} for k, v in cohorts.items()}},
        "columns": {
            1: {"id": "values", "t": "{*}"}, 2: {"id": "scope", "t": "{*}"}, 5: {"id": "cup", "t": "{*}"},
            6: {"id": "evidence", "t": "`config/external_inputs_A8.tsv`: {v}",
                "f": {"v": num(INPUTS, "cohort == '{c}'", "__line__", 0, meaning="line of the cited input table")}},
        },
    }


def s8() -> dict:
    cols = {}
    for j, (cid, metric) in enumerate(zip(METRIC_IDS, METRIC_COLS)):
        cols[2 + j] = {"id": cid, "t": "{v}", "f": {"v": {"fact": "abl_{c}_{t}_" + cid}}}
    col1 = {}
    for label, (tag, _) in ABL_MODELS.items():
        if "C = " in label:
            col1[tag] = {"t": label.split("C = ")[0] + "C = {v})", "f": {"v": {"const": "c_" + tag}}}
    model_map = {k: {"t": t, "mm": mm} for k, (t, mm) in ABL_MODELS.items() if k != "Unstandardized mixture model"}
    return {
        "id": "S8", "doc": "af1", "caption": "Table S8.",
        "header": ["Cohort", "Model", "n", "Top-1", "n at risk", "Host-attraction rate", "n native-truth",
                   "Native-truth top-1"],
        "key_cols": [0, 1], "row_id": "{c}/{t}",
        "maps": {0: {k: {"c": v} for k, v in ABL_COHORT.items()}, 1: model_map},
        "columns": {1: {"id": "model", "t": "{*}", "rows": {f"{c}/{t}": v for c in ABL_COHORT.values()
                                                          for t, v in col1.items()}}, **cols},
    }


def s9() -> dict:
    model_map = {k: {"t": t, "mm": mm} for k, (t, mm) in ABL_MODELS.items() if k != "Unstandardized mixture model"}
    col0 = {}
    for label, (tag, _) in ABL_MODELS.items():
        if "C = " in label:
            col0[tag] = {"t": label.split("C = ")[0] + "C = {v})", "f": {"v": {"const": "c_" + tag}}}
    return {
        "id": "S9", "doc": "af1", "caption": "Table S9.", "key_cols": [0], "row_id": "{t}",
        "maps": {0: model_map},
        "header_cells": {
            1: {"id": "h_top1", "t": "Top-1 (n = {v})", "f": {"v": {"fact": "ma_layer_n_evaluation"}}},
            2: {"id": "h_host", "t": "Host-attraction rate (n = {v})", "f": {"v": {"fact": "ma_layer_n_at_risk"}}},
            3: {"id": "h_native", "t": "Native-truth top-1 (n = {v})", "f": {"v": {"fact": "ma_layer_n_native_truth"}}},
        },
        "columns": {
            0: {"id": "model", "t": "{*}", "rows": col0},
            1: {"id": "top1", "t": "{v}", "f": {"v": {"fact": "ma_layer_{t}_top1"}}},
            2: {"id": "host", "t": "{v}", "f": {"v": {"fact": "ma_layer_{t}_host_rate"}}},
            3: {"id": "native", "t": "{v}", "f": {"v": {"fact": "ma_layer_{t}_native_truth_top1"}}},
            4: {"id": "esophagus", "t": "{v}", "f": {"v": {"fact": "ma_layer_{t}_n_esophagus"}}},
            5: {"id": "random_loss", "t": "{v}", "f": {"v": {"fact": "mg_{t}_random_mean_mean"}}},
        },
    }


def s7() -> dict:
    return {
        "id": "S7", "doc": "af1", "caption": "Table S7.", "key_cols": [0], "row_id": "{c0}",
        "maps": {0: {"re:.": {}}},
        "columns": {0: {"id": "name", "t": "{*}", "rows": {
            "Gene removal 20%": {"t": "Gene removal {q}%", "f": {"q": {"const": "if_q20"}}}}}},
    }


def card() -> list[dict]:
    first = {
        "id": "CARD1", "doc": "card", "caption": "Features and classes", "key_cols": [0], "row_id": "{m}",
        "header": ["Model", "Scaler", "C", "How C was chosen", "Rows"],
        "maps": {0: {"BASE-Z": {"m": "BASE-Z"}, "SA-Z": {"m": "SA-Z"}}},
        "columns": {
            2: {"id": "C", "t": "{v}", "f": {"v": num(PIPE, "model == '{m}'", "C", 0, form="raw")}},
            3: {"id": "how", "t": "{*}", "rows": {"SA-Z": {
                "t": "{k}-fold grouped cross-validation, macro-F1, grid {a}, {b}, {c}, smallest C within {e} of the best",
                "f": {"k": {"const": "cv_folds_sa"}, "a": {"const": "c_grid_sa_1"}, "b": {"const": "c_grid_sa_2"},
                      "c": {"const": "c_grid_sa_3"}, "e": {"const": "c_tie_tol"}}}}},
            4: {"id": "rows", "t": "{v}", "f": {"v": num(PIPE, "model == '{m}'", "n_train_rows", 0)}},
        },
    }
    cohort_rows = {
        "re:^MET500": {"c": "MET500", "ct": "MET500", "kind": "rna"},
        "re:^POG570": {"c": "POG570", "ct": "POG570", "kind": "rna"},
        "re:^Auxiliary RNA-seq": {"c": "aux_rnaseq", "ct": "Auxiliary RNA-seq", "kind": "rna"},
        "re:^Microarray": {"c": "layer", "ct": "Microarray", "kind": "ma"},
    }
    model_map = {k: {"t": t} for k, (t, _) in ABL_MODELS.items()
                 if k in ("Baseline", "Unstandardized mixture model", "HostMix-TOO")}
    rows0 = {}
    for c in ("MET500", "POG570", "aux_rnaseq"):
        rows0[f"{c}/base"] = {"t": "{*} (n = {n}, {r}, {v})",
                              "f": {"n": {"fact": f"abl_{c}_base_n"}, "r": {"fact": f"abl_{c}_base_nrisk"},
                                    "v": {"fact": f"abl_{c}_base_nnative"}}}
    rows0["layer/base"] = {"t": "{*} (n = {n}, {r}, {v})",
                           "f": {"n": {"fact": "ma_layer_n_evaluation"}, "r": {"fact": "ma_layer_n_at_risk"},
                                 "v": {"fact": "ma_layer_n_native_truth"}}}
    metric = {"top1": ("top1", "top1"), "host": ("host", "host_rate"), "native": ("native", "native_truth_top1")}
    cols = {}
    for j, (cid, (abl, ma)) in enumerate(metric.items()):
        rows = {f"layer/{t}": {"f": {"v": {"fact": f"ma_layer_{t}_{ma}"}}} for t in ("base", "mix", "sa")}
        cols[2 + j] = {"id": cid, "t": "{v}", "f": {"v": {"fact": "abl_{c}_{t}_" + abl}}, "rows": rows}
    second = {
        "id": "CARD2", "doc": "card", "caption": "Post hoc model", "key_cols": [0, 1], "carry": [0],
        "row_id": "{c}/{t}",
        "header": ["Cohort", "Model", "Top-1", "Host-attraction rate", "Native-truth top-1"],
        "maps": {0: cohort_rows, 1: model_map},
        "columns": {0: {"id": "cohort", "t": "{*}", "rows": rows0}, **cols},
    }
    return [first, second]


FS = "manuscript/bmc/figure_source"


def af2() -> list[dict]:
    def c(sheet, source, header_row=2, **params):
        return {"id": f"AF2_{sheet}", "doc": "af2", "sheet": sheet, "header_row": header_row, "builder": "copy",
                "params": {"source": source, **params}, "origin": "figure_source" if source.startswith(FS) else "results"}
    out = [
        c("Fig2_simulation", "results/stage5/sim_ext/full.tsv"),
        c("Fig2ab_means", f"{FS}/Fig2_means.tsv", rename={"host_pull_rate": "host_attraction_rate"},
          maps={"pool": {"in": "host pool", "out": "outside the host pool"}},
          add_map={"display": ["method", "METHOD_LABEL"]},
          sort={"by": ["pool", "method", "rho"], "ascending": [True, True, False]},
          columns=["pool", "rho", "method", "display", "n_tissues", "host_attraction_rate"]),
        c("Fig2c_leave_one_host_out", f"{FS}/Fig2_loho.tsv", rename={"hostmix": "hostmix_too"}),
        c("Fig3_Table2_confirmatory", "results/stage4/confirm/primary.tsv"),
        c("Fig3_auxiliary_primary", "results/stage7/confirm/tables/hypothesis_primary.tsv"),
        c("Fig3_auxiliary_secondary", "results/stage7/confirm/tables/hypothesis_secondary.tsv"),
        c("Fig4_ablation", f"{FS}/Fig4_ablation.tsv",
          columns=["panel", "cohort", "display", "method", "label", "metric", "value", "n", "source"]),
        c("S8_shares", "results/stage9/ablation/shares.tsv"),
        c("Fig5_Table3_published", "results/stage8/external/tables/metrics.tsv"),
        c("Fig5_paired", "results/stage8/external/tables/paired.tsv"),
        {"id": "AF2_Fig6a_cohorts", "doc": "af2", "sheet": "Fig6a_cohorts", "header_row": 2, "builder": "fig6a",
         "origin": "figure_source",
         "params": {"source": f"{FS}/Fig6a_cohorts.tsv", "files": {
             "stage8 metrics": "results/stage8/external/tables/metrics.tsv",
             "subcohort_metrics": "results/stage9/subcohort_metrics.tsv",
             "stage7 metrics": "results/stage7/confirm/tables/metrics.tsv",
             "post hoc aggregate of confirmatory predictions": "results/stage7/confirm/predictions.parquet"}}},
        c("Fig6b_POG570_strata", "results/stage4/confirm/secondary_by_site.tsv"),
        c("Fig6b_tumor_content", "results/stage4/confirm/secondary_tc.tsv"),
        c("Fig6c_aux_top1", f"{FS}/Fig6c_aux_top1.tsv"),
        c("Fig6d_removed_errors", f"{FS}/Fig6d_removed.tsv"),
        c("S1_heatmap", f"{FS}/S1_heatmap.tsv"),
        c("S2_sets", f"{FS}/S2_sets.tsv"),
        {"id": "AF2_S3_esophagus_masking", "doc": "af2", "sheet": "S3_esophagus_masking", "header_row": 2,
         "builder": "s3_masking", "origin": "results (values at 10 significant digits)",
         "params": {"source": "results/stage8/diagnostics/esophagus_absorption.tsv"}},
        c("S4_contributions", f"{FS}/S4_contributions.tsv",
          columns=["display", "method", "n_samples", "side", "order", "rank", "gene", "mean_contribution", "input"],
          keep_ids=True),
        c("S5_histology", f"{FS}/S5_histology.tsv"),
        c("S5_reconciliation", f"{FS}/S5_reconciliation.tsv"),
        {"id": "AF2_S6_meta_analysis", "doc": "af2", "sheet": "S6_meta_analysis", "header_row": 2, "builder": "s6_meta",
         "origin": "results and figure_source (ci_low, ci_high of cohort rows)",
         "params": {"meta": "results/stage7/confirm/tables/meta_analysis.tsv", "drawn": f"{FS}/S6_meta.tsv"}},
        c("Microarray_layer", "results/stage7/confirm/tables/microarray_effects.tsv"),
        c("Sensitivity", "results/stage7/confirm/tables/sensitivity.tsv"),
        c("Leave_one_cohort", "results/stage8/diagnostics/leave_one_cohort.tsv"),
        c("Conformal", "results/stage4/confirm/conformal.tsv"),
        c("Power", "results/stage6/power.tsv"),
        {"id": "AF2_S8_ablation", "doc": "af2", "sheet": "S8_ablation", "header_row": 3, "builder": "s8_display",
         "origin": "results (three-decimal display strings)", "params": {"source": "results/stage9/ablation/metrics.tsv"}},
        {"id": "AF2_S9_ablation_microarray", "doc": "af2", "sheet": "S9_ablation_microarray", "header_row": 2,
         "builder": "s9_microarray", "origin": "results",
         "params": {"metrics": "results/stage10/ablation_microarray/metrics.tsv",
                    "restriction": "results/stage10/ablation_microarray/restriction.tsv",
                    "shares": "results/stage10/ablation_microarray/shares.tsv",
                    "category": "results/stage10/ablation_microarray/category.tsv"}},
    ]
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    specs = [*s3(), s4(), s5(), *s6(), *prereg_tables(), *power_tables(), inputs_table(), s8(), s9(), s7(), *card()]
    for spec in af2():
        path = OUT / f"{spec['id']}.yaml"
        path.write_text(yaml.safe_dump(spec, sort_keys=False, allow_unicode=True, width=200))
    for spec in specs:
        path = OUT / f"{spec['id']}.yaml"
        path.write_text(yaml.safe_dump(spec, sort_keys=False, allow_unicode=True, width=200))
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
