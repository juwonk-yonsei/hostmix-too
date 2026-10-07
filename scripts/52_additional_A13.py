"""Build Additional files 1–3 for the stage-13 manuscript.

No new model fits. Tables are read from stored result files.
Internal model names appear only in Table S7 and in Additional files 2 and 3.
"""
from __future__ import annotations

import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pandas as pd
import yaml
from openpyxl import Workbook, load_workbook

ROOT = Path(__file__).resolve().parents[1]
BMC = ROOT / "manuscript" / "bmc"

DISPLAY = {
    "BASE-Z": "Baseline",
    "SA-Z": "HostMix-TOO",
    "SC-Z": "Site-specific mixtures",
    "NC-Z": "Normal classes",
    "LD-Z": "Linear deconvolution",
    "M1-Z": "Native-organ masking",
    "IF20-Z": "Gene removal 20%",
    "SC+IF20-Z": "Site-specific + removal 20%",
    "SC+IF40-Z": "Site-specific + removal 40%",
    "BASE-K": "Gene sets",
    "SA-K": "Gene sets + mixtures",
    "V0-K": "Gene sets + host correction",
    "SA-G": "Gated model",
    "SA-pool22": "22-tissue model",
    "SA-pool3": "3-tissue model",
    "BASE-MLP": "Perceptron",
    "SA-MLP": "Perceptron + mixtures",
    "SA-m0": "Scaled, pure tumors",
    "SA-m1": "One mixture per tumor",
    "SA-m8": "Eight mixtures per tumor",
    "SA-r30": "Mixture lower bound 0.3",
    "SA-r50": "Mixture lower bound 0.5",
    "SA-C0.01": "Mixtures, C = 0.01",
    "SA-C0.1": "Mixtures, C = 0.1",
    "PURE-Zs": "Scaled, pure (C = 0.03)",
    "PURE-Zs-w": "Scaled, pure (C = 0.15)",
    "MIX-Z0": "Unscaled, mixtures (C = 0.1)",
}


def show(method: str) -> str:
    if method in DISPLAY:
        return DISPLAY[method]
    if method.startswith("SA-LOHO-"):
        return "Leave-one-host-out (" + method[len("SA-LOHO-"):] + ")"
    if method.endswith("-K"):
        base = DISPLAY.get(method[:-2] + "-Z", method[:-2])
        return "Gene-set form of " + base
    return method


def r3(value) -> str:
    number = Decimal(str(value))
    return format(number.quantize(Decimal("0.001"), rounding=ROUND_HALF_UP), "f")


def md_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in frame.itertuples(index=False):
        cells = ["" if pd.isna(cell) else str(cell) for cell in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def host_pool() -> list[str]:
    plan = yaml.safe_load((ROOT / "config/analysis_plan_A2.yaml").read_text())
    return list(plan["host_pool_sorted"])


def section_tables() -> dict[str, str]:
    organs = read_tsv(ROOT / "config/mappings/project_to_organ.tsv")
    organs.columns = ["Tumor label", "Organ group"]
    sites = read_tsv(ROOT / "config/mappings/site_native.tsv")
    pool = set(host_pool())
    sites["In host pool"] = sites["gtex_tissues"].fillna("").map(
        lambda text: "yes" if any(part in pool for part in str(text).split("|") if part) else "no"
    )
    sites = sites.rename(columns={
        "standard_site": "Standard site",
        "gtex_tissues": "Host tissue or proxy",
        "native_organs": "Native organ set",
    })[["Standard site", "Host tissue or proxy", "Native organ set", "In host pool"]]

    split = read_tsv(ROOT / "config/split_gtex.tsv")
    ref = split[split["split"] == "GTEx-ref"]
    counts = (
        ref[ref["tissue"].isin(pool)]
        .groupby("tissue")
        .size()
        .reindex(host_pool())
        .rename("GTEx-ref samples")
        .reset_index()
        .rename(columns={"tissue": "Host-pool tissue"})
    )
    proxy = read_tsv(ROOT / "results/stage2/proxy_correlation.tsv")
    proxy = proxy.rename(columns={
        "hpa_query": "Query",
        "gtex_tissue": "Candidate tissue",
        "n_genes": "Genes compared",
        "n_gtex_ref": "GTEx-ref samples",
        "spearman": "Spearman",
    })
    proxy["Spearman"] = proxy["Spearman"].map(lambda value: format(Decimal(str(value)), "f"))
    proxy = proxy.drop(columns=["hpa_tissue"])

    rank = read_tsv(ROOT / "results/stage2/candidate_rank.tsv")
    rank["rank"] = range(1, len(rank) + 1)
    overall = read_tsv(ROOT / "results/stage2/met500_overall.tsv")
    overall["method_id"] = overall["method"] + "-" + overall["representation"]
    merged = overall.merge(rank[["method_id", "rank"]], on="method_id", how="left")
    s4 = pd.DataFrame({
        "Classifier": merged["method_id"].map(show),
        "Top-1": merged["top1"].map(r3),
        "Top-1 interval": merged["top1_ci"],
        "Host-attraction rate": merged["host_rate"].map(r3),
        "Host-rate interval": merged["host_rate_ci"],
        "Native-truth top-1": merged["native_truth_top1"].map(r3),
        "Candidate rank": merged["rank"].map(lambda value: "" if pd.isna(value) else str(int(value))),
    })

    s6_overall = read_tsv(ROOT / "results/stage4/confirm/secondary_overall.tsv")
    s6_site = read_tsv(ROOT / "results/stage4/confirm/secondary_by_site.tsv")
    s6_tc = read_tsv(ROOT / "results/stage4/confirm/secondary_tc.tsv")
    for frame in (s6_overall, s6_site, s6_tc):
        frame["Classifier"] = frame["method"].map(show)
        for column in ("top1", "host_rate", "native_truth_top1"):
            frame[column] = frame[column].map(r3)

    overall_table = md_table(s6_overall[["Classifier", "n", "top1", "n_at_risk", "host_rate", "n_native_truth", "native_truth_top1"]].rename(columns={
        "top1": "Top-1", "n_at_risk": "n at risk", "host_rate": "Host-attraction rate",
        "n_native_truth": "n native truth", "native_truth_top1": "Native-truth top-1",
    }))
    site_table = md_table(s6_site[["Classifier", "site_group", "n", "top1", "n_at_risk", "host_rate", "n_native_truth", "native_truth_top1"]].rename(columns={
        "site_group": "Site group", "top1": "Top-1", "n_at_risk": "n at risk",
        "host_rate": "Host-attraction rate", "n_native_truth": "n native truth",
        "native_truth_top1": "Native-truth top-1",
    }))
    tc_table = md_table(s6_tc[["Classifier", "tc_bin", "n", "top1", "n_at_risk", "host_rate"]].rename(columns={
        "tc_bin": "Tumor-content tertile", "top1": "Top-1", "n_at_risk": "n at risk",
        "host_rate": "Host-attraction rate",
    }))
    return {
        "s1": "### Table S1. Tumor labels and organ groups\n\nThirty-two tumor labels are summed to 26 organ groups after prediction. Lung combines lung adenocarcinoma and lung squamous cell carcinoma. Colorectal combines colon and rectal adenocarcinoma.\n\n" + md_table(organs),
        "s2": "### Table S2. Standard biopsy sites and native organ sets\n\nA biopsy is at risk when its true organ is outside the native set of its standard site. A native-truth biopsy is one whose true organ is inside that set. Lymph node uses spleen as a proxy, and bone marrow uses whole blood. Tissues listed here that are not in the ten-tissue host pool are outside the pool.\n\n" + md_table(sites),
        "s3": "### Table S3. Host pool and proxy correlations\n\nGTEx-ref sample counts for the ten host-pool tissues, then the eight Spearman correlations used to choose proxies for lymph node and bone marrow. The selected proxies are spleen for lymph node and whole blood for bone marrow.\n\n" + md_table(counts) + "\n\n" + md_table(proxy),
        "s4": "### Table S4. MET500 development comparison\n\nMET500 was the development cohort. Intervals are the bootstrap intervals from the development table. Candidate rank is the development ordering among the ten candidates; a blank rank means the classifier was not in that candidate list. This comparison is exploratory.\n\n" + md_table(s4),
        "s6": "### Table S6. POG570 alternatives and pre-specified strata\n\nThese outputs were named before the POG570 labels were unlocked. They are pre-specified descriptive results, not hypothesis tests.\n\nOverall.\n\n" + overall_table + "\n\nBy biopsy-site group.\n\n" + site_table + "\n\nBy tumor-content tertile.\n\n" + tc_table,
    }


def section_s5() -> str:
    pool = set(host_pool())
    sim = read_tsv(ROOT / "results/stage5/sim_ext/full.tsv")
    sim = sim[sim["rho"] == 0.6].copy()
    sim["pool"] = sim["tissue"].isin(pool)
    means = sim.groupby(["method", "pool"])["host_rate"].mean().unstack()
    tcga = read_tsv(ROOT / "results/stage5/tcga_test/variant_overall.tsv").set_index("method")
    met = read_tsv(ROOT / "results/stage5/posthoc/met500_overall.tsv").set_index("method")
    pog = read_tsv(ROOT / "results/stage5/posthoc/pog_overall.tsv").set_index("method")
    order = [
        "BASE-Z", "SA-Z", "SA-m0", "SA-m1", "SA-m8", "SA-r30", "SA-r50",
        "SA-C0.01", "SA-C0.1", "SA-pool3", "SA-pool22", "BASE-MLP", "SA-MLP",
    ]
    rows = []
    for method in order:
        rows.append({
            "Variant": show(method),
            "TCGA-test top-1": r3(tcga.loc[method, "top1"]),
            "In-pool host rate at rho 0.6": r3(means.loc[method, True]),
            "Pool-out host rate at rho 0.6": r3(means.loc[method, False]),
            "MET500 top-1": r3(met.loc[method, "top1"]),
            "MET500 host rate": r3(met.loc[method, "host_rate"]),
            "POG570 top-1": r3(pog.loc[method, "top1"]),
            "POG570 host rate": r3(pog.loc[method, "host_rate"]),
            "POG570 native-truth top-1": r3(pog.loc[method, "native_truth_top1"]),
        })
    loho = [method for method in pog.index if str(method).startswith("SA-LOHO-")]
    liver = pog.loc["SA-LOHO-Liver", "host_rate"]
    text = (
        "### Table S5. Exploratory variants\n\n"
        "These variants were fit between the two confirmations. Simulation rates are unweighted means across tissues at tumor RNA fraction 0.6. "
        f"The leave-one-host-out model that omitted liver had a POG570 host-attraction rate of {r3(liver)}. "
        f"Across the {len(loho)} left-out tissues, the POG570 host-attraction rate of the matching leave-one-host-out model ranged from "
        f"{r3(pog.loc[loho, 'host_rate'].min())} to {r3(pog.loc[loho, 'host_rate'].max())}. "
        "These rows are exploratory.\n\n"
    )
    return text + md_table(pd.DataFrame(rows))


def section_s7() -> str:
    rows = [
        ("Baseline", "BASE-Z", "scripts/08_train.py", "results/stage2/models/BASE_Z.joblib"),
        ("HostMix-TOO", "SA-Z", "scripts/08_train.py", "results/stage2/models/SA_Z.joblib"),
        ("Normal classes", "NC-Z", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Site-specific mixtures", "SC-Z", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Linear deconvolution", "LD-Z", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Native-organ masking", "M1-Z", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Gene removal 20%", "IF20-Z", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Gene sets", "BASE-K", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Gene sets + mixtures", "SA-K", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Gene sets + host correction", "V0-K", "scripts/08_train.py", "results/stage2/met500_overall.tsv"),
        ("Gated model", "SA-G", "scripts/24_aux_confirm.py", "results/stage7/confirm/tables/metrics.tsv"),
        ("22-tissue model", "SA-pool22", "scripts/14_stage5_train.py", "results/stage5/posthoc/pog_overall.tsv"),
        ("3-tissue model", "SA-pool3", "scripts/14_stage5_train.py", "results/stage5/posthoc/pog_overall.tsv"),
        ("POG570 hypotheses", "H1 H2 H3", "scripts/12_confirm.py", "results/stage4/confirm/primary.tsv"),
        ("Auxiliary hypotheses", "AH AS", "scripts/24_aux_confirm.py", "results/stage7/confirm/tables/hypothesis_primary.tsv"),
        ("Published classifiers", "SCOPE CUP-AI-Dx", "scripts/29_cup_wrapper_A8.py", "results/stage8/external/tables/metrics.tsv"),
        ("Ablation", "PURE-Zs PURE-Zs-w MIX-Z0", "scripts/45_ablation_A12.py", "results/stage9/ablation/metrics.tsv"),
        ("Figure source", "display names", "scripts/50_figures_A13.py", "manuscript/bmc/figure_source/"),
    ]
    frame = pd.DataFrame(rows, columns=["Manuscript name", "Internal name", "Script", "Output"])
    return "### Table S7. Correspondence between manuscript names and project files\n\nThis is the only table in Additional file 1 that uses internal names and file paths.\n\n" + md_table(frame)


def legends() -> str:
    return """## 13. Supplementary figures

**Figure S1. Accuracy and host attraction of the locked classifiers**

Each cell is top-1 organ accuracy or the host-attraction rate. Columns are MET500, POG570, the six confirmatory RNA-seq cohorts, and the microarray cohorts. MET500 cells come from the development table and POG570 cells from the confirmation table. The auxiliary column is a post hoc aggregate of the confirmatory predictions on the six RNA-seq cohorts, using the evaluation and at-risk sets of Table 3. Microarray cells are the combined microarray row (n = 556), which is larger than the one-patient descriptive layer (n = 536) reported in the text. Printed labels use two decimals; the source table keeps the full values. A cell that was never computed is marked with a dash. The figure is descriptive.

**Figure S2. Accuracy and host attraction in three analysis sets**

Top-1 accuracy and the host-attraction rate for the baseline, HostMix-TOO, linear deconvolution and normal classes in the native-truth set, the at-risk set and the pool-out at-risk set. MET500, POG570 and the auxiliary RNA-seq cohorts are shown separately. This decomposition was made after both confirmations and is post hoc.

**Figure S3. Esophagus calls and the effect of missing genes**

**a** Fraction of non-esophagus truths called esophagus, by cohort group and classifier. **b** HostMix-TOO predictions for the 183 GSE41258 native-truth colorectal tumors; 102 were called esophagus. **c** Change in top-1 accuracy after restricting profiles to the genes of a microarray platform, with the mean and range of ten random removals of the same number of classifier genes. The missing-gene simulation is post hoc. Models were not retrained.

**Figure S4. Genes contributing to esophagus calls in IMvigor210**

Mean contribution of each gene to the difference between the esophagus and bladder logits, for IMvigor210 bladder tumors in the native-truth set that were called esophagus. The contribution is the coefficient difference times the rank-normal score. HostMix-TOO was trained on standardized scores; this figure does not use that standardized input. The five largest baseline contributions are labeled. This analysis is post hoc and does not identify a mechanism.

**Figure S5. Assignment of esophagus calls to adenocarcinoma or squamous centroids**

Each row is one cohort, in the same order in every panel. The number at the end of a bar is that panel's number of esophagus calls, including zero when the classifier made none. Of 147 HostMix-TOO esophagus calls in GSE41258, 145 were closer to the adenocarcinoma centroid and 2 to the squamous centroid. The same count for HostMix-TOO in IMvigor210 is 105, not 147. This assignment is post hoc.

**Figure S6. Cohort-level differences in host-attraction rate**

Each bar is the difference, HostMix-TOO minus baseline, in one cohort. The interval is a Wald interval constructed from the variance recorded for that cohort difference. The dashed line is the DerSimonian–Laird estimate already computed across seven cohorts, including the descriptive cohort GSE209998, which is labeled as descriptive. FHCRC is absent because its at-risk count was below 10. This summary is post hoc and is not a replacement for the preregistered tests.
"""


def screening_counts() -> str:
    log = read_tsv(ROOT / "results/stage5/aux/screening_log.tsv")
    decision = log["decision"].astype(str)
    screened = int(decision.str.contains("선별").sum())
    excluded = int(decision.str.contains("제외").sum())
    return (
        f"The search log records {len(log)} expression studies. "
        f"{screened} passed the expression and site-annotation screen and {excluded} did not. "
        "A cohort entered the confirmatory or descriptive layer only when, after one sample per patient and after duplicate exclusion, its at-risk set contained at least 10 biopsies. "
        "The studies that met that rule are the cohorts in Table 1."
    )


def duplicate_counts() -> str:
    labels = read_tsv(ROOT / "config/aux_eval_labels.tsv")
    excluded = labels[labels["excluded"].astype(str).str.lower().isin(["true", "1"])]
    counts = excluded.groupby("cohort").size()
    parts = [f"{cohort}: {int(count)}" for cohort, count in counts.items()]
    return (
        "A sample was excluded when its Spearman correlation of log2(TPM + 1) with any MET500 or POG570 sample was high enough to suggest a duplicate, using the threshold fixed before the auxiliary plan. "
        "Excluded samples by cohort: " + "; ".join(parts) + ". "
        "Cohorts that do not appear in that list had no such exclusion."
    )


def write_additional_1() -> None:
    prereg_a3 = (ROOT / "config/prereg_A3_en.md").read_text().strip()
    prereg_a6 = (ROOT / "config/prereg_A6_en.md").read_text().strip()
    power = read_tsv(ROOT / "results/stage6/power.tsv")
    power["power"] = power["power"].map(r3)
    planted = read_tsv(ROOT / "results/stage7/planted/compare.tsv")
    mismatches = int((~planted["match"].astype(bool)).sum())
    tables = section_tables()
    text = f"""# Additional file 1. Supplementary methods, tables and figures

## 1. Data processing

The Toil matrix downloaded from the UCSC Xena hub is on the log2(TPM + 0.001) scale. Linear TPM was obtained by max(2^x − 0.001, 0), where x is the downloaded value. A negative downloaded value is the log representation of a TPM below 0.001, and the expression is then zero on the linear scale.

Gene identifiers were mapped to symbols with GENCODE version 23. The gene universe is the 17552 HGNC protein-coding symbols present in the Toil matrix, the MET500 matrix and the POG570 matrix. Previous HGNC symbols were updated to the current symbol. Aliases were not used.

Within each cohort, values were converted to a linear scale summing to 10^6 per sample. FPKM values were rescaled. Counts were divided by the GENCODE version 23 exon-union length and then rescaled. Log-scale values were exponentiated. PRINCE provides log2-transformed normalized values, which were converted by 2^v − 1 without a length correction.

Microarray probes matching more than one gene were removed. For each remaining gene the probe with the highest mean signal was kept. GSE71729 and GSE74685 are two-channel arrays. The sample channel was used and the reference channel was not.

## 2. Cohort search and selection

Public metastatic expression studies were sought in cBioPortal molecular profiles of type MRNA_EXPRESSION, excluding z-score profiles, and in GEO series with a biopsy-site annotation. {screening_counts()}

## 3. Duplicate screening

{duplicate_counts()}

## 4. Labels and sites

Table S1 gives the 32 tumor labels and the 26 organ groups. Table S2 gives the standard site, the host tissue or proxy, and the native organ set.

Fifty-nine MET500 specimens were excluded from the evaluation set because the annotated diagnosis did not map to an organ group. Fifty-eight POG570 samples were excluded for the same reason, leaving 512. In IMvigor210, 67 samples with biopsy site Kidney and sample type Metastasis are outside the bladder native-truth set; they remain in the evaluation set when the truth organ is bladder and the site is a standard site.

{tables["s1"]}

{tables["s2"]}

## 5. Host pool and proxies

The host pool is ten GTEx tissues. Lymph node and bone marrow have no GTEx tissue of the same name. The proxy was the candidate tissue with the highest Spearman correlation of mean expression against the Human Protein Atlas tissue query. Table S3 lists the pool counts and the eight correlations.

{tables["s3"]}

## 6. Mixture generation and seeds

Each HostMix-TOO training tumor contributes one unchanged profile and four mixtures. The tumor fraction of RNA is drawn from the uniform distribution on [0.15, 1.0). The host is a uniform draw of one pool tissue and then a uniform GTEx-ref sample of that tissue. The mixture seed recorded for HostMix-TOO is 890338470. Bootstrap intervals use 2000 patient resamples and seed 20261001. The random gene-removal control uses seed 20261023. Cross-validation uses seed 20261001.

## 7. Development comparison and exploratory variants

Table S4 is the MET500 comparison used during development. Table S5 is the set of variants fit between the two confirmations. Neither table is a hypothesis test.

{tables["s4"]}

{section_s5()}

## 8. POG570 alternatives and strata

{tables["s6"]}

## 9. Registered hypotheses

The following English text is a translation. The Korean plan files are authoritative.

{prereg_a3}

{prereg_a6}

## 10. Design power, synthetic labels, planted predictions and the independent checker

Design power for the auxiliary family was computed before unlock. The table below is that calculation. It is not a result of the confirmatory data.

{md_table(power)}

Before unlock, the confirmatory code was run on synthetic labels and on planted predictions with known expected results. An independent script that does not import the analysis code recomputed the auxiliary test statistics. That comparison has {len(planted)} rows and {mismatches} mismatches. No separate independent checker was recorded for the POG570 confirmation.

## 11. Published classifiers

SCOPE and CUP-AI-Dx were run from the authors' public code and trained models. Inputs were the locked rank-normal matrices, restricted to each tool's gene list. SCOPE calls of a normal-tissue class, marked by the authors' normal-sample suffix, were counted as incorrect in the primary comparison. A sensitivity count that mapped those calls to the corresponding organ is reported in the main text.

The CUP-AI-Dx example in the authors' metastatic metadata file was reproduced before the study cohorts were scored: 363 of 394 examples matched the authors' recorded call. The possible overlap between POG570 and the SCOPE development samples could not be checked, because the SCOPE development identifiers were not available to us. MET500 was the development cohort for HostMix-TOO and is not a confirmation.

## 12. Post hoc analyses

The ablation plan was committed before the three extra logistic models were fit. For each cohort it defined the gap between the baseline and HostMix-TOO host-attraction rates and the share of that gap attributable to the mixtures. The reading rule was that the mixtures would be considered the main source if all nine shares were at least 0.5. That rule was met. The ablation is post hoc relative to the preregistered hypotheses.

The missing-gene simulation restricted POG570 and MET500 profiles to the genes measured on each microarray platform and compared the change in top-1 accuracy with ten random removals of the same number of classifier genes. Models were not retrained.

The gene contribution in Figure S4 is the difference of the esophagus and bladder coefficients multiplied by the rank-normal score, averaged over the IMvigor210 native-truth tumors called esophagus. It does not multiply by the standardized value that HostMix-TOO used as its input.

Esophagus calls were assigned to the adenocarcinoma or squamous centroid of the TCGA esophageal carcinoma class by distance in the rank-normal space. Figure S5 uses one cohort order in every panel.

{legends()}

## 14. Repository correspondence

{section_s7()}
"""
    if re.search(r"[가-힣]", text):
        raise SystemExit("Hangul remains in Additional file 1")
    dest = BMC / "Additional_file_1.md"
    dest.write_text(text)
    print(dest, dest.stat().st_size)


def drop_ids(frame: pd.DataFrame) -> pd.DataFrame:
    banned = [column for column in frame.columns if re.search(r"sample|patient|specimen", column, re.I)]
    return frame.drop(columns=banned)


def translate_cell(value):
    if not isinstance(value, str):
        return value
    return {
        "계산": "computed",
        "마이크로어레이": "microarray",
        "사전 지정 서술, 가설 검정 아님": "pre-specified descriptive, not a hypothesis test",
        "사후": "post hoc",
        "확인 검정에 없음": "not in the confirmatory test",
    }.get(value, value)


def add_sheet(book: Workbook, name: str, description: str, frame: pd.DataFrame) -> None:
    frame = drop_ids(frame)
    sheet = book.create_sheet(name[:31])
    sheet.append([description])
    sheet.append(list(frame.columns))
    for row in frame.itertuples(index=False):
        values = []
        for value in row:
            if pd.isna(value):
                values.append(None)
            elif hasattr(value, "item"):
                values.append(translate_cell(value.item()))
            else:
                values.append(translate_cell(value))
        sheet.append(values)


def write_additional_2() -> None:
    book = Workbook()
    readme = book.active
    readme.title = "README"
    sheets = [
        ("Fig2_simulation", "Exploratory simulation. Tissue by rho by model. Not a hypothesis test.", ROOT / "results/stage5/sim_ext/full.tsv"),
        ("Fig3_Table2_confirmatory", "Preregistered POG570 contrasts.", ROOT / "results/stage4/confirm/primary.tsv"),
        ("Fig3_auxiliary_primary", "Preregistered auxiliary primary contrasts.", ROOT / "results/stage7/confirm/tables/hypothesis_primary.tsv"),
        ("Fig3_auxiliary_secondary", "Preregistered auxiliary secondary contrasts. Not supported.", ROOT / "results/stage7/confirm/tables/hypothesis_secondary.tsv"),
        ("Fig4_ablation", "Post hoc ablation metrics.", ROOT / "results/stage9/ablation/metrics.tsv"),
        ("Fig4_shares", "Post hoc mixture shares.", ROOT / "results/stage9/ablation/shares.tsv"),
        ("Fig5_Table3_published", "Pre-specified descriptive published-classifier metrics.", ROOT / "results/stage8/external/tables/metrics.tsv"),
        ("Fig5_paired", "Pre-specified descriptive paired differences.", ROOT / "results/stage8/external/tables/paired.tsv"),
        ("Fig6a_cohorts", "Cohort host-attraction rates drawn in Figure 6a.", BMC / "figure_source/Fig6a_cohorts.tsv"),
        ("Fig6b_POG570_strata", "Pre-specified descriptive POG570 site groups.", ROOT / "results/stage4/confirm/secondary_by_site.tsv"),
        ("Fig6b_tumor_content", "Pre-specified descriptive POG570 tumor-content tertiles.", ROOT / "results/stage4/confirm/secondary_tc.tsv"),
        ("Fig6c_aux_top1", "Post hoc auxiliary cohort top-1.", BMC / "figure_source/Fig6c_aux_top1.tsv"),
        ("Fig6d_removed_errors", "Post hoc fate of removed host-attraction errors. No sample identifiers.", BMC / "figure_source/Fig6d_removed.tsv"),
        ("S1_heatmap", "Locked-classifier heatmap source. Two-decimal labels in the figure.", BMC / "figure_source/S1_heatmap.tsv"),
        ("S2_sets", "Post hoc set metrics drawn in Figure S2.", BMC / "figure_source/S2_sets.tsv"),
        ("S3_esophagus_masking", "Post hoc esophagus-call fractions.", BMC / "figure_source/S3_esophagus_fraction.tsv"),
        ("S4_contributions", "Post hoc rank-normal contributions. Not the standardized input.", BMC / "figure_source/S4_contributions.tsv"),
        ("S5_histology", "Post hoc esophagus histology assignment.", BMC / "figure_source/S5_histology.tsv"),
        ("S5_reconciliation", "Figure S5 bar values checked against the stored counts.", BMC / "figure_source/S5_reconciliation.tsv"),
        ("S6_meta_analysis", "Post hoc cohort differences. Wald intervals from the recorded variance.", BMC / "figure_source/S6_meta.tsv"),
        ("Microarray_layer", "Pre-specified descriptive microarray effects. p values are stored and are not presented as tests.", ROOT / "results/stage7/confirm/tables/microarray_effects.tsv"),
        ("Sensitivity", "Auxiliary sensitivity rows.", ROOT / "results/stage7/confirm/tables/sensitivity.tsv"),
        ("Leave_one_cohort", "Pre-specified leave-one-cohort rows for the auxiliary family.", ROOT / "results/stage8/diagnostics/leave_one_cohort.tsv"),
        ("Conformal", "POG570 conformal summary. Aggregate only.", ROOT / "results/stage4/confirm/conformal.tsv"),
        ("Power", "Design power computed before auxiliary unlock.", ROOT / "results/stage6/power.tsv"),
    ]
    readme.append(["Sheet", "Status"])
    for name, description, path in sheets:
        readme.append([name, description])
        add_sheet(book, name, description, read_tsv(path))
    dest = BMC / "Additional_file_2.xlsx"
    book.save(dest)
    print(dest, dest.stat().st_size)


def write_additional_3() -> None:
    source = BMC / "Additional_file_3.xlsx"
    book = load_workbook(source)
    if "Deviations_en" in book.sheetnames:
        del book["Deviations_en"]
    sheet = book.create_sheet("Deviations_en")
    sheet.append(["English list of deviations from the committed plans. The history sheet is unchanged."])
    sheet.append(["item", "what the plan said", "what was done", "effect on the reported result"])
    rows = [
        (
            "Empty biopsy-site fields",
            "The first unlocked auxiliary confirmation was expected to write the hypothesis tables.",
            "That run stopped before writing a hypothesis table because empty site fields were not handled. A later commit changed only that handling. The second run is the reported result.",
            "The reported auxiliary hypotheses are the second run.",
        ),
        (
            "Eighteen GSE50760 liver metastases",
            "The input matrices hashed before the published classifiers were run were expected to contain every auxiliary at-risk biopsy.",
            "Those 18 liver-metastasis samples were scored afterwards, with the same input rule, so the at-risk denominator stayed 427.",
            "Top-1 accuracy uses 729 evaluation samples. The host-attraction rate uses 427 at-risk samples.",
        ),
        (
            "Classifier named in the comparison plan",
            "The comparison plan named a third published classifier.",
            "It was not run. No trained model or gene list was publicly available.",
            "Table 3 contains the baseline, HostMix-TOO, SCOPE and CUP-AI-Dx only.",
        ),
        (
            "MET500 prediction path in the stage-8 plan",
            "The plan names results/stage5/posthoc/met500_predictions.parquet.",
            "That file does not exist. Scoring used results/stage5/posthoc/met500_samples.tsv with results/stage5/posthoc/met500_BASE-Z_proba.npy and results/stage5/posthoc/met500_SA-Z_proba.npy. The plan file was not edited after the result.",
            "The MET500 rows in the published-classifier table use those files.",
        ),
    ]
    for row in rows:
        sheet.append(list(row))
    book.save(source)
    print(source, source.stat().st_size, "sheets", book.sheetnames)


def main() -> None:
    write_additional_1()
    write_additional_2()
    write_additional_3()


if __name__ == "__main__":
    main()
