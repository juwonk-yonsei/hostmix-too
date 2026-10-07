"""Single fact definitions for the A17 check, grouped by the place in the manuscript that needs them.

Each definition names a source file, a selector, a column and a transform. Derived facts name
other facts only. Called from 68_facts_A17.build().
"""
from __future__ import annotations

TCGA = "config/split_tcga.tsv"
GTEX = "config/split_gtex.tsv"
SET9 = "results/stage9/set_metrics.tsv"
EXT = "results/stage8/external/tables/metrics.tsv"


def define(add, g) -> None:
    P = {k: g[k] for k in g if k.startswith("P_")}
    SIM = g["P_SIM"]

    # ------------------------------------------------------------ data sets (Table 1, Methods: Data)
    add("tcga_train_n", TCGA, "split == 'TCGA-train'", "", "count", 0, cohort="TCGA-train", metric="n",
        set_="all", meaning="TCGA-train, number of primary tumors")
    add("tcga_test_n", TCGA, "split == 'TCGA-test'", "", "count", 0, cohort="TCGA-test", metric="n",
        set_="all", meaning="TCGA-test, number of primary tumors")
    add("tcga_learning_n", "results/audit/tcga_summary.json", "all", "n_learning", decimals=0, cohort="TCGA",
        metric="n", set_="all", meaning="TCGA primary tumors kept (one per patient, sample type 01 or 03 for LAML)")
    add("tcga_projects_n", TCGA, "split in ['TCGA-train', 'TCGA-test']", "", "nunique:project", 0, cohort="TCGA",
        metric="count", set_="all", meaning="number of TCGA projects among the kept primary tumors")
    add("n_labels", TCGA, "split == 'TCGA-train'", "", "nunique:learning_label", 0, cohort="TCGA-train",
        metric="count", set_="all", meaning="number of tumor labels (learning labels) in TCGA-train")
    add("n_organs", "config/mappings/project_to_organ.tsv", "all", "", "nunique:organ", 0, metric="count",
        meaning="number of organ groups the tumor labels map to")
    add("gtex_n", GTEX, "all", "", "count", 0, cohort="GTEx", metric="n", set_="all",
        meaning="GTEx samples from tissue donors kept for the study")
    add("gtex_nondonor_n", "results/audit/gtex_nondonor_ids.tsv", "all", "", "count", 0, cohort="GTEx",
        metric="n", meaning="GTEx cell-line samples without a donor, removed")
    add("gtexref_n", GTEX, "split == 'GTEx-ref'", "", "count", 0, cohort="GTEx", metric="n", set_="GTEx-ref",
        meaning="GTEx-ref, number of samples")
    add("gtexsim_n", GTEX, "split == 'GTEx-sim'", "", "count", 0, cohort="GTEx", metric="n", set_="GTEx-sim",
        meaning="GTEx-sim, number of samples")
    add("gtexref_pool_n", GTEX, "split == 'GTEx-ref' and tissue in @POOL", "", "count", 0, cohort="GTEx",
        metric="n", set_="GTEx-ref, host pool", meaning="GTEx-ref samples of the ten host-pool tissues")
    add("gtexsim_22_n", GTEX, "split == 'GTEx-sim' and (tissue in @POOL or tissue in @OUT)", "", "count", 0,
        cohort="GTEx", metric="n", set_="GTEx-sim, 22 tissues",
        meaning="GTEx-sim samples of the 22 simulation tissues")
    add("sim_tissues_n", SIM, "all", "", "nunique:tissue", 0, cohort="sim", metric="count",
        meaning="number of GTEx tissues in the simulation benchmark")
    add("sim_tumors_n", "results/stage5/sim_ext/tumors.tsv", "all", "", "count", 0, cohort="TCGA-test",
        metric="n", meaning="TCGA-test tumors drawn for the simulation benchmark")
    add("imv_kidney_n", "results/stage7/posthoc/imvigor_kidney_meta.json", "all", "n_kidney_in_predictions",
        decimals=0, cohort="IMvigor210", metric="n", set_="pool-out at risk",
        meaning="IMvigor210 metastases biopsied from the kidney (pool-out at-risk biopsies)")

    # ------------------------------------------------------------ Table 2 arms that are not stored in one row
    add("aux_sa_pooloutrisk_host", SET9, "analysis_cohort == 'aux_rnaseq' and set == 'pool_out_risk' and method == 'SA-Z'",
        "host_rate", cohort="aux", method="HostMix-TOO", metric="host-attraction rate", set_="pool-out at risk",
        meaning="auxiliary RNA-seq, HostMix-TOO, host-attraction rate, pool-out at-risk set (post hoc set table)")
    add("aux_gated_native", "", "", "", "derived: ext_aux_rnaseq_sa_native + AS1_diff", cohort="aux",
        method="gated", metric="top-1 accuracy", set_="native truth",
        meaning="auxiliary RNA-seq, gated model, top-1 accuracy, native-truth set (HostMix-TOO arm plus AS1 difference)")
    add("aux_gated_top1", "", "", "", "derived: ext_aux_rnaseq_sa_top1 + AS3_diff", cohort="aux",
        method="gated", metric="top-1 accuracy", set_="evaluation",
        meaning="auxiliary RNA-seq, gated model, top-1 accuracy, evaluation set (HostMix-TOO arm plus AS3 difference)")
    add("aux_pool22_pooloutrisk_host", "", "", "", "derived: aux_sa_pooloutrisk_host + AS2_diff", cohort="aux",
        method="22-tissue", metric="host-attraction rate", set_="pool-out at risk",
        meaning="auxiliary RNA-seq, 22-tissue model, host-attraction rate, pool-out at-risk set (HostMix-TOO arm plus AS2 difference)")
    for h, col in (("AH1", "n_favor"), ("AH1", "n_against"), ("AH2", "n_favor"), ("AH2", "n_against"),
                   ("AS2", "n_favor"), ("AS2", "n_against")):
        src = P["P_AHP"] if h.startswith("AH") else P["P_AHS"]
        add(f"{h}_{col}", src, f"hypothesis == '{h}'", col, decimals=0, cohort="aux",
            method="HostMix-TOO" if h.startswith("AH") else "22-tissue",
            metric="discordant pairs", set_={"AH1": "at risk", "AH2": "evaluation", "AS2": "pool-out at risk"}[h],
            meaning=f"{h}, auxiliary RNA-seq, discordant pairs {'in favor' if col == 'n_favor' else 'against'}")

    # ------------------------------------------------------------ Abstract
    add("pub_host_min_pct", "", "", "", "derived: min(ext_MET500_scope_host_pct, ext_POG570_scope_host_pct, "
        "ext_aux_rnaseq_scope_host_pct, ext_MET500_cup_host_pct, ext_POG570_cup_host_pct, ext_aux_rnaseq_cup_host_pct)",
        decimals=1, cohort="MET500 and POG570 and aux", method="SCOPE and CUP-AI-Dx", metric="host-attraction rate",
        set_="at risk", meaning="lowest host-attraction rate of SCOPE and CUP-AI-Dx over the three cohorts (percent)")
    add("pub_host_max_pct", "", "", "", "derived: max(ext_MET500_scope_host_pct, ext_POG570_scope_host_pct, "
        "ext_aux_rnaseq_scope_host_pct, ext_MET500_cup_host_pct, ext_POG570_cup_host_pct, ext_aux_rnaseq_cup_host_pct)",
        decimals=1, cohort="MET500 and POG570 and aux", method="SCOPE and CUP-AI-Dx", metric="host-attraction rate",
        set_="at risk", meaning="highest host-attraction rate of SCOPE and CUP-AI-Dx over the three cohorts (percent)")
    add("sa_host_min_pct", "", "", "", "derived: min(ext_MET500_sa_host_pct, ext_POG570_sa_host_pct, ext_aux_rnaseq_sa_host_pct)",
        decimals=1, cohort="MET500 and POG570 and aux", method="HostMix-TOO", metric="host-attraction rate",
        set_="at risk", meaning="lowest HostMix-TOO host-attraction rate over the three cohorts (percent)")
    add("sa_host_max_pct", "", "", "", "derived: max(ext_MET500_sa_host_pct, ext_POG570_sa_host_pct, ext_aux_rnaseq_sa_host_pct)",
        decimals=1, cohort="MET500 and POG570 and aux", method="HostMix-TOO", metric="host-attraction rate",
        set_="at risk", meaning="highest HostMix-TOO host-attraction rate over the three cohorts (percent)")
    for lo_hi, fn in (("min", "min"), ("max", "max")):
        add(f"pub_host_{lo_hi}", "", "", "", f"derived: {fn}(ext_MET500_scope_host, ext_POG570_scope_host, "
            "ext_aux_rnaseq_scope_host, ext_MET500_cup_host, ext_POG570_cup_host, ext_aux_rnaseq_cup_host)",
            cohort="MET500 and POG570 and aux", method="SCOPE and CUP-AI-Dx", metric="host-attraction rate",
            set_="at risk", meaning=f"{lo_hi} host-attraction rate of SCOPE and CUP-AI-Dx over the three cohorts")
        add(f"sa_host_{lo_hi}", "", "", "", f"derived: {fn}(ext_MET500_sa_host, ext_POG570_sa_host, ext_aux_rnaseq_sa_host)",
            cohort="MET500 and POG570 and aux", method="HostMix-TOO", metric="host-attraction rate", set_="at risk",
            meaning=f"{lo_hi} HostMix-TOO host-attraction rate over the three cohorts")

    # ------------------------------------------------------------ Additional file 1 tables
    pair = "results/stage2/met500_paired.tsv"
    psel = "method == 'SA' and representation == 'Z'"
    for col, what, set_ in (("n_host_fixed", "host-attraction errors under the baseline only", "at risk"),
                            ("n_host_new", "host-attraction errors under HostMix-TOO only", "at risk"),
                            ("n_method_only", "biopsies correct under HostMix-TOO only", "evaluation"),
                            ("n_base_only", "biopsies correct under the baseline only", "evaluation")):
        add(f"dev_pair_sa_{col}", pair, psel, col, decimals=0, cohort="MET500", method="HostMix-TOO",
            metric="discordant pairs", set_=set_, meaning=f"MET500 development, {what}")
    add("pw_h1_half_favor", "", "", "", "derived: dev_pair_sa_n_host_fixed / 2", decimals=0, cohort="MET500",
        method="HostMix-TOO", metric="discordant pairs", set_="at risk",
        meaning="half-effect H1 power scenario, numerator of p_favor (half of 34)")
    add("pw_h2_half_den", "", "", "", "derived: dev_base_n * 2", decimals=0, cohort="MET500", metric="n",
        set_="evaluation", meaning="half-effect H2 power scenario, denominator of p_favor (twice 437)")
    log = "results/stage5/variants/train_log.tsv"
    add("var_m8_rows_per_tumor", log, "model == 'SA-m8'", "n_per_tumor", decimals=0, method="SA-m8",
        metric="training rows per tumor", meaning="SA-m8 training rows per tumor (one pure profile plus mixtures)")
    add("var_m8_nmix", "", "", "", "derived: var_m8_rows_per_tumor - 1", decimals=0, method="SA-m8",
        metric="mixtures per tumor", meaning="SA-m8 mixtures per tumor")

    # ------------------------------------------------------------ figure texts (A17 section 4.5)
    pipe = "results/stage9/ablation/pipeline_table.tsv"
    add("sa_train_rows", pipe, "model == 'SA-Z'", "n_train_rows", decimals=0, method="HostMix-TOO", metric="n",
        meaning="HostMix-TOO training rows (pure profiles plus mixtures)")
    add("sa_c", pipe, "model == 'SA-Z'", "C", form="raw", method="HostMix-TOO", metric="C",
        meaning="regularization constant C of HostMix-TOO")
    add("base_c", pipe, "model == 'BASE-Z'", "C", form="raw", method="baseline", metric="C",
        meaning="regularization constant C of the baseline")
    add("sa_n_features", pipe, "model == 'SA-Z'", "n_features", decimals=0, method="HostMix-TOO", metric="genes",
        meaning="classifier genes of HostMix-TOO")
    add("aux_cohorts_n", "results/stage9/subcohort_metrics.tsv", "all", "", "nunique:cohort", 0, cohort="aux",
        metric="count", meaning="auxiliary RNA-seq cohorts of confirmation 2")
    add("sim_outside_n", SIM, "tissue in @OUT", "", "nunique:tissue", 0, cohort="sim", metric="count",
        meaning="simulated host tissues outside the training host pool")
    add("ma_cohorts_n", "results/stage7/confirm/tables/metrics.tsv", "layer == '\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774' and cohort != 'all'",
        "", "nunique:cohort", 0, cohort="microarray", metric="count", meaning="microarray cohorts of the descriptive layer")

    # ------------------------------------------------------------ Methods
    add("planted_compare_n", "results/stage7/planted/compare.tsv", "match == True", "", "count", 0, cohort="aux",
        metric="count", meaning="auxiliary-family comparisons of the independent checker on planted predictions that matched")
    add("split_test_pct", "", "", "", "derived: split_test_frac * 100", decimals=0, cohort="TCGA", metric="percent",
        meaning="TCGA-test share of the split (percent)")
    add("split_train_pct", "", "", "", "derived: 100 - split_test_frac * 100", decimals=0, cohort="TCGA",
        metric="percent", meaning="TCGA-train share of the split (percent)")
    add("met500_unmapped_n", "config/met500_eval_samples.tsv", "cohort in ['MISC', 'SECR']", "", "count", 0,
        cohort="MET500", metric="n", meaning="MET500 specimens of the unmapped cohorts MISC and SECR, excluded")
    add("pog570_profiles_n", "config/pog570_eval_labels.tsv", "all", "", "count", 0, cohort="POG570", metric="n",
        meaning="POG570 expression profiles (one per patient)")
    add("pog570_unmapped_n", "config/pog570_eval_labels.tsv", "organ == 'exclude'", "", "count", 0, cohort="POG570",
        metric="n", meaning="POG570 samples whose diagnosis could not be mapped to an organ group, excluded")
    add("su2c_dup_removed_n", "results/stage10/derived/duplicate_removed.tsv", "cohort == 'prad_su2c_2019'",
        "n_removed", decimals=0, cohort="SU2C", metric="n", note="post hoc recount, A17 section 4.7",
        meaning="SU2C samples removed by the duplicate screen (post hoc recount, both libraries)")
    add("gse50760_late_n", "results/stage8/external/scope/GSE50760_risk/SCOPE_topPredictions.txt", "all", "",
        "nunique:sample_name", 0, cohort="GSE50760", metric="n",
        meaning="GSE50760 at-risk liver-metastasis samples scored by the published classifiers afterwards")
    add("kidney_projects_n", "config/mappings/project_to_organ.tsv", "organ == 'Kidney'", "", "count", 0,
        metric="count", meaning="TCGA kidney cancer projects mapped to the kidney organ group")
    proxy = "results/stage2/proxy_correlation.tsv"
    add("proxy_ln_spleen", proxy, "hpa_query == 'lymph node' and gtex_tissue == 'Spleen'", "spearman",
        metric="Spearman correlation", meaning="Spearman correlation of spleen with the lymph node consensus profile")
    add("proxy_bm_blood", proxy, "hpa_query == 'bone marrow' and gtex_tissue == 'Whole Blood'", "spearman",
        metric="Spearman correlation", meaning="Spearman correlation of whole blood with the bone marrow consensus profile")
    add("sa_mix_rows", "", "", "", "derived: sa_train_rows - tcga_train_n", decimals=0, method="HostMix-TOO",
        metric="n", meaning="HostMix-TOO training mixtures")
    add("sa_rows_factor", "", "", "", "derived: sa_train_rows / tcga_train_n", decimals=0, method="HostMix-TOO",
        metric="ratio", meaning="HostMix-TOO training rows per pure tumor")
    add("var_c001", log, "model == 'SA-C0.01'", "C", form="raw", method="SA-C0.01", metric="C",
        meaning="C of the exploratory variant SA-C0.01")
    add("var_c01", log, "model == 'SA-C0.1'", "C", form="raw", method="SA-C0.1", metric="C",
        meaning="C of the exploratory variant SA-C0.1")
    add("var_loho_n", log, "model.str.startswith('SA-LOHO-')", "", "count", 0, method="LOHO", metric="count",
        meaning="leave-one-host-out variant models")
    add("pool22_outside_n", "", "", "", "derived: pool22_n - host_pool_n", decimals=0, method="22-tissue",
        metric="count", meaning="tissues of the 22-tissue host pool outside the ten-tissue host pool")
    add("cup_example_correct", "", "", "", "derived: round(cup_example_acc * cup_example_n / 100)", decimals=0,
        method="CUP-AI-Dx", metric="n", meaning="CUP-AI-Dx example, correctly classified metastatic samples")
    add("cup_example_acc_pct", "", "", "", "derived: cup_example_acc", decimals=1, method="CUP-AI-Dx",
        metric="percent", meaning="CUP-AI-Dx example, overall metastatic accuracy (percent)")
    add("abl_shares_total", "", "", "", "derived: abl_shares_n * abl_cohorts_n", decimals=0, metric="count",
        meaning="shares of the ablation reading rule (models times cohorts)")
    add("regen_coef_max", "results/stage9/ablation/regeneration_check.json", "all", "coef_max_abs", decimals=2,
        form="sci", method="HostMix-TOO", metric="max abs coefficient difference",
        meaning="largest absolute coefficient difference between the regenerated and locked HostMix-TOO")
    esca = "results/stage9/esca_histology_counts.tsv"
    add("esca_train_n", esca, "histology == 'adenocarcinoma'", "n_train", decimals=0, cohort="TCGA-train",
        metric="n", meaning="TCGA-train esophageal carcinomas")
    for h in ("adenocarcinoma", "squamous", "other"):
        add(f"esca_{h}_n", esca, f"histology == '{h}'", "n", decimals=0, cohort="TCGA-train", metric="n",
            meaning=f"TCGA-train esophageal carcinomas, {h}")


    # ------------------------------------------------------------ Results
    tissues = ["adiposesc", "adiposevis", "adrenal", "blood", "brain", "liver", "lung", "muscle", "skin", "spleen"]
    for arm, what in (("loho", "leave-one-host-out model on its left-out tissue"), ("base", "baseline"),
                      ("sa", "HostMix-TOO")):
        add(f"loho_{arm}_median_06", "", "", "",
            "derived: median(" + ", ".join(f"simt_{t}_06_{arm}_host" for t in tissues) + ")",
            method={"loho": "leave-one-host-out", "base": "baseline", "sa": "HostMix-TOO"}[arm],
            metric="host-attraction rate", meaning=f"simulation, rho 0.6, median host-attraction rate over the ten "
            f"host-pool tissues, {what}")
    for site, key in (("soft_tissue", "st"), ("remainder", "rem")):
        for arm in ("base", "sa"):
            add(f"pog_{key}_{arm}_ncorrect", "", "", "",
                f"derived: round(pogsite_{site}_{arm}_native * pogsite_{site}_{arm}_nnative)", decimals=0,
                cohort="POG570", method={"base": "baseline", "sa": "HostMix-TOO"}[arm], metric="n correct",
                set_="native truth", meaning=f"POG570 {site} biopsies, correct calls in the native-truth set")
    add("pog_native_net_loss", "", "", "", "derived: round((pog_base_native - pog_sa_native) * pog_base_nnative)",
        decimals=0, cohort="POG570", method="HostMix-TOO", metric="n correct", set_="native truth",
        meaning="POG570 native-truth set, correct calls lost by HostMix-TOO relative to the baseline (net)")
    shift = "results/stage8/diagnostics/error_shift_summary.tsv"
    add("prince_ah1_favor", shift, "cohort == 'paad_iatlas_prince_2022'", "n_favor", decimals=0, cohort="PRINCE",
        method="HostMix-TOO", metric="discordant pairs", set_="at risk",
        meaning="PRINCE, AH1 discordant pairs in favor of HostMix-TOO")
    add("su2c_ah1_favor", shift, "cohort == 'prad_su2c_2019'", "", "sum:n_favor", 0, cohort="SU2C",
        method="HostMix-TOO", metric="discordant pairs", set_="at risk",
        meaning="SU2C, AH1 discordant pairs in favor of HostMix-TOO (real and proxy host sites)")
    add("prince_risk_liver_n", "results/stage5/aux/cohorts/paad_iatlas_prince_2022/set_counts_after_flag.tsv",
        "standard_site == 'liver' and set == 'risk'", "n_patients", decimals=0, cohort="PRINCE", metric="n",
        set_="at risk", meaning="PRINCE at-risk biopsies taken from the liver (patients)")
    sens = "results/stage7/confirm/tables/sensitivity.tsv"
    aux6 = ["GSE50760", "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
            "paad_iatlas_prince_2022", "prad_su2c_2019"]
    for i, c in enumerate(aux6):
        add(f"sens_lo{i}_AH1_diff", sens, f"analysis == 'leave_out_{c}_AH1'", "diff", cohort="aux",
            method="HostMix-TOO", metric="host-attraction rate", set_="at risk",
            meaning=f"auxiliary RNA-seq without {c}, AH1 difference (pre-specified leave-one-cohort)")
    add("sens_lo_AH1_diff_max", "", "", "", "derived: max(" + ", ".join(f"sens_lo{i}_AH1_diff" for i in range(6)) + ")",
        cohort="aux", method="HostMix-TOO", metric="host-attraction rate", set_="at risk",
        meaning="auxiliary RNA-seq, leave-one-cohort AH1 differences, largest (closest to zero)")
    add("sens_lo_AH1_diff_min", "", "", "", "derived: min(" + ", ".join(f"sens_lo{i}_AH1_diff" for i in range(6)) + ")",
        cohort="aux", method="HostMix-TOO", metric="host-attraction rate", set_="at risk",
        meaning="auxiliary RNA-seq, leave-one-cohort AH1 differences, smallest")
    add("sens_lo_imv_AH2_diff", sens, "analysis == 'leave_out_blca_iatlas_imvigor210_2017_AH2'", "diff", cohort="aux",
        method="HostMix-TOO", metric="top-1 accuracy", set_="evaluation",
        meaning="auxiliary RNA-seq without IMvigor210, AH2 difference (pre-specified leave-one-cohort)")
    add("sens_lo_imv_AH2_n", sens, "analysis == 'leave_out_blca_iatlas_imvigor210_2017_AH2'", "n", decimals=0,
        cohort="aux", method="HostMix-TOO", metric="n", set_="evaluation",
        meaning="auxiliary RNA-seq without IMvigor210, AH2 n (pre-specified leave-one-cohort)")
    top = "results/stage8/diagnostics/imvigor_top_classes.tsv"
    for arm, m in (("base", "BASE-Z"), ("sa", "SA-Z")):
        add(f"imv_native_{arm}_esophagus", top, f"slice == 'native_truth' and method == '{m}' and predicted == 'Esophagus'",
            "n", decimals=0, cohort="IMvigor210", method={"base": "baseline", "sa": "HostMix-TOO"}[arm],
            metric="esophagus calls", set_="native truth",
            meaning="IMvigor210 bladder tumors biopsied from the bladder, esophagus calls")
    add("imv_native_n", top, "slice == 'native_truth' and method == 'BASE-Z' and predicted == 'Esophagus'", "n_slice",
        decimals=0, cohort="IMvigor210", metric="n", set_="native truth",
        meaning="IMvigor210 bladder tumors biopsied from the bladder")
    gate = "results/stage10/derived/gate_decisions.tsv"
    for set_name, key in (("evaluation", "eval"), ("native truth", "native")):
        add(f"gate_base_{key}_n", gate, f"set == '{set_name}'", "n_baseline", decimals=0, cohort="aux", method="gated",
            metric="n", set_=set_name, note="post hoc recount, A17 section 4.7",
            meaning=f"auxiliary RNA-seq, {set_name} samples sent to the baseline by the gate (post hoc count)")
    rates = "results/stage10/derived/microarray_rates.tsv"
    add("ma_rates_ld_top1", rates, "cohort == 'all' and method == 'LD-Z'", "top1", cohort="microarray",
        method="deconvolution", metric="top-1 accuracy", set_="evaluation",
        meaning="microarray layer, linear deconvolution, top-1 accuracy, evaluation set (post hoc)")
    add("ma_rates_ld_host", rates, "cohort == 'all' and method == 'LD-Z'", "host_rate", cohort="microarray",
        method="deconvolution", metric="host-attraction rate", set_="at risk",
        meaning="microarray layer, linear deconvolution, host-attraction rate, at-risk set (post hoc)")
    shares = [f"share_{c}_{m}" for c in ("MET500", "POG570", "aux_rnaseq") for m in ("mix", "pure03", "pure15")]
    add("share_min", "", "", "", "derived: min(" + ", ".join(shares) + ")", cohort="MET500 and POG570 and aux",
        metric="share", meaning="ablation, smallest share of the host-attraction gap (nine shares)")
    add("share_max", "", "", "", "derived: max(" + ", ".join(shares) + ")", cohort="MET500 and POG570 and aux",
        metric="share", meaning="ablation, largest share of the host-attraction gap (nine shares)")
    add("abl_ma_models_n", "results/stage10/ablation_microarray/metrics.tsv", "all", "", "nunique:model", 0,
        metric="count", meaning="models of the microarray ablation")
    var_hosts = ["php_mx1_host", "php_m8_host", "php_r30_host", "php_r50_host",
                 "php_c001_host", "php_c01_host", "php_samlp_host", "php_sa_host"]
    add("php_var_host_min", "", "", "", "derived: min(" + ", ".join(var_hosts) + ")", cohort="POG570",
        method="variants", metric="host-attraction rate", set_="at risk",
        meaning="POG570, lowest host-attraction rate over HostMix-TOO and the setting variants")
    add("php_var_host_max", "", "", "", "derived: max(" + ", ".join(var_hosts) + ")", cohort="POG570",
        method="variants", metric="host-attraction rate", set_="at risk",
        meaning="POG570, highest host-attraction rate over HostMix-TOO and the setting variants")
    ext = "results/stage8/external/tables/metrics.tsv"
    sel = "analysis_cohort == 'aux_rnaseq' and subset == 'full' and method == 'SCOPE'"
    add("scope_aux_ns_n", ext, sel, "n_scope_ns_top1", decimals=0, cohort="aux", method="SCOPE", metric="n",
        set_="evaluation", meaning="auxiliary RNA-seq, SCOPE top-1 calls of a normal-tissue class")
    add("scope_aux_top1_ns_organ", ext, sel, "top1_ns_as_organ", cohort="aux", method="SCOPE", metric="top-1 accuracy",
        set_="evaluation", meaning="auxiliary RNA-seq, SCOPE top-1 accuracy with normal-tissue calls counted as the organ")
    tot = "results/stage8/diagnostics/error_shift_totals.tsv"
    for c, key in (("POG570", "POG570"), ("MET500", "MET500"), ("aux_rnaseq", "aux")):
        for col in ("n_favor", "n_sa_correct"):
            add(f"shift_{c}_{col}", tot, f"analysis == '{c}'", col, decimals=0, cohort=key, method="HostMix-TOO",
                metric="n", set_="at risk", meaning=f"{c}, baseline host-attraction errors removed by HostMix-TOO, {col}")
    cup = "results/stage9/cup_correct_sa_wrong.tsv"
    add("cupsa_total", cup, "field == 'cohort' and level == 'blca_iatlas_imvigor210_2017'", "n_cup_correct_sa_wrong",
        decimals=0, cohort="aux", method="CUP-AI-Dx and HostMix-TOO", metric="n", set_="evaluation",
        meaning="auxiliary RNA-seq, samples correct under CUP-AI-Dx and not under HostMix-TOO")
    for level, key, ck in (("blca_iatlas_imvigor210_2017", "imv", "IMvigor210"), ("brca_iatlas_anders_2022", "anders", "Anders")):
        add(f"cupsa_{key}", cup, f"field == 'cohort' and level == '{level}'", "n", decimals=0, cohort=ck,
            method="CUP-AI-Dx and HostMix-TOO", metric="n", set_="evaluation",
            meaning=f"{ck}, samples correct under CUP-AI-Dx and not under HostMix-TOO")
    add("cupsa_sa_esophagus", cup, "field == 'SA-Z' and level == 'Esophagus'", "n", decimals=0, cohort="aux",
        method="HostMix-TOO", metric="esophagus calls", set_="evaluation",
        meaning="auxiliary RNA-seq, CUP-correct HostMix-TOO-wrong samples that HostMix-TOO called esophagus")
    colon = "results/stage8/diagnostics/gse41258_colon_confusion.tsv"
    add("gse41258_sa_esophagus", colon, "method == 'SA-Z' and predicted == 'Esophagus'", "n", decimals=0,
        cohort="GSE41258", method="HostMix-TOO", metric="esophagus calls",
        meaning="GSE41258 primary colorectal tumors, HostMix-TOO esophagus calls")
    add("gse41258_colon_n", colon, "method == 'SA-Z' and predicted == 'Esophagus'", "n_slice", decimals=0,
        cohort="GSE41258", metric="n", meaning="GSE41258 primary colorectal tumors")
    hist = "results/stage9/esophagus_histology_summary.tsv"
    for c, key, ck in (("GSE41258", "gse41258", "GSE41258"), ("paad_iatlas_prince_2022", "prince", "PRINCE"),
                       ("blca_iatlas_imvigor210_2017", "imv", "IMvigor210")):
        for col in ("n", "n_closer_adenocarcinoma", "n_closer_squamous"):
            add(f"hist_{key}_sa_{col}", hist, f"method == 'SA-Z' and cohort == '{c}'", col, decimals=0, cohort=ck,
                method="HostMix-TOO", metric="esophagus calls",
                meaning=f"{ck}, HostMix-TOO esophagus calls, {col}")
    add("gpl96_pog570_missing_n", "results/stage9/mask_feature_info.tsv",
        "mask == 'GPL96' and cohort == 'POG570' and rep == '0'", "n_B0_dropped", decimals=0, cohort="POG570",
        metric="n", meaning="classifier genes not measured on GPL96 (POG570 restriction)")
    for arm in ("base", "sa"):
        add(f"mg_{arm}_gpl96_pog570_extra", "", "", "",
            f"derived: mg_{arm}_gpl96_pog570_restriction_drop - mg_{arm}_gpl96_pog570_random_mean", cohort="POG570",
            method={"base": "baseline", "sa": "HostMix-TOO"}[arm], metric="top-1 loss",
            meaning="POG570, GPL96 restriction loss minus the mean random-removal loss")
    prince = [f"sub_prince_{m}_top1" for m in ("base", "sa", "scope", "cup")]
    add("prince_top1_min", "", "", "", "derived: min(" + ", ".join(prince) + ")", cohort="PRINCE",
        method="any", metric="top-1 accuracy", set_="evaluation", meaning="PRINCE, lowest top-1 accuracy of the four methods")
    add("prince_top1_max", "", "", "", "derived: max(" + ", ".join(prince) + ")", cohort="PRINCE",
        method="any", metric="top-1 accuracy", set_="evaluation", meaning="PRINCE, highest top-1 accuracy of the four methods")

    # ------------------------------------------------------------ Discussion, legends, table notes
    for lo_hi, fn in (("min", "min"), ("max", "max")):
        add(f"sa_base_host_ratio_{lo_hi}", "", "", "",
            f"derived: {fn}(ext_MET500_sa_host / ext_MET500_base_host, ext_POG570_sa_host / ext_POG570_base_host, "
            "ext_aux_rnaseq_sa_host / ext_aux_rnaseq_base_host)", decimals=3, cohort="MET500 and POG570 and aux",
            method="HostMix-TOO", metric="host-attraction rate ratio", set_="at risk", approx="0.01",
            meaning=f"{lo_hi} over the three RNA-seq cohorts of the HostMix-TOO to baseline host-attraction ratio")
    add("aux6_nc_host", "results/stage10/derived/aux_six_rates.tsv", "method == 'NC-Z'", "host_rate", cohort="aux",
        method="normal classes", metric="host-attraction rate", set_="at risk",
        meaning="auxiliary RNA-seq (six cohorts pooled), normal classes, host-attraction rate (post hoc)")
    nc_scope = "dev_nc_host, pog_nc_host, aux6_nc_host, ext_MET500_scope_host, ext_POG570_scope_host, ext_aux_rnaseq_scope_host"
    for lo_hi, fn in (("min", "min"), ("max", "max")):
        add(f"nc_scope_host_{lo_hi}_pct", "", "", "", f"derived: {fn}({nc_scope}) * 100", decimals=0,
            cohort="MET500 and POG570 and aux", method="normal classes and SCOPE", metric="host-attraction rate",
            set_="at risk", meaning=f"{lo_hi} host-attraction rate of the normal-class model and SCOPE over the "
            "three RNA-seq cohorts (percent)")
    add("pog_liver_share", "", "", "", "derived: pogsite_liver_base_n / pog_base_n", decimals=3, cohort="POG570",
        metric="share", set_="evaluation", approx="0.05",
        meaning="POG570, share of evaluation biopsies taken from the liver ('about a third')")
    add("geo_held_n", "results/stage10/derived/screening_decisions.tsv", "path == 'GEO' and decision == '\ubcf4\ub958'",
        "n_studies", decimals=0, metric="n", note="post hoc recount, A17 section 4.7",
        meaning="GEO series held after a title-only screen (last decision per study)")

    # ------------------------------------------------------------ Additional file 1
    slog, sdec = "results/stage5/aux/screening_log.tsv", "results/stage10/derived/screening_decisions.tsv"
    add("screen_rows_n", slog, "all", "", "count", 0, metric="n", meaning="rows of the cohort search log")
    add("screen_ids_n", slog, "all", "", "nunique:id", 0, metric="n", meaning="study identifiers of the cohort search log")
    add("screen_extra_rows_n", "", "", "", "derived: screen_rows_n - screen_ids_n", decimals=0, metric="n",
        meaning="search-log rows that update an identifier already in the log")
    for fid, sel, what in (("screen_excluded_n", "decision == '\uc81c\uc678'", "excluded"),
                           ("screen_selected_n", "decision == '\uc120\ubcc4'", "passed the biopsy-site screen"),
                           ("screen_included_n", "decision == '\ud3ec\ud568'", "marked included"),
                           ("screen_forward_n", "decision in ['\ud3ec\ud568', '\uc120\ubcc4']", "taken forward")):
        add(fid, sdec, sel, "", "sum:n_studies", 0, metric="n", note="post hoc recount, A17 section 4.7",
            meaning=f"studies {what} (last decision per study, GEO and cBioPortal)")
    for crit in (2, 3, 4, 6, 7):
        add(f"screen_crit{crit}_n", "results/stage10/derived/screening_criteria.tsv", f"failed_criterion == {crit}",
            "n_studies", decimals=0, metric="n", note="post hoc recount, A17 section 4.7",
            meaning=f"excluded studies that failed eligibility criterion {crit} (last decision per study)")
    labels = "config/aux_eval_labels.tsv"
    for fid, cohort, reason, label in (
            ("excl8_aurora_n", "GSE209998", "organ_rule:8", "GSE209998"), ("excl8_gse41258_n", "GSE41258", "organ_rule:8", "GSE41258"),
            ("excl8_gse50760_n", "GSE50760", "organ_rule:8", "GSE50760"), ("excl8_gse71729_n", "GSE71729", "organ_rule:8", "GSE71729"),
            ("excl2_gse41258_n", "GSE41258", "organ_rule:2", "GSE41258"), ("excl2_gse71729_n", "GSE71729", "organ_rule:2", "GSE71729"),
            ("excl1_su2c_n", "prad_su2c_2019", "organ_rule:1", "SU2C"), ("excl1_gse74685_n", "GSE74685", "organ_rule:1", "GSE74685"),
            ("exclc_imv_n", "blca_iatlas_imvigor210_2017", "organ_rule:cohort", "IMvigor210"),
            ("exclc_prince_n", "paad_iatlas_prince_2022", "organ_rule:cohort", "PRINCE"),
            ("exclc_anders_n", "brca_iatlas_anders_2022", "organ_rule:cohort", "Anders")):
        add(fid, labels, f"cohort == '{cohort}' and exclude_reason == '{reason}'", "", "count", 0, cohort=label,
            metric="n", meaning=f"{label}, samples excluded with reason {reason}")
    add("exclc_total_n", labels, "exclude_reason == 'organ_rule:cohort'", "", "count", 0, cohort="aux", metric="n",
        meaning="auxiliary samples that kept the cohort organ label without a standard biopsy site")
    add("proxy_corr_n", proxy, "all", "", "count", 0, metric="count",
        meaning="Spearman correlations computed to choose the lymph node and bone marrow proxies")
    add("sa_mix_seed", "results/stage2/augmentation_seeds.json", "all", "seeds.SA", decimals=0, method="HostMix-TOO",
        metric="seed", meaning="mixture seed recorded for HostMix-TOO")
    add("dev_candidates_n", "results/stage2/candidate_rank.tsv", "all", "", "count", 0, cohort="MET500",
        metric="count", meaning="classifiers in the development candidate ranking")
    loho = ", ".join(f"php_loho{t}_host" for t in ("adiposesc", "adiposevis", "adrenal", "blood", "brain", "liver",
                                                     "lung", "muscle", "skin", "spleen"))
    for lo_hi in ("min", "max"):
        add(f"loho_host_{lo_hi}", "", "", "", f"derived: {lo_hi}({loho})", cohort="POG570", method="LOHO",
            metric="host-attraction rate", set_="at risk",
            meaning=f"{lo_hi} POG570 host-attraction rate over the ten leave-one-host-out models")
    add("power_alpha", "results/stage3/power.json", "all", "alpha", decimals=2, metric="alpha",
        meaning="one-sided alpha of the design power simulation")
    add("power_nsim", "results/stage3/power.json", "all", "n_sim", decimals=0, metric="n",
        meaning="simulations per scenario of the design power simulation")
    add("power_seed", "results/stage3/power.json", "all", "seed", decimals=0, metric="seed",
        meaning="seed of the design power simulation")
    add("planted_rows_n", "results/stage7/planted/compare.tsv", "all", "", "count", 0, cohort="aux", metric="n",
        meaning="rows of the independent-checker comparison on planted predictions")
    add("planted_mismatch_n", "results/stage7/planted/compare.tsv", "match == False", "", "count", 0, cohort="aux",
        metric="n", meaning="mismatches of the independent-checker comparison on planted predictions")
    scale = "results/stage8/external/scope/scale_check.tsv"
    for tag, cohort in (("pog", "POG570"), ("imv", "blca_iatlas_imvigor210_2017"), ("anders", "brca_iatlas_anders_2022"),
                        ("prince", "paad_iatlas_prince_2022"), ("dfci", "mel_dfci_2019")):
        add(f"scale_{tag}_frac", scale, f"cohort == '{cohort}'", "top1_changed_fraction", cohort=cohort, method="SCOPE",
            metric="fraction", meaning=f"{cohort}, fraction of samples whose SCOPE top-1 class changed in the doubled run")
    add("scale_dfci_n", scale, "cohort == 'mel_dfci_2019'", "n", decimals=0, cohort="DFCI melanoma", method="SCOPE",
        metric="n", meaning="DFCI melanoma, samples of the SCOPE scale check")
    add("scale_dfci_changed", "", "", "", "derived: round(scale_dfci_frac * scale_dfci_n)", decimals=0,
        cohort="DFCI melanoma", method="SCOPE", metric="n",
        meaning="DFCI melanoma, samples whose SCOPE top-1 class changed in the doubled run")
    add("scale_four_max", "", "", "", "derived: max(scale_pog_frac, scale_imv_frac, scale_anders_frac, scale_prince_frac)",
        decimals=0, cohort="POG570, IMvigor210, Anders, PRINCE", method="SCOPE", metric="fraction",
        meaning="largest changed fraction of the four cohorts other than DFCI melanoma")
    add("m1_host_max", "", "", "", "derived: max(dev_m1_host, pog_m1_host)", decimals=0, method="masking",
        metric="host-attraction rate", set_="at risk", meaning="host-attraction rate of native-organ masking")
    add("set9_sets_n", SET9, "all", "", "nunique:set", 0, metric="count",
        meaning="analysis sets of the post hoc set metrics (Figure S2)")
    add("ma_eval_samples_n", "results/stage7/confirm/tables/metrics.tsv",
        "layer == '\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774' and cohort == 'all' and standard_site == 'all' and method == 'SA-Z'",
        "n", decimals=0, cohort="microarray", metric="n", set_="evaluation",
        meaning="microarray evaluation samples, one per patient within each cohort")
    over = "results/stage10/derived/patient_overlap.tsv"
    add("overlap_eval_n", over, "set == 'evaluation'", "n_both", decimals=0, cohort="GSE74685 and FHCRC", metric="n",
        note="post hoc recount, A17 section 4.7", meaning="patients evaluated in both GSE74685 and FHCRC")
    add("overlap_risk_n", over, "set == 'at risk'", "n_both", decimals=0, cohort="GSE74685 and FHCRC", metric="n",
        note="post hoc recount, A17 section 4.7", meaning="at-risk patients of GSE74685 also at risk in FHCRC")
    meta = "results/stage7/confirm/tables/meta_analysis.tsv"
    add("meta_gse74685_n", meta, "cohort == 'GSE74685'", "n", decimals=0, cohort="GSE74685", metric="n",
        set_="at risk", meaning="GSE74685 at-risk n of the cohort meta-analysis")
    add("meta_rna_k", meta, "status == '\uacc4\uc0b0' and layer == 'RNA-seq'", "k", decimals=0, metric="count",
        meaning="RNA-seq cohorts of the meta-analysis")
    add("meta_ma_k", meta, "status == '\uacc4\uc0b0' and layer == '\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774'", "k",
        decimals=0, metric="count", meaning="microarray cohorts of the meta-analysis")
    contrib = "results/stage9/imvigor_contribution_fix1.tsv"
    for tag, method, label in (("base", "BASE-Z", "baseline"), ("sa", "SA-Z", "HostMix-TOO")):
        add(f"contrib_{tag}_n", contrib, f"method == '{method}' and rank == 1", "n_samples", decimals=0,
            cohort="IMvigor210", method=label, metric="n", set_="native truth",
            meaning=f"IMvigor210 native-truth bladder tumors called esophagus by the {label} (gene contributions)")
    add("abl_rna_models_n", "results/stage9/ablation/metrics.tsv", "all", "", "nunique:method", 0, metric="count",
        meaning="logistic models of the post hoc ablation")
    check = "manuscript/checks/posthoc_model_check.json"
    add("mixmodel_check_n", check, "all", "n_met500", decimals=0, cohort="MET500", method="MIX-Z0", metric="n",
        meaning="MET500 biopsies of the released unstandardized mixture model check")
    add("mixmodel_check_mismatch", check, "all", "joblib_top1_mismatch_vs_stage9", decimals=0, cohort="MET500",
        method="MIX-Z0", metric="n", meaning="top-1 mismatches of the released unstandardized mixture model check")
    for lo_hi in ("min", "max"):
        add(f"pub_host_{lo_hi}_pct0", "", "", "", f"derived: pub_host_{lo_hi} * 100", decimals=0,
            cohort="MET500 and POG570 and aux", method="SCOPE and CUP-AI-Dx", metric="host-attraction rate",
            set_="at risk", meaning=f"{lo_hi} host-attraction rate of the published classifiers (percent, integer)")
