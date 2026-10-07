"""Check the numbers of the A16 replacement text against the stored results.

Every value is read by key. The baseline and HostMix-TOO microarray values are also read from
the stage-7 confirmation predictions, with the one-patient-per-layer sets of the stage-10 code,
and their paired differences and discordant counts are compared with the stage-7 microarray table.
"""
from __future__ import annotations

from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S10 = ROOT / "results" / "stage10" / "ablation_microarray"
S9 = ROOT / "results" / "stage9" / "ablation" / "metrics.tsv"
OUT = ROOT / "manuscript" / "checks" / "numbers_A16.tsv"
MICROARRAY = "\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774"
NAMES = {"BASE-Z": "baseline", "PURE-Zs": "standardized pure C = 0.03", "PURE-Zs-w": "standardized pure C = 0.15",
         "MIX-Z0": "unstandardized mixture", "SA-Z": "HostMix-TOO"}


def show(value: float, digits: int = 3) -> str:
    text = str(Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP))
    return text.replace("-", "\u2212")


def one(frame: pd.DataFrame, **keys) -> pd.Series:
    mask = np.ones(len(frame), dtype=bool)
    for column, value in keys.items():
        mask &= frame[column].astype(str).eq(str(value)).to_numpy()
    if mask.sum() != 1:
        raise SystemExit(f"{keys}: {mask.sum()} rows")
    return frame.loc[mask].iloc[0]


def as_bool(series: pd.Series) -> pd.Series:
    return series.map(lambda value: str(value).strip().lower() in {"true", "1"})


def dedupe(frame: pd.DataFrame) -> pd.DataFrame:
    if frame["patient_id"].is_unique:
        return frame
    part = frame.assign(_rank=np.where(frame["library"].fillna("").eq("polyA"), 0, 1))
    part = part.sort_values(["_rank", "sample_id", "cohort"], kind="mergesort")
    return part.drop_duplicates("patient_id", keep="first").drop(columns="_rank")


def stage7_microarray() -> tuple[dict, list[dict]]:
    columns = ["cohort", "patient_id", "sample_id", "library", "organ", "native_organs", "layer", "layer_test",
               "selected_eval", "selected_risk", "selected_native", "BASE-Z__pred", "SA-Z__pred"]
    pred = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet", columns=columns)
    part = pred.loc[pred["layer"].eq(MICROARRAY)].copy()
    for column in ("sample_id", "patient_id", "library", "organ", "native_organs"):
        part[column] = part[column].fillna("").astype(str)
    tested = part.loc[as_bool(part["layer_test"])]
    sets = {name: dedupe(tested.loc[as_bool(tested[flag])])
            for name, flag in (("evaluation", "selected_eval"), ("at_risk", "selected_risk"), ("native_truth", "selected_native"))}
    values, checks = {}, []
    effects = pd.read_csv(ROOT / "results/stage7/confirm/tables/microarray_effects.tsv", sep="\t").set_index("hypothesis")
    for name, hypothesis in (("evaluation", "AH2"), ("at_risk", "AH1"), ("native_truth", "AH3")):
        sub = sets[name]
        truth = sub["organ"].to_numpy(dtype=object)
        if name == "at_risk":
            natives = [set(filter(None, text.split("|"))) for text in sub["native_organs"]]
            hit = {m: np.array([p in natives[i] for i, p in enumerate(sub[f"{m}__pred"])]) for m in ("BASE-Z", "SA-Z")}
            favor, against = int((hit["BASE-Z"] & ~hit["SA-Z"]).sum()), int((~hit["BASE-Z"] & hit["SA-Z"]).sum())
        else:
            hit = {m: sub[f"{m}__pred"].to_numpy(dtype=object) == truth for m in ("BASE-Z", "SA-Z")}
            favor, against = int((~hit["BASE-Z"] & hit["SA-Z"]).sum()), int((hit["BASE-Z"] & ~hit["SA-Z"]).sum())
        for method in ("BASE-Z", "SA-Z"):
            values[(name, method)] = (float(hit[method].mean()), len(sub))
        stored = effects.loc[hypothesis]
        diff = values[(name, "SA-Z")][0] - values[(name, "BASE-Z")][0]
        stored_diff = stored["diff_sa_minus_comparator"]
        checks.append({"hypothesis": hypothesis, "n": len(sub), "stored_n": int(stored["n"]), "diff": diff,
                       "stored_diff": float(stored_diff), "favor": favor, "stored_favor": int(stored["n_favor"]),
                       "against": against, "stored_against": int(stored["n_against"])})
    return values, checks


def main() -> None:
    metrics = pd.read_csv(S10 / "metrics.tsv", sep="\t")
    calls = pd.read_csv(S10 / "esophagus_calls.tsv", sep="\t")
    restriction = pd.read_csv(S10 / "restriction.tsv", sep="\t")
    shares = pd.read_csv(S10 / "shares.tsv", sep="\t")
    sim = pd.read_csv(ROOT / "results/stage5/sim_ext/full.tsv", sep="\t")
    rna = pd.read_csv(S9, sep="\t")
    primary = pd.read_csv(ROOT / "results/stage4/confirm/primary.tsv", sep="\t")
    contrib = pd.read_csv(ROOT / "results/stage9/imvigor_contribution_fix1.tsv", sep="\t")
    stage7 = pd.read_csv(ROOT / "results/stage7/confirm/tables/metrics.tsv", sep="\t")
    rows = []

    def add(group: str, meaning: str, value: float, expected: str, digits: int = 3, source: str = "", n=None):
        shown = show(value, digits) if digits >= 0 else str(int(round(value)))
        rows.append({"group": group, "meaning": meaning, "stored": value, "n": n, "shown": shown,
                     "expected": expected, "match": shown == expected, "source": source})

    expected = {
        ("evaluation", "top1"): {"BASE-Z": "0.724", "SA-Z": "0.455", "PURE-Zs": "0.360", "PURE-Zs-w": "0.332", "MIX-Z0": "0.776"},
        ("at_risk", "host_rate"): {"BASE-Z": "0.256", "SA-Z": "0.101", "MIX-Z0": "0.092"},
        ("native_truth", "native_truth_top1"): {"BASE-Z": "0.871", "SA-Z": "0.554", "MIX-Z0": "0.891"},
    }
    for (set_name, metric), wanted in expected.items():
        for model, text in wanted.items():
            row = one(metrics, model=model, set=set_name, cohort="all", metric=metric)
            add("microarray", f"{set_name} {metric}, {NAMES[model]}", row["value"], text,
                source="stage10 metrics.tsv", n=int(row["n"]))
    for model, text in (("BASE-Z", "72.4"), ("SA-Z", "45.5"), ("MIX-Z0", "77.6")):
        add("abstract", f"microarray top-1 %, {NAMES[model]}", 100 * one(metrics, model=model, set="evaluation", cohort="all", metric="top1")["value"], text, 1, "stage10 metrics.tsv")
    for model, text in (("BASE-Z", "25.6"), ("MIX-Z0", "9.2")):
        add("abstract", f"microarray host-attraction %, {NAMES[model]}", 100 * one(metrics, model=model, set="at_risk", cohort="all", metric="host_rate")["value"], text, 1, "stage10 metrics.tsv")
    for model, text in (("BASE-Z", "30"), ("SA-Z", "147"), ("PURE-Zs", "183"), ("PURE-Zs-w", "191"), ("MIX-Z0", "11")):
        row = one(calls, model=model, cohort="all")
        add("esophagus calls", f"evaluation set, {NAMES[model]}", row["n_esophagus"], text, -1, "stage10 esophagus_calls.tsv", int(row["n"]))
    for model, text in (("BASE-Z", "0.006"), ("MIX-Z0", "0.006"), ("PURE-Zs", "0.049"), ("PURE-Zs-w", "0.049"), ("SA-Z", "0.031")):
        part = restriction.loc[restriction["model"].eq(model)]
        if len(part) != 8:
            raise SystemExit(f"restriction rows {model}: {len(part)}")
        add("random removal", f"mean drop over 8 platform-cohort pairs, {NAMES[model]}", part["random_mean"].mean(), text, source="stage10 restriction.tsv", n=len(part))
    for contrast, name, text in (("accuracy", "s_mix", "\u22120.194"), ("accuracy", "s_std", "1.354"), ("accuracy", "s_stdw", "1.458"),
                                 ("random_removal", "s2_mix", "\u22120.028"), ("random_removal", "s2_std", "1.711"), ("random_removal", "s2_stdw", "1.704")):
        add("shares", f"{contrast} {name}", one(shares, contrast=contrast, name=name)["value"], text, source="stage10 shares.tsv")
    for method, text in (("SA-LOHO-Brain - Cortex", "0.758"), ("BASE-Z", "0.068")):
        row = one(sim, tissue="Brain - Cortex", rho="0.6", method=method)
        add("simulation", f"brain cortex, rho 0.6, {method}", row["host_rate"], text, source="stage5 sim_ext/full.tsv", n=int(row["n_at_risk"]))
    for method, text in (("MIX-Z0", "0.846"), ("PURE-Zs", "0.857"), ("PURE-Zs-w", "0.868"), ("SA-Z", "0.802"), ("BASE-Z", "0.868")):
        row = one(rna, cohort="POG570", method=method)
        add("POG570 native truth", NAMES[method], row["native_truth_top1"], text, source="stage9 ablation/metrics.tsv", n=int(row["n_native_truth"]))
    for cohort, method, text in (("POG570", "MIX-Z0", "0.034"), ("aux_rnaseq", "MIX-Z0", "0.082"), ("POG570", "SA-Z", "0.024"), ("aux_rnaseq", "SA-Z", "0.056")):
        row = one(rna, cohort=cohort, method=method)
        add("RNA-seq host attraction", f"{cohort}, {NAMES[method]}", row["host_rate"], text, source="stage9 ablation/metrics.tsv", n=int(row["n_at_risk"]))
    h3 = one(primary, hypothesis="H3")
    add("existing", "POG570 native-truth difference, points", -100 * h3["diff"], "6.6", 1, "stage4 confirm/primary.tsv", int(h3["n"]))
    h1 = one(primary, hypothesis="H1")
    add("existing", "POG570 discordant pairs in favor (H1)", h1["n_favor"], "77", -1, "stage4 confirm/primary.tsv")
    add("existing", "POG570 discordant pairs against (H1)", h1["n_against"], "0", -1, "stage4 confirm/primary.tsv")
    for method, text in (("BASE-Z", "145"), ("SA-Z", "65")):
        add("existing", f"IMvigor210 native truth called esophagus, {NAMES[method]}", contrib.loc[contrib["method"].eq(method), "n_samples"].unique()[0], text, -1, "stage9 imvigor_contribution_fix1.tsv")
    imv = one(stage7, layer="RNA-seq", cohort="blca_iatlas_imvigor210_2017", standard_site="all", method="SA-Z")
    add("existing", "IMvigor210 native-truth n", imv["n_native_truth"], "194", -1, "stage7 confirm/tables/metrics.tsv")

    values, checks = stage7_microarray()
    for (set_name, method), text in ((("evaluation", "BASE-Z"), "0.724"), (("evaluation", "SA-Z"), "0.455"),
                                     (("at_risk", "BASE-Z"), "0.256"), (("at_risk", "SA-Z"), "0.101"),
                                     (("native_truth", "BASE-Z"), "0.871"), (("native_truth", "SA-Z"), "0.554")):
        value, n = values[(set_name, method)]
        add("stage-7 microarray layer", f"{set_name}, {NAMES[method]}", value, text, source="stage7 predictions, layer sets", n=n)

    table = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT, sep="\t", index=False)
    with pd.option_context("display.width", 220, "display.max_colwidth", 70, "display.max_rows", 200):
        print(table[["group", "meaning", "n", "stored", "shown", "expected", "match"]].to_string(index=False))
        print(pd.DataFrame(checks).to_string(index=False))
    bad_checks = [c for c in checks if c["n"] != c["stored_n"] or abs(c["diff"] - c["stored_diff"]) > 1e-12
                  or c["favor"] != c["stored_favor"] or c["against"] != c["stored_against"]]
    print(f"numbers {len(table)}, mismatches {int((~table['match']).sum())}, stage-7 effect mismatches {len(bad_checks)}")
    if (~table["match"]).any() or bad_checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
