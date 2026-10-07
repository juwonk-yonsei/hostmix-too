"""Run the downloaded SCOPE ensemble once per cohort. No plots."""
from __future__ import annotations

import argparse
from pathlib import Path

import cancerscope
import pandas as pd
from cancerscope.scope_ensemble_functions import get_ensemble_score

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
INPUTS = ROOT / "results/stage8/external/inputs"
OUT = ROOT / "results/stage8/external/scope"
SCALE = {"POG570", "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019", "paad_iatlas_prince_2022"}


def average_probabilities(pred_dict, samples: list[str]) -> pd.DataFrame:
    totals: dict[str, dict[str, float]] = {sample: {} for sample in samples}
    counts: dict[str, dict[str, int]] = {sample: {} for sample in samples}
    for arr in pred_dict.values():
        for sample, pairs in zip(samples, arr):
            for label, value in pairs:
                totals[sample][str(label)] = totals[sample].get(str(label), 0.0) + float(value)
                counts[sample][str(label)] = counts[sample].get(str(label), 0) + 1
    labels = sorted({label for sample in totals for label in totals[sample]})
    rows = []
    for sample in samples:
        row = {"sample": sample}
        for label in labels:
            n = counts[sample].get(label, 0)
            row[label] = totals[sample].get(label, 0.0) / n if n else 0.0
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample")


def run_file(scope_obj, path: Path, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix, samples, features, code = scope_obj.load_data(str(path))
    pred_dict = scope_obj.predict(
        X=matrix, x_features=features, x_features_genecode=code,
        x_sample_names=list(samples), get_preds_dict=True,
    )
    if len(pred_dict) != 5:
        raise SystemExit(f"{path.name} models {len(pred_dict)}")
    ensemble = get_ensemble_score(pred_dict)
    ensemble["sample_name"] = [samples[int(i)] for i in ensemble["sample_ix"].tolist()]
    ensemble.to_csv(out_dir / "SCOPE_topPredictions.txt", sep="\t", index=False)
    averaged = average_probabilities(pred_dict, list(samples))
    averaged.to_csv(out_dir / "SCOPE_meanProbabilities.tsv", sep="\t", float_format="%.8e")
    print(f"wrote {out_dir} samples {len(samples)} models {len(pred_dict)}", flush=True)
    return ensemble


def doubled(path: Path, dest: Path) -> None:
    frame = pd.read_csv(path, sep="\t")
    value_cols = frame.columns[1:]
    frame[value_cols] = frame[value_cols] * 2
    frame.to_csv(dest, sep="\t", index=False, float_format="%.8e")


def top1(path: Path) -> pd.Series:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    top = frame.loc[frame["rank_pred"].astype(str) == "1"]
    if top["sample_name"].duplicated().any():
        raise SystemExit(f"duplicate rank 1 in {path}")
    return top.set_index("sample_name")["label"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", action="append")
    args = parser.parse_args()
    names = args.cohort or [path.name.split(".scope")[0] for path in sorted(INPUTS.glob("*.scope.tsv.gz"))]
    scope_obj = cancerscope.scope()
    changed = []
    for name in names:
        src = INPUTS / f"{name}.scope.tsv.gz"
        primary = run_file(scope_obj, src, OUT / name)
        if name not in SCALE:
            continue
        temp = OUT / f"{name}.x2.tsv"
        doubled(src, temp)
        run_file(scope_obj, temp, OUT / f"{name}_x2")
        temp.unlink()
        a = top1(OUT / name / "SCOPE_topPredictions.txt")
        b = top1(OUT / f"{name}_x2" / "SCOPE_topPredictions.txt")
        both = a.index.intersection(b.index)
        fraction = float((a.loc[both] != b.loc[both]).mean()) if len(both) else float("nan")
        changed.append({"cohort": name, "n": int(len(both)), "top1_changed_fraction": fraction})
        print(f"scale {name} changed {fraction}", flush=True)
    if changed:
        pd.DataFrame(changed).to_csv(OUT / "scale_check.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
