"""Post hoc recount of the Figure S3a microarray group on the microarray layer set (A18 §2.1).

The set is the one of Table 1 and results/stage10/derived/microarray_rates.tsv: microarray layer, layer test
rows, evaluation selection, one sample per patient across the layer (polyA library first, then sample id,
then cohort). Among its biopsies whose true organ is not esophagus, the fraction called esophagus is counted
for the baseline, HostMix-TOO and linear deconvolution, as in results/stage8/diagnostics/esophagus_absorption.tsv.
Nothing is refit or rescored; predictions are read from results/stage7/confirm/predictions.parquet.

Output: results/stage10/derived/s3a_microarray_layer.tsv (post hoc aggregate).
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
D17 = importlib.import_module("64_derived_facts_A17")

ROOT = D17.ROOT
OUT = D17.OUT / "s3a_microarray_layer.tsv"
METHODS = ["BASE-Z", "SA-Z", "LD-Z"]
ESOPHAGUS = "Esophagus"


def main() -> None:
    columns = ["cohort", "patient_id", "sample_id", "library", "organ", "layer", "layer_test", "selected_eval",
               *[f"{m}__pred" for m in METHODS]]
    pred = pd.read_parquet(D17.PRED, columns=columns)
    array = pred.loc[pred["layer"].eq(D17.MICROARRAY) & D17.flag(pred["layer_test"])]
    evaluation = D17.one_per_patient(array.loc[D17.flag(array["selected_eval"])])
    keep = evaluation["organ"].ne(ESOPHAGUS)
    rows = []
    for method in METHODS:
        calls = evaluation[f"{method}__pred"].eq(ESOPHAGUS)
        rows.append({"note": "post hoc", "analysis": "microarray", "cohort": "all", "set": "layer evaluation",
                     "method": method, "n_evaluation": int(len(evaluation)),
                     "n_truth_esophagus": int((~keep).sum()), "n_truth_not_esophagus": int(keep.sum()),
                     "n_esophagus_calls_all": int(calls.sum()),
                     "n_esophagus_calls_not_esophagus": int((calls & keep).sum()),
                     "esophagus_fraction": float((calls & keep).sum() / keep.sum())})
    frame = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT, sep="\t", index=False)
    print(frame.to_string(index=False))
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "rows": len(frame)}))


if __name__ == "__main__":
    main()
