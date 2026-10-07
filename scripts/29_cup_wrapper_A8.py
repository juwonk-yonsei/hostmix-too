"""CUP-AI-Dx wrapper. Same preprocessing as external_test.py. Does not edit the tool."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CUP_ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too/results/stage7/external_tools/CUP-AI-Dx")
sys.path.insert(0, str(CUP_ROOT))

from datasets import load_dataset  # noqa: E402
from keras.models import load_model  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from utils import reshape_data_1d  # noqa: E402
import json  # noqa: E402

TARGET = "92.13197969543148"


def encodings(data_dir: Path):
    raw = json.loads((data_dir / "label_encoding.json").read_text())
    int_to_label = {int(key): value for key, value in raw.items()}
    return int_to_label


def scale_like_external_test(frame: pd.DataFrame) -> np.ndarray:
    features = pd.read_csv(CUP_ROOT / "data/features_791.csv", header=None, index_col=0).index
    missing = [column for column in features if column not in frame.columns]
    if missing:
        raise SystemExit(f"input is missing {len(missing)} CUP features")
    ordered = frame.loc[:, features]
    scaler = StandardScaler()
    scaler_input = ordered.transpose().to_numpy(dtype=np.float64)
    scaler.fit(scaler_input)
    return np.transpose(scaler.transform(scaler_input))


def predict(matrix: np.ndarray, model):
    probabilities = model.predict(reshape_data_1d(matrix))
    labels = encodings(CUP_ROOT / "data")
    order = [labels[index] for index in range(len(labels))]
    pred = [labels[int(index)] for index in np.argmax(probabilities, axis=1)]
    return pred, probabilities, order


def check_example(model) -> None:
    matrix, truth = load_dataset(str(CUP_ROOT / "data"), "metastatic")
    pred, _probabilities, _order = predict(matrix, model)
    # external_test adjusts labels only for PDX. Metastatic keeps the raw class.
    correct = int(np.sum(np.asarray(pred, dtype=object) == truth.to_numpy(dtype=object)))
    accuracy = correct / len(truth) * 100
    printed = f"{accuracy}"
    print(f"Overall metastatic accuracy: {printed}%")
    print(f"correct {correct} n {len(truth)}")
    if printed != TARGET:
        raise SystemExit(f"wrapper accuracy {printed} != {TARGET}")


def predict_file(path: Path, out: Path, model) -> None:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    matrix = scale_like_external_test(frame)
    pred, probabilities, order = predict(matrix, model)
    table = pd.DataFrame(probabilities, index=frame.index, columns=order)
    table.insert(0, "pred", pred)
    table.index.name = "sample"
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, sep="\t", float_format="%.8e")
    print(f"wrote {out} n {len(table)}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-example", action="store_true")
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    model = load_model(CUP_ROOT / "models/inception_net_1d.h5")
    if args.check_example:
        check_example(model)
    if args.input:
        predict_file(Path(args.input), Path(args.output), model)
    if args.input_dir:
        out_dir = Path(args.output_dir)
        for path in sorted(Path(args.input_dir).glob("*.cup.tsv.gz")):
            predict_file(path, out_dir / f"{path.name.split('.cup')[0]}.tsv", model)


if __name__ == "__main__":
    main()
