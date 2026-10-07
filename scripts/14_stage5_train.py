#!/usr/bin/env python3
"""Train the stage-5 SA variants and write the shared simulation mixes.

Does not set a thread cap. Does not run scripts/12_confirm.py or scripts/13_posthoc.py.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

from stage2_fit import fit_logit
from stage5_lib import (
    POOL10,
    POOL12,
    SIM_RHOS,
    STAGE2_SIM_TISSUES,
    c_for,
    draw_mixes,
    dump_json,
    gtex_ref_ids,
    gtex_sim_ids,
    holdout_log_loss,
    load_gene_pack,
    load_paths,
    load_tpm,
    materialize,
    organ_of_labels,
    predict_organ,
    rank_b0,
    simulation_tumor_index,
    stage5_dir,
    sum_1e6,
    test_frame,
    train_frame,
    variant_names,
    variant_seeds,
    write_disk,
)
from util import SEED, sha256_file

LABELS = {}


def register_labels(y: np.ndarray) -> None:
    LABELS["y"] = y


def train_labels() -> np.ndarray:
    if "y" not in LABELS:
        raise SystemExit("training labels were not registered")
    return LABELS["y"]


def records_equal(left: pd.DataFrame, right: list[dict]) -> dict:
    other = pd.DataFrame(right)
    out = {"n_left": int(len(left)), "n_right": int(len(other))}
    if len(left) != len(other):
        out["field_match"] = False
        return out
    for column in ("tumor", "host", "tissue"):
        out[f"{column}_mismatch"] = int((left[column].astype(str).to_numpy() != other[column].astype(str).to_numpy()).sum())
    rho_diff = np.abs(left["rho"].to_numpy(dtype=np.float64) - other["rho"].to_numpy(dtype=np.float64))
    out["rho_max_abs"] = float(rho_diff.max()) if len(rho_diff) else None
    out["field_match"] = (
        all(out[f"{column}_mismatch"] == 0 for column in ("tumor", "host", "tissue"))
        and out["rho_max_abs"] == 0.0
    )
    return out


def stack_pure_mix(pure: np.ndarray, mixes: np.ndarray, n_mix: int) -> np.ndarray:
    n, p = pure.shape
    stride = 1 + n_mix
    out = np.empty((n * stride, p), dtype=np.float32)
    out[0::stride] = pure
    for j in range(n_mix):
        out[1 + j::stride] = mixes[j::n_mix]
    return out


def y_repeat(n_each: int) -> np.ndarray:
    return np.repeat(train_labels(), n_each)


def loho_records(base_records, left_out, ids, rng):
    remaining = [tissue for tissue in POOL10 if tissue != left_out]
    out = []
    for row in base_records:
        if row["tissue"] != left_out:
            out.append(dict(row))
            continue
        tissue = remaining[int(rng.integers(0, len(remaining)))]
        group = ids[tissue]
        host = group[int(rng.integers(0, len(group)))]
        out.append({"tumor": row["tumor"], "host": host, "tissue": tissue, "rho": float(row["rho"])})
    return out


def redraw_rho(base_records, low, rng):
    return [{**row, "rho": float(rng.uniform(low, 1.0))} for row in base_records]


def fit_mlp(x_raw, y, patients, seed: int, loss_path: Path) -> dict:
    labels = np.array(sorted(set(y.tolist())))
    frame = pd.DataFrame({"patient": patients, "label": y}).drop_duplicates("patient")
    rest, hold = train_test_split(
        frame["patient"].to_numpy(), test_size=0.1, stratify=frame["label"].to_numpy(), random_state=seed,
    )
    rest_set, hold_set = set(rest.tolist()), set(hold.tolist())
    if rest_set & hold_set:
        raise SystemExit("MLP patient split overlaps")
    hold_mask = np.array([patient in hold_set for patient in patients])
    rest_mask = ~hold_mask
    scaler_search = StandardScaler().fit(x_raw[rest_mask])
    x_rest = scaler_search.transform(x_raw[rest_mask])
    x_hold = scaler_search.transform(x_raw[hold_mask])
    y_rest, y_hold = y[rest_mask], y[hold_mask]
    clf = new_mlp(seed)
    rng = np.random.default_rng(seed)
    losses = []
    best = np.inf
    best_epoch = 0
    wait = 0
    stopped = False
    label_list = labels.tolist()
    for epoch in range(1, 101):
        order = rng.permutation(len(y_rest))
        clf.partial_fit(x_rest[order], y_rest[order], classes=labels)
        loss = holdout_log_loss(clf, x_hold, y_hold, label_list)
        losses.append({"epoch": epoch, "holdout_log_loss": loss})
        print("mlp", loss_path.stem, "epoch", epoch, "loss", loss, flush=True)
        if loss < best:
            best = loss
            best_epoch = epoch
            wait = 0
        else:
            wait += 1
            if wait >= 10:
                stopped = True
                break
    pd.DataFrame(losses).to_csv(loss_path, sep="\t", index=False)
    scaler = StandardScaler().fit(x_raw)
    x_all = scaler.transform(x_raw)
    final = new_mlp(seed)
    rng = np.random.default_rng(seed)
    for _epoch in range(best_epoch):
        order = rng.permutation(len(y))
        final.partial_fit(x_all[order], y[order], classes=labels)
    return {
        "clf": final,
        "scaler": scaler,
        "C": None,
        "n_iter": int(best_epoch),
        "max_iter_hit": bool((not stopped) and len(losses) == 100),
        "classes": [str(c) for c in final.classes_],
        "E_star": int(best_epoch),
        "search_epochs": int(len(losses)),
        "patience_stopped": bool(stopped),
        "n_holdout_patients": int(len(hold_set)),
        "n_fit_patients": int(len(rest_set)),
    }


def new_mlp(seed: int) -> MLPClassifier:
    return MLPClassifier(
        hidden_layer_sizes=(512,), activation="relu", solver="adam", alpha=1e-4,
        batch_size=256, learning_rate_init=1e-3, random_state=seed, shuffle=False,
        early_stopping=False, max_iter=1,
    )


def design_matrix(name, z_pure, z_sa, base_records, tumors, train_pos, tumor_norm, bank, ids, b0_idx, seed, z_cache):
    if name in ("SA-m0", "BASE-MLP"):
        return z_pure, y_repeat(1), 1
    if name == "SA-m1":
        x = np.vstack([z_pure, z_sa[0::4]])
        labels = train_labels()
        return x, np.concatenate([labels, labels]), 2
    if name in ("SA-C0.01", "SA-C0.1", "SA-MLP"):
        return stack_pure_mix(z_pure, z_sa, 4), y_repeat(5), 5
    rng = np.random.default_rng(seed)
    if name == "SA-m8":
        extra = draw_mixes(tumors, POOL10, ids, rng, False, None)
        z_new = rank_b0(materialize(tumor_norm, train_pos, bank, extra), b0_idx)
        n = len(tumors)
        mixes = np.empty((n * 8, z_pure.shape[1]), dtype=np.float32)
        for j in range(4):
            mixes[j::8] = z_sa[j::4]
            mixes[4 + j::8] = z_new[j::4]
        return stack_pure_mix(z_pure, mixes, 8), y_repeat(9), 9
    if name == "SA-r30":
        records = redraw_rho(base_records, 0.30, rng)
    elif name == "SA-r50":
        records = redraw_rho(base_records, 0.50, rng)
    elif name == "SA-pool3":
        records = draw_mixes(tumors, ["Liver", "Lung", "Spleen"], ids, rng, False, None)
    elif name == "SA-pool22":
        records = draw_mixes(tumors, POOL10 + POOL12, ids, rng, False, None)
    elif name.startswith("SA-LOHO-"):
        records = loho_records(base_records, name[len("SA-LOHO-"):], ids, rng)
    else:
        raise SystemExit(f"unknown variant {name}")
    if name not in z_cache:
        z_cache[name] = rank_b0(materialize(tumor_norm, train_pos, bank, records), b0_idx)
    return stack_pure_mix(z_pure, z_cache[name], 4), y_repeat(5), 5


def write_sim_mixes(paths, root: Path) -> pd.DataFrame:
    test = test_frame(paths)
    take = simulation_tumor_index(test)
    tumor_ids = test["sample"].to_numpy()[take]
    tissues = POOL10 + POOL12
    sim_ids = gtex_sim_ids(paths, tissues)
    dropped = [{"tissue": tissue, "n_gtex_sim": int(len(sim_ids[tissue]))} for tissue in tissues if len(sim_ids[tissue]) < 3]
    pd.DataFrame(dropped).to_csv(root / "sim_ext" / "dropped_sim_tissues.tsv", sep="\t", index=False)
    rows = []
    rng = np.random.default_rng(SEED)
    for tissue in STAGE2_SIM_TISSUES:
        ids = sim_ids[tissue]
        if len(ids) < 3:
            raise SystemExit(f"stage-2 simulation tissue has fewer than 3 GTEx-sim samples: {tissue}")
        for rho in SIM_RHOS:
            picked = rng.integers(0, len(ids), size=len(tumor_ids))
            for i, host_i in enumerate(picked):
                rows.append({
                    "tissue": tissue, "rho": rho, "tumor": tumor_ids[i], "host": ids[int(host_i)],
                    "draw": "stage2_seed_20261001",
                })
    rest = [tissue for tissue in tissues if tissue not in STAGE2_SIM_TISSUES and len(sim_ids[tissue]) >= 3]
    rng_new = np.random.default_rng(20261005)
    for tissue in rest:
        ids = sim_ids[tissue]
        for rho in SIM_RHOS:
            picked = rng_new.integers(0, len(ids), size=len(tumor_ids))
            for i, host_i in enumerate(picked):
                rows.append({
                    "tissue": tissue, "rho": rho, "tumor": tumor_ids[i], "host": ids[int(host_i)],
                    "draw": "stage5_seed_20261005",
                })
    frame = pd.DataFrame(rows)
    frame.to_parquet(root / "sim_ext" / "mixes.parquet", index=False)
    pd.DataFrame({"tumor": tumor_ids}).to_csv(root / "sim_ext" / "tumors.tsv", sep="\t", index=False)
    return frame


def verify_stage2_sim(paths, root: Path, mixes: pd.DataFrame, b0_idx) -> None:
    test = test_frame(paths)
    take = simulation_tumor_index(test)
    tumor_ids = test["sample"].to_numpy()[take]
    truth = organ_of_labels(test["learning_label"].to_numpy()[take])
    needed = set(tumor_ids.tolist())
    part = mixes.loc[mixes["draw"] == "stage2_seed_20261001"]
    needed.update(part["host"].tolist())
    print("verify sim load", len(needed), flush=True)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    needed_list = list(needed)
    stacked = sum_1e6(tpm[[row[sample] for sample in needed_list]])
    block = {sample: stacked[i] for i, sample in enumerate(needed_list)}
    del tpm, stacked
    bundles = {
        "BASE": joblib.load(paths["results"] / "stage2" / "models" / "BASE_Z.joblib"),
        "SA": joblib.load(paths["results"] / "stage2" / "models" / "SA_Z.joblib"),
    }
    stored = pd.read_csv(paths["results"] / "stage2" / "simulation.tsv", sep="\t")
    stored = stored.loc[stored["representation"] == "Z"]
    rows = []
    tumor_block = np.stack([block[sample] for sample in tumor_ids])
    for tissue in STAGE2_SIM_TISSUES:
        for rho in SIM_RHOS:
            sub = part.loc[(part["tissue"] == tissue) & (np.isclose(part["rho"], rho))]
            if sub["tumor"].tolist() != tumor_ids.tolist():
                raise SystemExit(f"sim mix tumor order mismatch {tissue} {rho}")
            hosts = np.stack([block[sample] for sample in sub["host"]])
            rho_v = np.float32(rho)
            mixed = rho_v * tumor_block + (np.float32(1.0) - rho_v) * hosts
            x = rank_b0(mixed, b0_idx)
            for method, bundle in bundles.items():
                _organ, pred = predict_organ(bundle, x)
                top1 = float(np.mean(pred == truth))
                hit = stored.loc[
                    (stored["tissue"] == tissue) & np.isclose(stored["rho"], rho) & (stored["method"] == method)
                ]
                expected = float(hit["top1"].iloc[0])
                rows.append({
                    "tissue": tissue, "rho": rho, "method": method,
                    "top1": top1, "stored_top1": expected, "abs_diff": abs(top1 - expected),
                })
                print("sim check", tissue, rho, method, top1, expected, flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(root / "sim_ext" / "stage2_reproduction.tsv", sep="\t", index=False)
    max_diff = float(table["abs_diff"].max())
    if max_diff != 0.0:
        raise SystemExit(f"stage-2 simulation mixes were not reproduced; max abs top1 diff {max_diff}")


def main() -> None:
    paths = load_paths()
    root = stage5_dir(paths)
    write_disk(paths["logs"] / "disk_stage5_start.txt")
    names = variant_names()
    seeds = variant_seeds(names)
    dump_json(root / "variants" / "variant_seeds.json", {
        "generator": "numpy Generator(20261005).integers(0, 2**31-1) in sorted-name order",
        "names_sorted": names,
        "seeds": seeds,
        "C_rule": "variants that change the mix count keep C * rows_per_tumor = 0.15",
    })
    aug = json.loads((paths["results"] / "stage2" / "augmentation_seeds.json").read_text())
    if list(aug["host_pool_sorted"]) != POOL10:
        raise SystemExit("stored host pool is not the stage-5 pool of 10")
    _genes, b0_idx, _off, _ind, _set_names = load_gene_pack(paths)
    print("simulation mixes", flush=True)
    mixes = write_sim_mixes(paths, root)
    verify_stage2_sim(paths, root, mixes, b0_idx)

    train = train_frame(paths)
    tumors = train["sample"].tolist()
    train_pos = {sample: i for i, sample in enumerate(tumors)}
    y = train["learning_label"].to_numpy()
    patients = train["patient"].to_numpy()
    register_labels(y)
    print("load tpm", flush=True)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    tumor_norm = sum_1e6(tpm[[row[sample] for sample in tumors]])
    tissues = POOL10 + POOL12
    ids = gtex_ref_ids(paths, tissues)
    dropped = [{"tissue": tissue, "n_gtex_ref": int(len(ids[tissue]))} for tissue in tissues if len(ids[tissue]) < 3]
    pd.DataFrame(dropped).to_csv(root / "variants" / "dropped_ref_tissues.tsv", sep="\t", index=False)
    if dropped:
        raise SystemExit(f"GTEx-ref tissues below 3 samples: {dropped}")
    bank = {}
    host_norm = {}
    for tissue in tissues:
        host_norm[tissue] = sum_1e6(tpm[[row[sample] for sample in ids[tissue]]])
        for i, sample in enumerate(ids[tissue]):
            bank[sample] = host_norm[tissue][i]
    del tpm
    print("pure Z", flush=True)
    z_pure = rank_b0(tumor_norm, b0_idx)
    np.save(root / "variants" / "train_Z_b0.npy", z_pure)
    np.save(root / "variants" / "train_mean_Z.npy", z_pure.mean(axis=0).astype(np.float32))
    pd.DataFrame({"sample": tumors, "patient": patients, "learning_label": y}).to_csv(
        root / "variants" / "train_rows.tsv", sep="\t", index=False
    )
    means = np.vstack([rank_b0(host_norm[tissue], b0_idx).mean(axis=0) for tissue in tissues]).astype(np.float32)
    np.save(root / "variants" / "gtex_ref_mean_Z.npy", means)
    (root / "variants" / "gtex_ref_mean_tissues.json").write_text(json.dumps(tissues))

    print("reproduce SA", flush=True)
    saved = pd.read_parquet(paths["results"] / "stage2" / "mixes_SA.parquet")
    redraw = draw_mixes(tumors, POOL10, ids, np.random.default_rng(int(aug["seeds"]["SA"])), False, None)
    rec = records_equal(saved, redraw)
    base_records = saved.to_dict("records")
    seen = []
    for item in base_records:
        if not seen or seen[-1] != item["tumor"]:
            seen.append(item["tumor"])
    if seen != tumors:
        raise SystemExit("SA mix tumor order is not the train order")
    z_saved = rank_b0(materialize(tumor_norm, train_pos, bank, base_records), b0_idx)
    z_redraw = rank_b0(materialize(tumor_norm, train_pos, bank, redraw), b0_idx)
    rec["feature_max_abs"] = float(np.max(np.abs(z_saved.astype(np.float64) - z_redraw.astype(np.float64))))
    rec["n_mix"] = int(len(saved))
    dump_json(root / "variants" / "sa_reproduce.json", rec)
    print("SA reproduce", rec, flush=True)
    if not rec["field_match"] or rec["feature_max_abs"] != 0.0:
        raise SystemExit("SA training rows were not reproduced")
    del z_redraw, redraw

    model_dir = root / "variants" / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    log_path = root / "variants" / "train_log.tsv"
    rows = pd.read_csv(log_path, sep="\t").to_dict("records") if log_path.exists() else []
    done = {row["model"] for row in rows}
    z_cache = {}
    for name in names:
        dest = model_dir / f"{name}.joblib"
        if dest.exists() and name in done:
            print("skip", name, flush=True)
            continue
        print("fit", name, flush=True)
        t0 = time.perf_counter()
        x, y_fit, n_per = design_matrix(
            name, z_pure, z_saved, base_records, tumors, train_pos, tumor_norm, bank, ids, b0_idx, seeds[name], z_cache,
        )
        if len(y_fit) != len(x):
            raise SystemExit(f"{name} label length {len(y_fit)} != {len(x)}")
        if name in ("BASE-MLP", "SA-MLP"):
            bundle = fit_mlp(x, y_fit, np.repeat(patients, n_per), seeds[name], root / "variants" / f"mlp_loss_{name}.tsv")
        else:
            bundle = fit_logit(x, y_fit, c_for(name), scale=True)
            bundle["E_star"] = None
            bundle["search_epochs"] = None
            bundle["patience_stopped"] = None
        elapsed = time.perf_counter() - t0
        rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        joblib.dump(bundle, dest)
        row = {
            "model": name, "n_rows": int(x.shape[0]), "n_per_tumor": int(n_per),
            "C": None if bundle["C"] is None else repr(float(bundle["C"])),
            "seconds": elapsed, "max_rss_kbytes_after": rss, "n_iter": bundle["n_iter"],
            "max_iter_hit": bundle["max_iter_hit"], "E_star": bundle.get("E_star"),
            "search_epochs": bundle.get("search_epochs"), "patience_stopped": bundle.get("patience_stopped"),
            "sha256": sha256_file(dest), "seed": seeds[name],
        }
        rows = [item for item in rows if item.get("model") != name]
        rows.append(row)
        pd.DataFrame(rows).to_csv(log_path, sep="\t", index=False)
        print("saved", name, "rows", x.shape[0], "sec", round(elapsed, 1), "rss", rss, flush=True)
        del x
    print("train done", flush=True)


if __name__ == "__main__":
    main()
