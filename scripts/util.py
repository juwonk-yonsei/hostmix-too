"""Shared paths, seeds, and logging. Import only after thread env vars are set."""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import yaml

SEED = 20261001
N_THREADS = 16


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_paths() -> dict:
    root = project_root()
    with open(root / "config" / "paths.yaml") as handle:
        raw = yaml.safe_load(handle)
    out = {"root": root, "seed": int(raw.get("seed", SEED)), "n_threads": int(raw.get("n_threads", N_THREADS))}
    for key, value in raw.items():
        if key in ("seed", "n_threads"):
            continue
        out[key] = root / value
    return out


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def append_tsv(path: Path, header: str, row: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(header + "\n")
    with open(path, "a") as handle:
        handle.write(row + "\n")


def assert_disjoint(left, right, message: str) -> None:
    overlap = set(left) & set(right)
    assert not overlap, f"{message}: {len(overlap)} overlapping ids, examples {list(overlap)[:5]}"
