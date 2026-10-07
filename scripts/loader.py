"""Common loader. Stage-1 scripts read cohort files through this module.

POG570 expression is locked. Requesting the matrix without unlock=True raises.
Allowed without unlock: structure summary, gene-id list, metadata frequencies,
and the biopsy-site by diagnosis table.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from util import load_paths, sha256_file

_PATHS = load_paths()


class POG570Locked(RuntimeError):
    pass


def _raw(*parts: str) -> Path:
    return _PATHS["data_raw"].joinpath(*parts)


def read_phenotype() -> pd.DataFrame:
    path = _raw("toil", "TcgaTargetGTEX_phenotype.txt.gz")
    # The file is not valid UTF-8 (byte 0xCA). latin-1 keeps every byte.
    return pd.read_csv(path, sep="\t", dtype=str, encoding="latin-1")


def read_category() -> pd.DataFrame:
    return pd.read_csv(_raw("toil", "TCGA_GTEX_category.txt"), sep="\t", dtype=str)


def read_probemap() -> pd.DataFrame:
    return pd.read_csv(_raw("toil", "gencode.v23.annotation.gene.probemap"), sep="\t", dtype=str)


def read_hgnc() -> pd.DataFrame:
    return pd.read_csv(_raw("hgnc", "hgnc_complete_set.txt"), sep="\t", dtype=str, low_memory=False)


def read_met500_meta() -> pd.DataFrame:
    return pd.read_csv(_raw("met500", "M.meta.plus.txt"), sep="\t", dtype=str)


def read_skcm_clinical() -> pd.DataFrame:
    return pd.read_csv(_raw("xena", "SKCM_clinicalMatrix"), sep="\t", dtype=str)


def read_pog570_table(which: str) -> pd.DataFrame:
    name = {"s1": "Table_S1_Demographics.xlsx", "s2": "Table_S2_Treatment.xlsx"}[which]
    # dtype=str keeps trailing spaces that are part of the raw labels.
    return pd.read_excel(_raw("POG570", name), dtype=str)


def pog570_gene_ids() -> list[str]:
    """Gene identifiers only. Expression values are not returned."""
    path = _raw("POG570", "POG570_TPM_expression.txt.gz")
    ids = []
    with gzip.open(path, "rt") as handle:
        header = handle.readline()
        if not header:
            return ids
        for line in handle:
            gene = line.split("\t", 1)[0].strip()
            if gene:
                ids.append(gene)
    return ids


def pog570_expression_structure() -> dict:
    """Sample count, gene count, id format, and exact min/median/max.

    Values are used only for this scale summary, then deleted. They are not
    returned, saved, or joined to diagnosis or biopsy site.
    """
    path = _raw("POG570", "POG570_TPM_expression.txt.gz")
    parts = []
    n_genes = 0
    example = None
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        for line in handle:
            gene, _, rest = line.partition("\t")
            if example is None:
                example = gene
            parts.append(np.fromstring(rest, sep="\t", dtype=np.float32))
            n_genes += 1
    values = np.concatenate(parts) if parts else np.array([], dtype=np.float32)
    del parts
    finite = values[np.isfinite(values)]
    summary = {
        "n_samples": len(sample_ids),
        "n_genes": n_genes,
        "gene_id_example": example,
        "gene_id_format": "Ensembl gene id without version",
        "sample_id_example": sample_ids[0] if sample_ids else None,
        "n_sample_id_lengths": sorted({len(x) for x in sample_ids}),
        "raw_min": float(finite.min()) if finite.size else None,
        "raw_max": float(finite.max()) if finite.size else None,
        "raw_median": float(np.median(finite)) if finite.size else None,
        "n_negative_values": int(np.sum(values < 0)),
        "n_nonfinite": int(values.size - finite.size),
        "n_values": int(values.size),
        "readme_unit": "TPM",
    }
    del values, finite
    return summary


def load_pog570_expression(unlock: bool = False):
    """Return the expression matrix only when unlock=True. Stage 1 must not call this."""
    if not unlock:
        raise POG570Locked(
            "POG570 expression is locked until stage 3. "
            "Refusing to return the expression matrix without unlock=True."
        )
    path = _raw("POG570", "POG570_TPM_expression.txt.gz")
    return pd.read_csv(path, sep="\t", index_col=0)


def write_pog570_lock_hashes() -> Path:
    """Record SHA-256 of every file in data/raw/POG570."""
    folder = _raw("POG570")
    out = _PATHS["config"] / "locked_hashes.tsv"
    rows = ["path\tsize_bytes\tsha256"]
    for path in sorted(folder.iterdir()):
        if not path.is_file():
            continue
        rows.append(f"data/raw/POG570/{path.name}\t{path.stat().st_size}\t{sha256_file(path)}")
    out.write_text("\n".join(rows) + "\n")
    return out


def verify_pog570_lock_hashes() -> None:
    table = pd.read_csv(_PATHS["config"] / "locked_hashes.tsv", sep="\t", dtype=str)
    for row in table.itertuples(index=False):
        path = _PATHS["root"] / row.path
        digest = sha256_file(path)
        if digest != row.sha256:
            raise POG570Locked(f"POG570 file hash changed: {row.path}")
