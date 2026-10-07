#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gzip
import re
import subprocess
import time

import requests

from util import append_tsv, load_paths, now_iso, sha256_file

# URLs that responded. .gz was tried first for category and probemap and returned 403,
# so the uncompressed objects are the ones stored.
FILES = [
    ("toil", "data/raw/toil/TcgaTargetGtex_rsem_gene_tpm.gz", "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz"),
    ("toil", "data/raw/toil/TcgaTargetGTEX_phenotype.txt.gz", "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"),
    ("toil", "data/raw/toil/TCGA_GTEX_category.txt", "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_GTEX_category.txt"),
    ("toil", "data/raw/toil/gencode.v23.annotation.gene.probemap", "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap/gencode.v23.annotation.gene.probemap"),
    ("met500", "data/raw/met500/M.mx.txt.gz", "https://ucscpublic.xenahubs.net/download/MET500/geneExpression/M.mx.txt.gz"),
    ("met500", "data/raw/met500/M.mx.log2.txt.gz", "https://ucscpublic.xenahubs.net/download/MET500/geneExpression/M.mx.log2.txt.gz"),
    ("met500", "data/raw/met500/M.meta.plus.txt", "https://ucscpublic.xenahubs.net/download/MET500/geneExpression/M.meta.plus.txt"),
    ("POG570", "data/raw/POG570/POG570_TPM_expression.txt.gz", "https://www.bcgsc.ca/downloads/POG570/POG570_TPM_expression.txt.gz"),
    ("POG570", "data/raw/POG570/Table_S1_Demographics.xlsx", "https://www.bcgsc.ca/downloads/POG570/Table_S1_Demographics.xlsx"),
    ("POG570", "data/raw/POG570/Table_S2_Treatment.xlsx", "https://www.bcgsc.ca/downloads/POG570/Table_S2_Treatment.xlsx"),
    ("POG570", "data/raw/POG570/POG570_README.txt", "https://www.bcgsc.ca/downloads/POG570/POG570_README.txt"),
    ("POG570", "data/raw/pog570_data_release_policy.html", "https://bcgsc.ca/data-release-policy-open-access"),
    ("hgnc", "data/raw/hgnc/hgnc_complete_set.txt", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"),
    ("hpa", "data/raw/hpa/rna_tissue_consensus.tsv.zip", "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"),
    ("hpa", "data/raw/hpa/rna_tissue_hpa.tsv.zip", "https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip"),
    ("hpa", "data/raw/hpa/rna_tissue_hpa_samples.tsv.zip", "https://www.proteinatlas.org/download/tsv/rna_tissue_hpa_samples.tsv.zip"),
    ("xena", "data/raw/xena/SKCM_clinicalMatrix", "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.SKCM.sampleMap/SKCM_clinicalMatrix"),
]

GEO = {
    "GSE41258": "GSE41nnn",
    "GSE14018": "GSE14nnn",
    "GSE74685": "GSE74nnn",
    "GSE50760": "GSE50nnn",
}


def write_environment(root: Path) -> None:
    dest = root / "logs/environment.txt"
    if dest.exists() and dest.stat().st_size > 0:
        return
    py = root / ".venv/bin/python"
    pip = root / ".venv/bin/pip"
    parts = [now_iso(), subprocess.check_output([str(py), "-V"], text=True), subprocess.check_output([str(pip), "freeze"], text=True)]
    lscpu = subprocess.check_output(["lscpu"], text=True)
    keep = []
    for line in lscpu.splitlines():
        if line.startswith(("Architecture", "CPU(s)", "Thread(s)", "Core(s)", "Socket(s)", "Model name", "NUMA")):
            keep.append(line)
    parts.append("==== lscpu ====\n" + "\n".join(keep) + "\n")
    parts.append("==== free -g ====\n" + subprocess.check_output(["free", "-g"], text=True))
    dest.write_text("\n".join(parts))


def log_file(root: Path, source: str, rel: str, url: str) -> None:
    path = root / rel
    log = root / "logs/download_log.tsv"
    header = "source\tfile\turl\tsize_bytes\tsha256\tdatetime\tstatus"
    if not path.exists() or path.stat().st_size == 0:
        append_tsv(log, header, f"{source}\t{rel}\t{url}\t\t\t{now_iso()}\tmissing")
        return
    digest = sha256_file(path)
    append_tsv(log, header, f"{source}\t{rel}\t{url}\t{path.stat().st_size}\t{digest}\t{now_iso()}\tok")


def stream_geo_header(url: str, dest: Path, deadline: float) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "skipped_exists"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        response.raw.decode_content = False
        decoder = gzip.GzipFile(fileobj=response.raw)
        with open(dest, "w") as handle:
            while True:
                if time.time() > deadline:
                    return "stopped_time_limit"
                line = decoder.readline()
                if not line:
                    break
                text = line.decode("utf-8", "replace")
                if text.startswith("!series_matrix_table_begin"):
                    break
                handle.write(text)
    return "ok_header_only"


def geo_summarize(root: Path) -> None:
    """Count characteristic keys from the saved series-matrix headers. Expression is not read."""
    import pandas as pd
    from collections import Counter

    dest = root / "results/audit/geo_metadata_summary.tsv"
    rows = []
    for path in sorted((root / "data/raw/geo").glob("*_header.txt")):
        fields: dict[str, list[str]] = {}
        with open(path) as handle:
            for line in handle:
                if not line.startswith("!Sample_"):
                    continue
                parts = line.rstrip("\n").split("\t")
                fields.setdefault(parts[0], []).extend(part.strip().strip('"') for part in parts[1:])
        n_samples = len(fields.get("!Sample_geo_accession", []))
        platforms = sorted(set(fields.get("!Sample_platform_id", [])))
        keyed: dict[str, Counter] = {}
        for name, values in fields.items():
            if "characteristics" not in name and name not in ("!Sample_source_name_ch1", "!Sample_title"):
                continue
            for value in values:
                if ":" in value:
                    key, rest = value.split(":", 1)
                    key = key.strip()
                    rest = rest.strip()
                else:
                    key, rest = name, value
                keyed.setdefault(key, Counter())[rest] += 1
        accession = path.name.split("__", 1)[0]
        for key, counter in keyed.items():
            for value, count in counter.most_common():
                rows.append({
                    "accession": accession,
                    "n_samples": n_samples,
                    "platform": ",".join(platforms),
                    "field": key,
                    "value": value,
                    "count": count,
                })
    pd.DataFrame(rows).to_csv(dest, sep="\t", index=False)


def geo_metadata(root: Path) -> None:
    deadline = time.time() + 60 * 60
    out_dir = root / "data/raw/geo"
    out_dir.mkdir(parents=True, exist_ok=True)
    log = root / "logs/download_log.tsv"
    header = "source\tfile\turl\tsize_bytes\tsha256\tdatetime\tstatus"
    for gse, folder in GEO.items():
        if time.time() > deadline:
            break
        index_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{folder}/{gse}/matrix/"
        page = requests.get(index_url, timeout=60)
        page.raise_for_status()
        names = sorted(set(re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', page.text)))
        for name in names:
            url = index_url + name
            dest = out_dir / f"{gse}__{name.replace('.txt.gz','')}_header.txt"
            status = stream_geo_header(url, dest, deadline)
            size = dest.stat().st_size if dest.exists() else 0
            digest = sha256_file(dest) if size else ""
            append_tsv(log, header, f"geo\t{dest.relative_to(root)}\t{url}\t{size}\t{digest}\t{now_iso()}\t{status}")
            print(gse, name, status, size, flush=True)


def main() -> None:
    paths = load_paths()
    root = paths["root"]
    (root / "logs").mkdir(exist_ok=True)
    if not (root / "logs/disk_start.txt").exists():
        (root / "logs/disk_start.txt").write_text(subprocess.check_output(["df", "-h", str(root)], text=True) + now_iso() + "\n")
    write_environment(root)
    # Rewrite the download log from current files so reruns do not duplicate rows.
    log = root / "logs/download_log.tsv"
    if log.exists():
        log.unlink()
    for source, rel, url in FILES:
        print("hash", rel, flush=True)
        log_file(root, source, rel, url)
    onc = root / "data/reference/ONCOfind"
    if (onc / ".git").exists():
        commit = subprocess.check_output(["git", "-C", str(onc), "rev-parse", "HEAD"], text=True).strip()
        append_tsv(
            log,
            "source\tfile\turl\tsize_bytes\tsha256\tdatetime\tstatus",
            f"oncfind\tdata/reference/ONCOfind\thttps://github.com/yeonuk-Jeong/ONCOfind\t\t{commit}\t{now_iso()}\tgit_clone",
        )
    geo_metadata(root)
    print("download log written")


if __name__ == "__main__":
    main()
