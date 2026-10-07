"""Gene-id crosswalk and the common protein-coding set G.

Symbol source for Toil and MET500 is the gencode v23 probemap.
A symbol is eligible only when it is also an HGNC protein-coding symbol.
POG570 ids are unversioned Ensembl ids. They are matched to the probemap
after stripping the version. If one unversioned id maps to more than one
symbol, that id is dropped and counted.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from loader import pog570_gene_ids, read_hgnc, read_probemap
from util import load_paths


def protein_coding_symbols() -> set[str]:
    hgnc = read_hgnc()
    pc = hgnc.loc[hgnc["locus_group"] == "protein-coding gene", "symbol"]
    return set(pc.dropna().astype(str))


def probemap_symbol() -> pd.DataFrame:
    """Versioned Ensembl id -> gencode symbol, restricted to HGNC protein-coding symbols."""
    symbols = protein_coding_symbols()
    probe = read_probemap()
    probe = probe.loc[probe["gene"].isin(symbols), ["id", "gene"]].drop_duplicates()
    return probe


def entrez_to_symbol() -> dict[str, str]:
    hgnc = read_hgnc()
    pc = hgnc.loc[hgnc["locus_group"] == "protein-coding gene", ["entrez_id", "symbol", "status"]]
    pc = pc.dropna(subset=["entrez_id", "symbol"])
    pc["entrez_id"] = pc["entrez_id"].astype(str).str.replace(r"\.0$", "", regex=True)
    # Prefer Approved if an entrez id is repeated.
    pc["_rank"] = (pc["status"] != "Approved").astype(int)
    pc = pc.sort_values(["entrez_id", "_rank"]).drop_duplicates("entrez_id", keep="first")
    return dict(zip(pc["entrez_id"], pc["symbol"]))


def matrix_gene_ids(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    ids = []
    with opener(path, "rt") as handle:
        handle.readline()
        for line in handle:
            gene = line.split("\t", 1)[0].strip()
            if gene:
                ids.append(gene)
    return ids


def build_common_genes(toil_ids: list[str], met_ids: list[str]) -> dict:
    probe = probemap_symbol()
    id_to_symbol = dict(zip(probe["id"], probe["gene"]))
    # Unversioned Ensembl -> symbols. Drop ids whose symbols disagree.
    by_base: dict[str, set[str]] = {}
    for ens, sym in id_to_symbol.items():
        base = ens.split(".")[0]
        by_base.setdefault(base, set()).add(sym)
    unique_base = {base: next(iter(syms)) for base, syms in by_base.items() if len(syms) == 1}
    n_base_conflict = sum(1 for syms in by_base.values() if len(syms) > 1)

    toil_symbols = {id_to_symbol[i] for i in toil_ids if i in id_to_symbol}
    met_symbols = {id_to_symbol[i] for i in met_ids if i in id_to_symbol}
    pog_ids = pog570_gene_ids()
    pog_symbols = {unique_base[i.split(".")[0]] for i in pog_ids if i.split(".")[0] in unique_base}
    pog_unmapped = sum(1 for i in pog_ids if i.split(".")[0] not in unique_base)

    common = sorted(toil_symbols & met_symbols & pog_symbols)
    return {
        "genes": common,
        "n_toil_symbols": len(toil_symbols),
        "n_met500_symbols": len(met_symbols),
        "n_pog570_symbols": len(pog_symbols),
        "n_pog570_ids": len(pog_ids),
        "n_pog570_unmapped_ids": pog_unmapped,
        "n_unversioned_symbol_conflicts": n_base_conflict,
        "n_met500_ids": len(met_ids),
        "n_met500_ids_unmapped": sum(1 for i in met_ids if i not in id_to_symbol),
        "n_toil_ids": len(toil_ids),
        "n_toil_ids_mapped": sum(1 for i in toil_ids if i in id_to_symbol),
        "id_to_symbol": id_to_symbol,
        "unique_base": unique_base,
    }


def load_gene_sets(genes: list[str]):
    """Return set names and member indices into `genes`, keeping sets with at least 5 members in G."""
    paths = load_paths()
    gmt = paths["data_reference"] / "ONCOfind/GeneSetFeaturization/GeneSetFeaturization/Hc2c6c8_merged.CGfilt.F8300.v2023.2.Hs.entrez.gmt"
    entrez = entrez_to_symbol()
    index = {gene: i for i, gene in enumerate(genes)}
    names = []
    members = []
    n_raw = 0
    with open(gmt) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            n_raw += 1
            idxs = []
            seen = set()
            for token in parts[2:]:
                sym = entrez.get(token)
                if sym is None or sym not in index or sym in seen:
                    continue
                seen.add(sym)
                idxs.append(index[sym])
            if len(idxs) >= 5:
                names.append(parts[0])
                members.append(idxs)
    return names, members, n_raw
