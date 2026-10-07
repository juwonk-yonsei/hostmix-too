"""Vectorized two-sample KS score matching ONCOfind GetES under scipy's asymptotic p-values.

GetES compares ks_2samp p-values, not the statistics. For |G| > 10000, scipy
method='auto' uses the asymptotic one-sided formula, which is strictly
decreasing in the statistic and identical for both alternatives. Therefore
p_greater < p_less iff D+ > D-, and a tie returns the negative less-statistic.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit(cache=True)
def _ks_one(right_row, tie_row, members, n_genes):
    m = members.shape[0]
    r = np.empty(m, dtype=np.int32)
    t = np.empty(m, dtype=np.int32)
    for i in range(m):
        g = members[i]
        r[i] = right_row[g]
        t[i] = tie_row[g]
    order = np.argsort(r)
    max_diff = 0.0
    min_diff = 0.0
    i = 0
    while i < m:
        j = i + 1
        ri = r[order[i]]
        while j < m and r[order[j]] == ri:
            j += 1
        ti = t[order[i]]
        c1_at = float(ri)
        c1_before = float(ri - ti)
        diff_before = c1_before / n_genes - float(i) / m
        diff_at = c1_at / n_genes - float(j) / m
        if diff_before > max_diff:
            max_diff = diff_before
        if diff_at > max_diff:
            max_diff = diff_at
        if diff_before < min_diff:
            min_diff = diff_before
        if diff_at < min_diff:
            min_diff = diff_at
        i = j
    dplus = max_diff
    dminus = -min_diff
    if dminus < 0.0:
        dminus = 0.0
    elif dminus > 1.0:
        dminus = 1.0
    if dplus > dminus:
        return dplus
    return -dminus


@njit(parallel=True, cache=True)
def ks_matrix(right, tie, offsets, indices):
    """right, tie: (n_samples, n_genes) int32. offsets/indices pack gene-set members."""
    n_samples = right.shape[0]
    n_genes = right.shape[1]
    n_sets = offsets.shape[0] - 1
    out = np.empty((n_samples, n_sets), dtype=np.float32)
    for s in prange(n_samples):
        for g in range(n_sets):
            a = offsets[g]
            b = offsets[g + 1]
            out[s, g] = _ks_one(right[s], tie[s], indices[a:b], n_genes)
    return out


@njit(parallel=True, cache=True)
def right_and_tie(values):
    """values: (n_samples, n_genes) float32. right = count of genes with value <= this gene."""
    n_samples, n_genes = values.shape
    right = np.empty((n_samples, n_genes), dtype=np.int32)
    tie = np.empty((n_samples, n_genes), dtype=np.int32)
    for s in prange(n_samples):
        order = np.argsort(values[s])
        i = 0
        while i < n_genes:
            j = i + 1
            vi = values[s, order[i]]
            while j < n_genes and values[s, order[j]] == vi:
                j += 1
            width = j - i
            for k in range(i, j):
                gene = order[k]
                right[s, gene] = j
                tie[s, gene] = width
            i = j
    return right, tie


def pack_sets(member_lists):
    offsets = np.zeros(len(member_lists) + 1, dtype=np.int64)
    chunks = []
    for i, members in enumerate(member_lists):
        arr = np.asarray(members, dtype=np.int32)
        chunks.append(arr)
        offsets[i + 1] = offsets[i] + arr.size
    indices = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int32)
    return offsets, indices


def original_getes(value: np.ndarray, gene_idx: np.ndarray) -> float:
    """ONCOfind GetES, using scipy.stats.ks_2samp defaults (method='auto')."""
    from scipy.stats import ks_2samp

    subset = value[gene_idx]
    ks_g = ks_2samp(value, subset, alternative="greater")
    ks_l = ks_2samp(value, subset, alternative="less")
    if ks_g.pvalue < ks_l.pvalue:
        return float(ks_g.statistic)
    return float(ks_l.statistic) * -1.0
