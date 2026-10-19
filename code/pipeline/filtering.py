"""
Cell and gene filtering steps.

QC-metric-dependent filters (like %mt) must run AFTER qc.calculate_qc().
The pipeline orchestrator handles that ordering; see run_pipeline.py.
"""
from __future__ import annotations

import logging
from typing import Iterable

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def filter_doublets(adata: ad.AnnData) -> ad.AnnData:
    """Drop cells whose Donor_ID is literally 'doublet'."""
    before = adata.n_obs
    if "Donor_ID" not in adata.obs.columns:
        log.warning("Donor_ID column missing — skipping doublet filter.")
        return adata
    adata = adata[~(adata.obs["Donor_ID"] == "doublet")].copy()
    log.info(f"Doublet filter: {before} -> {adata.n_obs} cells")
    return adata


def filter_bad_villages(adata: ad.AnnData, bad_villages: Iterable[str]) -> ad.AnnData:
    """
    Drop cells assigned to unexpected villages.

    Keep a cell if it has a bad-village hit BUT also has a valid non-bad village
    in either village_id or village_id_2 (i.e. the "saving" village). Otherwise
    drop it.
    """
    before = adata.n_obs
    if "village_id" not in adata.obs.columns or "village_id_2" not in adata.obs.columns:
        log.warning("village_id / village_id_2 missing — skipping village filter.")
        return adata

    bad = list(bad_villages)
    has_bad = (
        adata.obs["village_id"].isin(bad) | adata.obs["village_id_2"].isin(bad)
    )

    null_tokens = ["", "nan", "NaN", "NONE", "Unknown"]

    def _valid(col):
        return (
            ~col.isin(bad)
            & col.notna()
            & ~col.isin(null_tokens)
        )

    has_other = _valid(adata.obs["village_id"]) | _valid(adata.obs["village_id_2"])
    keep = ~(has_bad & ~has_other)
    adata = adata[keep].copy()
    log.info(f"Bad-village filter: {before} -> {adata.n_obs} cells")
    return adata


def filter_min_genes(adata: ad.AnnData, min_genes: int) -> ad.AnnData:
    before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    log.info(f"min_genes>={min_genes}: {before} -> {adata.n_obs} cells")
    return adata


def filter_min_cells(adata: ad.AnnData, min_cells: int) -> ad.AnnData:
    before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    log.info(f"min_cells>={min_cells}: {before} -> {adata.n_vars} genes")
    return adata

def filter_min_counts(adata: ad.AnnData, min_counts: int) -> ad.AnnData:
    before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    log.info(f"min_counts>={min_counts}: {before} -> {adata.n_obs} cells")
    return adata


def filter_mt(adata: ad.AnnData, max_pct: float) -> ad.AnnData:
    """Drop cells with pct_counts_mt > max_pct. Requires calculate_qc first."""
    if "pct_counts_mt" not in adata.obs.columns:
        raise RuntimeError(
            "filter_mt called before calculate_qc_metrics; pct_counts_mt is not set."
        )
    before = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct, :].copy()
    log.info(f"pct_counts_mt<={max_pct}: {before} -> {adata.n_obs} cells")
    return adata
