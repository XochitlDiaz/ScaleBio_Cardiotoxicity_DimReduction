"""
QC-metric computation.
"""
from __future__ import annotations

import logging
from typing import Mapping

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def annotate_mt(adata: ad.AnnData, mt_ensembl_ids: Mapping[str, str]) -> ad.AnnData:
    """Flag MT genes in adata.var['mt'] using their Ensembl IDs."""
    ids = set(mt_ensembl_ids.values())
    adata.var["mt"] = adata.var_names.isin(ids)
    n_found = int(adata.var["mt"].sum())
    log.info(f"MT genes present in adata: {n_found}/{len(ids)}")
    return adata


def calculate_qc(adata: ad.AnnData) -> ad.AnnData:
    """Compute standard QC metrics with %counts in top 20 genes and log1p versions."""
    if "mt" not in adata.var.columns:
        raise RuntimeError("annotate_mt() must be called before calculate_qc().")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True,
        percent_top=[20],
        log1p=True,
    )
    log.info("QC metrics computed.")
    return adata
