"""
Normalization, log-transformation, and HVG selection.
"""
from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def normalize_and_log(adata: ad.AnnData) -> ad.AnnData:
    """Median-normalize (target_sum=None) then store log1p into layers['log1p_norm']."""
    scaled = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scaled["X"], copy=True)
    log.info("Stored normalized log1p matrix in adata.layers['log1p_norm'].")
    return adata


def find_hvg(adata: ad.AnnData) -> ad.AnnData:
    """Highly variable genes based on the log1p_norm layer."""
    sc.pp.highly_variable_genes(adata, layer="log1p_norm")
    n_hvg = int(adata.var["highly_variable"].sum())
    log.info(f"Highly variable genes: {n_hvg}")
    return adata
