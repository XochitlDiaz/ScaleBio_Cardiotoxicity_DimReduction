"""
PCA, neighbor graph, UMAP, and Leiden clustering.
"""
from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def run_pca(adata: ad.AnnData) -> ad.AnnData:
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable", layer="log1p_norm")
    log.info(f"PCA computed: {adata.obsm['X_pca'].shape[1]} components.")
    return adata


def run_neighbors(adata: ad.AnnData, n_pcs: int, n_neighbors: int) -> ad.AnnData:
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    log.info(f"Neighbor graph built (n_pcs={n_pcs}, n_neighbors={n_neighbors}).")
    return adata


def run_umap(adata: ad.AnnData) -> ad.AnnData:
    sc.tl.umap(adata)
    log.info("UMAP embedding computed.")
    return adata


def leiden_key(resolution: float, prefix: str = "leiden_r") -> str:
    """Canonical column name for a Leiden clustering at a given resolution."""
    return f"{prefix}{resolution}"


def run_leiden(adata: ad.AnnData, resolution: float, key_added: str) -> ad.AnnData:
    """Run Leiden and store labels in adata.obs[key_added] as a categorical.

    Uses flavor='igraph' (python-igraph's built-in Leiden) so we don't need
    the standalone `leidenalg` package. n_iterations=2 is required because
    the igraph flavor defaults to a single iteration.
    """
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        flavor="igraph",
        n_iterations=2,
    )
    adata.obs[key_added] = adata.obs[key_added].astype("category")
    n_clusters = adata.obs[key_added].nunique()
    log.info(
        f"Leiden clustering @ resolution={resolution}: "
        f"{n_clusters} clusters -> adata.obs['{key_added}']"
    )
    return adata


def run_leiden_sweep(
    adata: ad.AnnData,
    resolutions: list[float],
    prefix: str = "leiden_r",
    skip_existing: bool = True,
) -> ad.AnnData:
    """
    Run Leiden at each resolution in the sweep.

    Each clustering is stored under adata.obs[f"{prefix}{res}"]. If
    skip_existing=True, resolutions whose column already exists are not
    recomputed — useful in redo_plots for exploring new resolutions without
    re-clustering the ones already in the checkpoint.
    """
    log.info(f"Leiden sweep across resolutions: {resolutions}")
    for res in resolutions:
        key = leiden_key(res, prefix=prefix)
        if skip_existing and key in adata.obs.columns:
            log.info(f"  {key} already present — skipping.")
            continue
        run_leiden(adata, resolution=res, key_added=key)
    return adata