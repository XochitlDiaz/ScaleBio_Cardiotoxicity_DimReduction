"""
CSV / h5ad export routines.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def export_pca_and_qc_csvs(adata: ad.AnnData, csv_dir: Path) -> None:
    """Dump PCA loadings & coordinates, QC metrics, and HVG table to CSV."""
    csv_dir = Path(csv_dir)
    csv_dir.mkdir(parents=True, exist_ok=True)

    # PCA loadings (var x PCs) — key is usually 'PCs' but scanpy has used 'gsl' historically
    loadings_key = "PCs" if "PCs" in adata.varm else "gsl"
    if loadings_key in adata.varm:
        loadings = pd.DataFrame(
            adata.varm[loadings_key],
            index=adata.var_names,
            columns=[f"PC_{i+1}" for i in range(adata.varm[loadings_key].shape[1])],
        )
        loadings.to_csv(csv_dir / "pca_loadings.csv")

    # PCA coordinates (obs x PCs)
    if "X_pca" in adata.obsm:
        coords = pd.DataFrame(
            adata.obsm["X_pca"],
            index=adata.obs_names,
            columns=[f"PC_{i+1}" for i in range(adata.obsm["X_pca"].shape[1])],
        )
        coords.to_csv(csv_dir / "pca_coordinates.csv")

    # HVG table
    hvg_cols = [
        c for c in adata.var.columns
        if "highly_variable" in c or c in ("means", "dispersions", "dispersions_norm")
    ]
    if hvg_cols:
        adata.var[hvg_cols].to_csv(csv_dir / "highly_variable_genes.csv")

    # QC metrics per cell / per gene
    qc_var_cols = [
        c for c in adata.var.columns
        if c.startswith("n_cells_by_counts")
        or c.startswith("mean_counts")
        or c.startswith("pct_dropout_by_counts")
    ]
    qc_obs_cols = [
        c for c in adata.obs.columns
        if c.startswith("n_genes_by_counts")
        or c.startswith("total_counts")
        or c.startswith("pct_counts")
    ]
    if qc_obs_cols:
        adata.obs[qc_obs_cols].to_csv(csv_dir / "cell_qc_metrics.csv")
    if qc_var_cols:
        adata.var[qc_var_cols].to_csv(csv_dir / "gene_qc_metrics.csv")

    log.info(f"CSV exports written to: {csv_dir}")


def save_per_batch(
    adata: ad.AnnData,
    output_dir: Path,
    batch_col: Optional[str],
) -> None:
    """Write one .h5ad per unique batch value."""
    if batch_col is None:
        log.warning("save_per_batch: no batch column resolved — skipping.")
        return
    if batch_col not in adata.obs.columns:
        log.warning(f"save_per_batch: '{batch_col}' not in adata.obs — skipping.")
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Per-batch export dir: {output_dir}")

    skip_tokens = {"Unknown", "nan", "NaN"}
    unique_batches = sorted(adata.obs[batch_col].astype(str).unique())

    for batch in unique_batches:
        if batch in skip_tokens:
            continue
        log.info(f"Processing batch: {batch}")
        sub = adata[adata.obs[batch_col].astype(str) == batch].copy()
        safe = str(batch).replace(" ", "_").replace("/", "-")
        out = output_dir / f"batch_subset_{safe}.h5ad"
        sub.write_h5ad(out)
        log.info(f"  -> {sub.n_obs} cells written to {out}")

    log.info("All batch subsets exported.")