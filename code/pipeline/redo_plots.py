"""
Recompute plots from the checkpoint h5ad written at the end of the DEG step.

Use this after fixing a plot bug (e.g. gene-symbol mapping) so you don't have
to re-run the 16-hour PCA/UMAP/Leiden/DEG prefix.

Usage:
    python -m pipeline.redo_plots
"""
from __future__ import annotations

import logging
import sys

import matplotlib

matplotlib.use("Agg")

import scanpy as sc

from pipeline import config as cfg
from pipeline import plotting, metadata, dim_reduction
from pipeline.plotting import try_plot


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger("redo_plots")


def main() -> None:
    checkpoint = cfg.SAVE_PATH / "full_adata_after_deg.h5ad"
    if not checkpoint.exists():
        log.error(f"Checkpoint not found: {checkpoint}")
        sys.exit(1)

    log.info(f"Loading checkpoint: {checkpoint}")
    adata = sc.read_h5ad(checkpoint)
    log.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Optionally attach cell type labels from a CSV (barcode -> cell_type).
    if cfg.CELL_TYPE_CSV is not None:
        try:
            metadata.attach_cell_type_csv(
                adata,
                csv_path=cfg.CELL_TYPE_CSV,
                barcode_col=cfg.CELL_TYPE_BARCODE_COL,
                label_col=cfg.CELL_TYPE_LABEL_COL,
                obs_col="cell_type",
            )
        except Exception as exc:  # noqa: BLE001
            log.warning(f"cell_type merge FAILED: {type(exc).__name__}: {exc}")
    else:
        log.info("CELL_TYPE_CSV is None — skipping cell-type merge.")

    # Rebuild gene-symbol mapping using the current config.
    sym_to_ens = plotting.ensure_gene_symbols(
        adata,
        var_symbol_col_candidates=cfg.VAR_SYMBOL_COL_CANDIDATES,
        mapping_tsv=cfg.ENSEMBL_TO_SYMBOL_TSV,
    )

    # Ensure every resolution in the sweep is present in adata.obs.
    # Any resolution already computed (from the pipeline run) is skipped.
    # Newly-added resolutions in config get computed here without touching
    # the expensive PCA/UMAP steps.
    dim_reduction.run_leiden_sweep(
        adata, resolutions=cfg.LEIDEN_RESOLUTIONS, prefix=cfg.LEIDEN_KEY_PREFIX,
    )

    # UMAPs (per Leiden resolution + covariates)
    for res in cfg.LEIDEN_RESOLUTIONS:
        key = dim_reduction.leiden_key(res, prefix=cfg.LEIDEN_KEY_PREFIX)
        try_plot(f"umap_{key}", plotting.plot_umap_leiden, adata, cfg.FIGDIR, key)
    for cov in cfg.UMAP_COVARIATES:
        try_plot(f"umap_{cov}", plotting.plot_umap_covariate, adata, cfg.FIGDIR, cov)

    # Marker gene plots
    try_plot(
        "markers_umap",
        plotting.plot_marker_groups_umap,
        adata, cfg.MARKER_GROUPS, sym_to_ens, cfg.FIGDIR,
    )
    try_plot(
        "markers_dotplot",
        plotting.plot_marker_dotplot_by_leiden,
        adata, cfg.MARKER_GROUPS, sym_to_ens, cfg.LEIDEN_KEY, cfg.FIGDIR,
    )

    # DEG plots (results are already in adata.uns from the previous run)
    try_plot(
        "deg_ranked_genes",
        plotting.plot_deg_ranked_genes,
        adata, cfg.FIGDIR, cfg.DEG_KEY, cfg.DEG_TOP_N_FOR_SUMMARY,
    )
    try_plot(
        "deg_dotplot",
        plotting.plot_deg_dotplot,
        adata, cfg.FIGDIR, cfg.DEG_KEY, cfg.DEG_TOP_N_FOR_DOTPLOT,
    )
    try_plot(
        "deg_heatmap",
        plotting.plot_deg_heatmap,
        adata, cfg.FIGDIR, cfg.DEG_KEY, cfg.DEG_TOP_N_FOR_DOTPLOT,
    )

    # Composition plots: one bar per batch, sections coloured by covariate.
    # Both proportion (normalized) and count versions are saved for each covariate.
    batch_col = adata.uns.get("batch_col") or "Batch_meta"
    log.info(f"Stacked-bar plots against batch column: {batch_col}")
    for cov in cfg.REPLOT_COVARIATES:
        try_plot(
            f"stacked_prop_{cov}",
            plotting.plot_stacked_bar_by_batch,
            adata, batch_col, cov, cfg.FIGDIR, True,
        )
        try_plot(
            f"stacked_count_{cov}",
            plotting.plot_stacked_bar_by_batch,
            adata, batch_col, cov, cfg.FIGDIR, False,
        )

    log.info("Plot regeneration complete.")


if __name__ == "__main__":
    main()