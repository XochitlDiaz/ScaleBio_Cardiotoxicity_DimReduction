"""
Master pipeline runner.

Usage:
    python -m pipeline.run_pipeline

Every plot is called via `try_plot`, so a single failing plot logs a warning
and the rest of the pipeline continues.
"""
from __future__ import annotations

import logging
import sys

import matplotlib

matplotlib.use("Agg")  # headless — must come before pyplot is imported anywhere else

import scanpy as sc

from pipeline import config as cfg
from pipeline import (
    io_utils,
    metadata,
    filtering,
    qc,
    normalization,
    dim_reduction,
    plotting,
    export,
    deg,
    subset,
    cell_cycle,
)
from pipeline.plotting import try_plot


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger("run_pipeline")


def main() -> None:
    sc.settings.n_jobs = cfg.N_JOBS

    # -------------------------------------------------------------------------
    # 1. Load & concatenate
    # -------------------------------------------------------------------------
    log.info("=== Loading & concatenating samples ===")
    adata_dict = io_utils.load_h5ad_dict(cfg.H5AD_LIST_PATH)
    adata = io_utils.concat_samples(adata_dict)
    log.info(f"After concat: {adata.n_obs} cells x {adata.n_vars} genes")

    # -------------------------------------------------------------------------
    # 2. Attach metadata & clean obs strings
    # -------------------------------------------------------------------------
    log.info("=== Attaching metadata ===")
    adata = metadata.attach_metadata(
        adata,
        tsv_path_1=cfg.METADATA_TSV_1,
        tsv_path_2=cfg.METADATA_TSV_2,
        tsv_2_cols=cfg.METADATA_TSV_2_COLS,
    )
    adata = metadata.clean_obs_strings(adata)
    batch_col = metadata.resolve_batch_column(adata, cfg.BATCH_COL_CANDIDATES)

    # -------------------------------------------------------------------------
    # 3. Categorical / metadata-based filters (no QC metrics needed yet)
    # -------------------------------------------------------------------------
    log.info("=== Categorical filters ===")
    adata = filtering.filter_doublets(adata)
    adata = filtering.filter_bad_villages(adata, cfg.BAD_VILLAGES)

    # -------------------------------------------------------------------------
    # 4. Basic count-based filters
    # -------------------------------------------------------------------------
    log.info("=== Count-based filters ===")
    adata = filtering.filter_min_genes(adata, cfg.MIN_GENES_PER_CELL)
    adata = filtering.filter_min_cells(adata, cfg.MIN_CELLS_PER_GENE)
    adata = filtering.filter_min_counts(adata, cfg.MIN_COUNTS_PER_CELL)

    # -------------------------------------------------------------------------
    # 5. Annotate MT genes -> compute QC metrics -> filter on %MT
    #    (this ordering is important — the previous script had a bug here)
    # -------------------------------------------------------------------------
    log.info("=== QC metrics ===")
    adata = qc.annotate_mt(adata, cfg.MT_GENES_TO_ENSEMBL)
    adata = qc.calculate_qc(adata)

    # QC plots BEFORE MT filtering so the pre-filter distribution is preserved.
    try_plot("qc_scatter_disease",   plotting.plot_qc_scatter_disease, adata, cfg.FIGDIR)
    try_plot("qc_scatter_mt",        plotting.plot_qc_scatter_mt,      adata, cfg.FIGDIR)
    try_plot("qc_displot_counts",    plotting.plot_total_counts_hist,  adata, cfg.FIGDIR)
    try_plot("qc_violins",           plotting.plot_qc_violins,         adata, cfg.FIGDIR)

    adata = filtering.filter_mt(adata, cfg.MAX_PCT_MT)

    # -------------------------------------------------------------------------
    # 6. Normalization, log1p, HVG
    # -------------------------------------------------------------------------
    log.info("=== Normalization / HVG ===")
    adata = normalization.normalize_and_log(adata)
    adata = normalization.find_hvg(adata)

    try_plot("hvg_scatter", plotting.plot_hvg_scatter, adata, cfg.FIGDIR)

    # Build the symbol -> Ensembl mapping now so downstream steps (cell cycle
    # scoring, marker plots) can share it. Also populates adata.var["gene_symbol"].
    sym_to_ens = plotting.ensure_gene_symbols(
        adata,
        var_symbol_col_candidates=cfg.VAR_SYMBOL_COL_CANDIDATES,
        mapping_tsv=cfg.ENSEMBL_TO_SYMBOL_TSV,
    )

    # -------------------------------------------------------------------------
    # 7. Cell cycle scoring (Tirosh gene lists) + UMAP by phase / scores
    # -------------------------------------------------------------------------
    log.info("=== Cell cycle scoring ===")
    try:
        cell_cycle.score_cell_cycle(
            adata,
            s_genes_symbols=cfg.S_GENES,
            g2m_genes_symbols=cfg.G2M_GENES,
            sym_to_ens=sym_to_ens,
            layer=cfg.CELL_CYCLE_LAYER,
        )
    except Exception as exc:  # noqa: BLE001 — non-fatal, keep pipeline running
        log.warning(f"Cell cycle scoring FAILED: {type(exc).__name__}: {exc}")

    # -------------------------------------------------------------------------
    # 8. PCA + variance plots
    # -------------------------------------------------------------------------
    log.info("=== PCA ===")
    adata = dim_reduction.run_pca(adata)

    try_plot("pca_variance", plotting.plot_pca_variance, adata, cfg.FIGDIR, 50)
    try_plot("pca_by_village", plotting.plot_pca_by, adata, cfg.FIGDIR, "village_id")

    # -------------------------------------------------------------------------
    # 9. Neighbors + UMAP + Leiden
    # -------------------------------------------------------------------------
    log.info("=== Neighbors / UMAP / Leiden ===")
    adata = dim_reduction.run_neighbors(
        adata, n_pcs=cfg.N_PCS_FOR_NEIGHBORS, n_neighbors=cfg.N_NEIGHBORS,
    )
    adata = dim_reduction.run_umap(adata)
    adata = dim_reduction.run_leiden_sweep(
        adata, resolutions=cfg.LEIDEN_RESOLUTIONS, prefix=cfg.LEIDEN_KEY_PREFIX,
    )

    # Save a random 1/N subset immediately once UMAP+Leiden are done so we
    # have a lightweight file for quick exploration even if a later step fails.
    try:
        subset.save_mini_adata(
            adata,
            save_path=cfg.SAVE_PATH,
            denominator=cfg.SUBSET_DENOMINATOR,
            seed=cfg.RANDOM_SEED,
        )
    except Exception as exc:  # noqa: BLE001
        log.warning(f"Mini adata save FAILED: {type(exc).__name__}: {exc}")

    # UMAP for each Leiden resolution in the sweep (side-by-side comparison)
    for res in cfg.LEIDEN_RESOLUTIONS:
        key = dim_reduction.leiden_key(res, prefix=cfg.LEIDEN_KEY_PREFIX)
        try_plot(
            f"umap_{key}",
            plotting.plot_umap_leiden, adata, cfg.FIGDIR, key,
        )

    # UMAP by each covariate
    for cov in cfg.UMAP_COVARIATES:
        try_plot(f"umap_{cov}", plotting.plot_umap_covariate, adata, cfg.FIGDIR, cov)

    # Cell-cycle UMAP (phase + S / G2M scores were computed pre-UMAP in step 7)
    try_plot("umap_cell_cycle", plotting.plot_cell_cycle_umap, adata, cfg.FIGDIR)

    # -------------------------------------------------------------------------
    # 10. Marker gene plots (UMAP panels + dotplot vs. Leiden)
    # -------------------------------------------------------------------------
    log.info("=== Marker gene plots ===")

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

    # -------------------------------------------------------------------------
    # 11. Differential expression per Leiden cluster
    # -------------------------------------------------------------------------
    log.info("=== DEG per Leiden cluster ===")
    try:
        deg.compute_deg(
            adata,
            groupby=cfg.LEIDEN_KEY,
            method=cfg.DEG_METHOD,
            layer="log1p_norm",
            key_added=cfg.DEG_KEY,
        )
        deg_df = deg.deg_results_to_dataframe(
            adata, key=cfg.DEG_KEY, symbol_col="gene_symbol",
        )
        deg.save_deg_tables(
            deg_df, out_dir=cfg.DEG_DIR,
            top_n_per_cluster=cfg.DEG_TOP_N_FOR_SUMMARY,
        )
    except Exception as exc:  # noqa: BLE001 — non-fatal, keep pipeline running
        log.warning(f"DEG computation FAILED: {type(exc).__name__}: {exc}")

    # DEG plots (each wrapped independently — one failing plot doesn't kill the rest)
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

    # Checkpoint: save the fully-processed adata now (PCA, neighbors, UMAP,
    # Leiden, DEG all populated). If a downstream step needs to be re-run,
    # you can start from this file instead of redoing the 16-hour prefix.
    try:
        checkpoint = cfg.SAVE_PATH / "full_adata_after_deg.h5ad"
        adata.write_h5ad(checkpoint)
        log.info(f"Checkpoint written: {checkpoint}")
    except Exception as exc:  # noqa: BLE001
        log.warning(f"Checkpoint save FAILED: {type(exc).__name__}: {exc}")

    # -------------------------------------------------------------------------
    # 12. CSV exports & disk saves
    # -------------------------------------------------------------------------
    log.info("=== CSV exports ===")
    export.export_pca_and_qc_csvs(adata, cfg.CSV_EXPORT_PATH)

    log.info("=== Per-batch export ===")
    export.save_per_batch(adata, cfg.PER_BATCH_DIR, batch_col=batch_col)

    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()