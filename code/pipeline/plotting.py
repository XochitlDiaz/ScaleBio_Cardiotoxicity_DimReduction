"""
All plotting functions live here.

Every function is small and side-effecting (it draws to the current matplotlib
figure and saves it). The orchestrator calls each one via `try_plot`, which
catches any exception, logs it, closes the figure, and moves on — so a broken
plot never crashes the pipeline.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Callable, Dict, Iterable, Mapping, Optional

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

log = logging.getLogger(__name__)


# =============================================================================
# Wrappers
# =============================================================================
def try_plot(name: str, fn: Callable, *args, **kwargs) -> None:
    """Run a plotting function; log-and-continue on any exception."""
    log.info(f"[plot] {name}")
    try:
        fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 — deliberate broad catch
        log.warning(f"[plot] {name} FAILED: {type(exc).__name__}: {exc}")
        plt.close("all")


def _save(figdir: Path, name: str, dpi: int = 200) -> None:
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    path = figdir / f"{name}.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    log.info(f"[plot]   -> {path}")


# =============================================================================
# Gene symbol helpers
# =============================================================================
def ensure_gene_symbols(
    adata: ad.AnnData,
    var_symbol_col_candidates: Iterable[str],
    mapping_tsv: Optional[Path] = None,
    out_col: str = "gene_symbol",
) -> Dict[str, str]:
    """
    Ensure adata.var[out_col] contains a gene symbol for every Ensembl ID
    (falling back to the Ensembl ID itself if no symbol is known). Returns
    a symbol -> Ensembl dict for convenience.

    Sources tried in order:
      1. an existing symbol column in adata.var (from var_symbol_col_candidates)
      2. mapping_tsv (Ensembl ID in first col, symbol in second)
      3. identity fallback
    """
    # 1. existing var column
    src_col = None
    for cand in var_symbol_col_candidates:
        if cand in adata.var.columns:
            src_col = cand
            log.info(f"Using existing var column '{src_col}' as gene symbols.")
            break

    if src_col is not None:
        adata.var[out_col] = adata.var[src_col].astype(str)
    elif mapping_tsv is not None and Path(mapping_tsv).exists():
        log.info(f"Loading Ensembl -> symbol mapping from {mapping_tsv}")
        # usecols=[0, 1] handles both 2-column and 3-column features files
        # (Cell Ranger features.tsv[.gz] has 3 cols: ensembl_id, symbol, feature_type).
        # Pandas auto-detects gzip from the .gz extension.
        df = pd.read_csv(
            mapping_tsv,
            sep="\t",
            header=None,
            usecols=[0, 1],
            names=["ensembl", "symbol"],
            dtype=str,
        )
        # Deduplicate on ensembl_id in case the mapping came from concatenated
        # per-sample features files (last symbol wins).
        df = df.drop_duplicates(subset="ensembl", keep="last")
        m = dict(zip(df["ensembl"], df["symbol"]))
        adata.var[out_col] = pd.Series(adata.var_names, index=adata.var_names).map(m)
        adata.var[out_col] = adata.var[out_col].fillna(adata.var_names.to_series())
        n_mapped = int((adata.var[out_col] != adata.var_names).sum())
        log.info(f"  mapped {n_mapped}/{adata.n_vars} genes to symbols")
    else:
        log.warning(
            "No gene symbol source found — falling back to Ensembl IDs as symbols."
        )
        adata.var[out_col] = adata.var_names.astype(str)

    # Build reverse mapping (symbol -> ensembl). If duplicates exist, the last wins.
    sym_to_ens = dict(zip(adata.var[out_col].astype(str), adata.var_names.astype(str)))
    return sym_to_ens


# =============================================================================
# QC plots
# =============================================================================
def plot_qc_scatter_disease(adata: ad.AnnData, figdir: Path) -> None:
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="ct_cardiovascular_disease", show=False,
    )
    _save(figdir, "qc_scatter_counts_vs_genes_by_disease")


def plot_qc_scatter_mt(adata: ad.AnnData, figdir: Path) -> None:
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="pct_counts_mt", show=False,
    )
    _save(figdir, "qc_scatter_counts_vs_genes_by_mt")


def plot_total_counts_hist(adata: ad.AnnData, figdir: Path) -> None:
    g = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    g.set(xlabel="total_counts", ylabel="cells")
    g.savefig(Path(figdir) / "qc_displot_total_counts.png",
              dpi=200, bbox_inches="tight")
    plt.close("all")


def plot_qc_violins(adata: ad.AnnData, figdir: Path) -> None:
    """One combined figure with the three main QC violins."""
    sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        multi_panel=True, show=False,
    )
    _save(figdir, "qc_violins_combined")


# =============================================================================
# HVG plot
# =============================================================================
def plot_hvg_scatter(adata: ad.AnnData, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=adata.var, x="means", y="dispersions",
                    hue="highly_variable" if "highly_variable" in adata.var else None,
                    s=5, ax=ax, edgecolor=None)
    ax.set_xscale("log")
    ax.set_title("HVG selection")
    _save(figdir, "hvg_means_vs_dispersions")


# =============================================================================
# PCA plots
# =============================================================================
def plot_pca_by(adata: ad.AnnData, figdir: Path, color: str) -> None:
    sc.pl.pca_scatter(adata, color=color, show=False)
    _save(figdir, f"pca_scatter_{color}")


def plot_pca_variance(adata: ad.AnnData, figdir: Path, n_pcs: int = 50) -> None:
    with plt.rc_context({"figure.figsize": (8, 4)}):
        sc.pl.pca_variance_ratio(adata, log=False, n_pcs=n_pcs, show=False)
        _save(figdir, "pca_variance_ratio_linear")
    with plt.rc_context({"figure.figsize": (8, 4)}):
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, show=False)
        _save(figdir, "pca_variance_ratio_log")


# =============================================================================
# UMAP plots
# =============================================================================
def plot_umap_covariate(adata: ad.AnnData, figdir: Path, color: str) -> None:
    """UMAP coloured by one covariate. Skips gracefully if column missing."""
    if color not in adata.obs.columns:
        log.info(f"  skipping UMAP by '{color}' — column not in adata.obs")
        return
    sc.pl.umap(
        adata, color=color, ncols=1, wspace=0.4,
        show=False, layer="log1p_norm", vmax="p99",
    )
    safe_name = color.replace("/", "-").replace(" ", "_")
    _save(figdir, f"umap_{safe_name}")


def plot_umap_leiden(adata: ad.AnnData, figdir: Path, leiden_key: str) -> None:
    if leiden_key not in adata.obs.columns:
        raise KeyError(f"Leiden key '{leiden_key}' not in adata.obs")
    sc.pl.umap(adata, color=leiden_key, legend_loc="on data",
               legend_fontsize=8, show=False)
    _save(figdir, f"umap_{leiden_key}")


# =============================================================================
# Marker-gene plots
# =============================================================================
def _resolve_group(
    genes: Iterable[str],
    sym_to_ens: Mapping[str, str],
    var_names: pd.Index,
) -> tuple[list[str], list[str]]:
    """For a list of gene symbols, return (ensembl_ids_present, symbols_present)."""
    ens, syms = [], []
    for sym in genes:
        eid = sym_to_ens.get(sym)
        if eid is not None and eid in var_names:
            ens.append(eid)
            syms.append(sym)
    return ens, syms


def plot_marker_groups_umap(
    adata: ad.AnnData,
    marker_groups: Mapping[str, Iterable[str]],
    sym_to_ens: Mapping[str, str],
    figdir: Path,
) -> None:
    """One UMAP-panel figure per marker group; symbol titles, Ensembl lookup."""
    for group_name, genes in marker_groups.items():
        ens, syms = _resolve_group(genes, sym_to_ens, adata.var_names)
        if not ens:
            log.info(f"  marker group '{group_name}': no genes present in adata")
            continue

        sc.pl.umap(
            adata, color=ens, title=syms, ncols=3, wspace=0.4,
            show=False, layer="log1p_norm",
        )
        fig = plt.gcf()
        try:
            fig.set_layout_engine("none")  # matplotlib >= 3.6
        except AttributeError:
            fig.set_tight_layout(False)
        fig.suptitle(group_name, fontsize=16, fontweight="bold", y=0.995)
        fig.subplots_adjust(top=0.92)

        safe = group_name.replace(" ", "_").replace("/", "-")
        _save(figdir, f"umap_markers_{safe}")


def plot_marker_dotplot_by_leiden(
    adata: ad.AnnData,
    marker_groups: Mapping[str, Iterable[str]],
    sym_to_ens: Mapping[str, str],
    leiden_key: str,
    figdir: Path,
    symbol_col: str = "gene_symbol",
) -> None:
    """
    Grouped dotplot of marker groups vs. Leiden clusters.

    var_names are passed as {group_label: [ensembl_ids]} so scanpy draws the
    group-bracket annotations across the top of the plot. If a gene-symbol
    column is present in adata.var, symbols are used as axis labels.
    """
    if leiden_key not in adata.obs.columns:
        raise KeyError(f"Leiden key '{leiden_key}' not in adata.obs")

    var_names_grouped: Dict[str, list[str]] = {}
    for group_name, genes in marker_groups.items():
        ens, _ = _resolve_group(genes, sym_to_ens, adata.var_names)
        if ens:
            var_names_grouped[group_name] = ens

    if not var_names_grouped:
        log.warning("Marker dotplot: no genes matched — skipping.")
        return

    kwargs = dict(
        adata=adata,
        var_names=var_names_grouped,
        groupby=leiden_key,
        layer="log1p_norm",
        standard_scale="var",
        dendrogram=True,
        show=False,
    )
    if symbol_col in adata.var.columns:
        kwargs["gene_symbols"] = symbol_col

    sc.pl.dotplot(**kwargs)
    _save(figdir, f"dotplot_markers_by_{leiden_key}")


# =============================================================================
# Cell cycle plots
# =============================================================================
def plot_cell_cycle_umap(adata: ad.AnnData, figdir: Path) -> None:
    """UMAP coloured by cell-cycle phase and the S / G2M continuous scores."""
    required = ("phase", "S_score", "G2M_score")
    missing = [k for k in required if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"Cell cycle columns missing in adata.obs: {missing}")

    sc.pl.umap(
        adata,
        color=list(required),
        title=["Cell Cycle Phase", "S Phase Score", "G2M Phase Score"],
        color_map="viridis",
        frameon=False,
        wspace=0.3,
        show=False,
    )
    _save(figdir, "umap_cell_cycle")


# =============================================================================
# DEG plots (top markers per Leiden cluster)
# =============================================================================
def plot_deg_ranked_genes(
    adata: ad.AnnData,
    figdir: Path,
    key: str = "rank_genes_groups",
    n_genes: int = 25,
) -> None:
    """One panel per cluster showing the top N genes ranked by score."""
    sc.pl.rank_genes_groups(
        adata, n_genes=n_genes, sharey=False, key=key, show=False,
    )
    _save(figdir, f"deg_ranked_genes_top{n_genes}")


def plot_deg_dotplot(
    adata: ad.AnnData,
    figdir: Path,
    key: str = "rank_genes_groups",
    n_genes: int = 5,
    symbol_col: str = "gene_symbol",
) -> None:
    """Dotplot of the top N DEGs per cluster across all clusters."""
    kwargs = dict(
        adata=adata,
        n_genes=n_genes,
        key=key,
        standard_scale="var",
        show=False,
    )
    if symbol_col in adata.var.columns:
        kwargs["gene_symbols"] = symbol_col
    sc.pl.rank_genes_groups_dotplot(**kwargs)
    _save(figdir, f"deg_dotplot_top{n_genes}_per_cluster")


def plot_deg_heatmap(
    adata: ad.AnnData,
    figdir: Path,
    key: str = "rank_genes_groups",
    n_genes: int = 5,
    symbol_col: str = "gene_symbol",
) -> None:
    """Heatmap of top N DEGs per cluster."""
    kwargs = dict(
        adata=adata,
        n_genes=n_genes,
        key=key,
        layer="log1p_norm",
        standard_scale="var",
        show=False,
    )
    if symbol_col in adata.var.columns:
        kwargs["gene_symbols"] = symbol_col
    sc.pl.rank_genes_groups_heatmap(**kwargs)
    _save(figdir, f"deg_heatmap_top{n_genes}_per_cluster")