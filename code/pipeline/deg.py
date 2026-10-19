"""
Differential expression per Leiden cluster (each cluster vs. rest).

- Uses sc.tl.rank_genes_groups on adata.layers['log1p_norm'].
- Saves one combined CSV with all clusters and one CSV per cluster.
- Attaches gene symbols (if a symbol column is present in adata.var).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import anndata as ad
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def compute_deg(
    adata: ad.AnnData,
    groupby: str,
    method: str = "wilcoxon",
    layer: str = "log1p_norm",
    key_added: str = "rank_genes_groups",
) -> ad.AnnData:
    """Run rank_genes_groups on the given clustering; each group vs. rest."""
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' not in adata.obs")
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not in adata.layers")

    n_groups = adata.obs[groupby].nunique()
    if n_groups < 2:
        raise ValueError(
            f"Cannot run DE with only {n_groups} group(s) in '{groupby}'."
        )

    log.info(
        f"Computing DEGs on '{groupby}' ({n_groups} groups) "
        f"method={method}, layer={layer}"
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        reference="rest",
        method=method,
        layer=layer,
        use_raw=False,
        pts=True,
        key_added=key_added,
    )
    log.info(f"DEGs stored in adata.uns['{key_added}']")
    return adata


def deg_results_to_dataframe(
    adata: ad.AnnData,
    key: str = "rank_genes_groups",
    symbol_col: Optional[str] = "gene_symbol",
) -> pd.DataFrame:
    """
    Flatten adata.uns[key] into a long-format DataFrame with one row per
    (cluster, gene). Columns: cluster, rank, gene, symbol (optional), score,
    logfc, pval, pval_adj, pct_in, pct_out.
    """
    if key not in adata.uns:
        raise KeyError(f"adata.uns['{key}'] not found — run compute_deg first.")

    result = adata.uns[key]
    groups: List[str] = list(result["names"].dtype.names)

    # pts / pts_rest are pandas DataFrames in modern scanpy:
    #   index = gene id, columns = cluster label. Older versions may omit them.
    pts_df = result.get("pts", None)
    pts_rest_df = result.get("pts_rest", None)

    frames = []
    for group in groups:
        genes = list(result["names"][group])
        df = pd.DataFrame({
            "cluster":  group,
            "rank":     range(1, len(genes) + 1),
            "gene":     genes,
            "score":    list(result["scores"][group]),
            "logfc":    list(result["logfoldchanges"][group]),
            "pval":     list(result["pvals"][group]),
            "pval_adj": list(result["pvals_adj"][group]),
        })

        if pts_df is not None and group in pts_df.columns:
            df["pct_in"] = df["gene"].map(pts_df[group])
        if pts_rest_df is not None and group in pts_rest_df.columns:
            df["pct_out"] = df["gene"].map(pts_rest_df[group])

        frames.append(df)

    long_df = pd.concat(frames, ignore_index=True)

    # Attach gene symbols if the column exists on adata.var.
    if symbol_col is not None and symbol_col in adata.var.columns:
        symbol_map = adata.var[symbol_col].astype(str).to_dict()
        long_df.insert(3, "symbol", long_df["gene"].map(symbol_map))

    return long_df


def save_deg_tables(
    long_df: pd.DataFrame,
    out_dir: Path,
    top_n_per_cluster: Optional[int] = None,
) -> None:
    """
    Write:
      out_dir/deg_all_clusters.csv                  (long format, all rows)
      out_dir/per_cluster/deg_cluster_<X>.csv       (one file per cluster)
      out_dir/deg_top<N>_per_cluster.csv            (only if top_n_per_cluster set)
    """
    out_dir = Path(out_dir)
    per_cluster_dir = out_dir / "per_cluster"
    per_cluster_dir.mkdir(parents=True, exist_ok=True)

    combined_path = out_dir / "deg_all_clusters.csv"
    long_df.to_csv(combined_path, index=False)
    log.info(f"Combined DEG table: {combined_path}")

    for cluster, sub in long_df.groupby("cluster", observed=True):
        safe = str(cluster).replace(" ", "_").replace("/", "-")
        path = per_cluster_dir / f"deg_cluster_{safe}.csv"
        sub.to_csv(path, index=False)
    log.info(f"Per-cluster DEG tables: {per_cluster_dir}")

    if top_n_per_cluster is not None:
        top = (
            long_df.sort_values(["cluster", "rank"])
                   .groupby("cluster", observed=True)
                   .head(top_n_per_cluster)
        )
        top_path = out_dir / f"deg_top{top_n_per_cluster}_per_cluster.csv"
        top.to_csv(top_path, index=False)
        log.info(f"Top-{top_n_per_cluster} DEG summary: {top_path}")