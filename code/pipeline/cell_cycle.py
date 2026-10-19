"""
Cell cycle scoring using the Tirosh et al. 2016 S / G2M gene lists.

`sc.tl.score_genes_cell_cycle` reads from `adata.X`, but in this pipeline `X`
holds raw counts and the log-normalized matrix lives in `adata.layers["log1p_norm"]`.
We temporarily swap `X` for the scoring call so the scores are computed on
log-normalized values, then restore `X`.
"""
from __future__ import annotations

import logging
from typing import Iterable, Mapping

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def _match_ensembl(
    symbols: Iterable[str],
    sym_to_ens: Mapping[str, str],
    var_names,
) -> list[str]:
    return [
        sym_to_ens[s] for s in symbols
        if s in sym_to_ens and sym_to_ens[s] in var_names
    ]


def score_cell_cycle(
    adata: ad.AnnData,
    s_genes_symbols: Iterable[str],
    g2m_genes_symbols: Iterable[str],
    sym_to_ens: Mapping[str, str],
    layer: str | None = "log1p_norm",
) -> ad.AnnData:
    """
    Compute S_score, G2M_score, and phase (G1/S/G2M) in `adata.obs`.

    Symbols are mapped to Ensembl IDs via `sym_to_ens` (from
    `plotting.ensure_gene_symbols`) so this works with Ensembl var_names.
    """
    s_ens = _match_ensembl(s_genes_symbols, sym_to_ens, adata.var_names)
    g2m_ens = _match_ensembl(g2m_genes_symbols, sym_to_ens, adata.var_names)
    log.info(
        f"Matched {len(s_ens)}/{len(list(s_genes_symbols))} S-phase and "
        f"{len(g2m_ens)}/{len(list(g2m_genes_symbols))} G2M-phase genes."
    )

    if not s_ens or not g2m_ens:
        raise RuntimeError(
            "Cell cycle scoring aborted — no matching genes for one of the sets."
        )

    if layer is not None and layer in adata.layers:
        original_X = adata.X
        adata.X = adata.layers[layer]
        try:
            sc.tl.score_genes_cell_cycle(
                adata, s_genes=s_ens, g2m_genes=g2m_ens,
            )
        finally:
            adata.X = original_X
    else:
        sc.tl.score_genes_cell_cycle(
            adata, s_genes=s_ens, g2m_genes=g2m_ens,
        )

    log.info(
        "Cell cycle phase distribution:\n"
        f"{adata.obs['phase'].value_counts().to_string()}"
    )
    return adata
