"""
Attach per-cell metadata TSVs to adata.obs and normalize missing values.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

import anndata as ad
import pandas as pd

log = logging.getLogger(__name__)


def _load_metadata(tsv_path: Path, usecols: Optional[Iterable[int]] = None) -> pd.DataFrame:
    """Load a metadata TSV and set '_index' as the row index, stripping any
    '<something>_anndata_' prefix that upstream tools add to the barcodes."""
    read_kwargs = {"sep": "\t"}
    if usecols is not None:
        read_kwargs["usecols"] = list(usecols)
    df = pd.read_csv(tsv_path, **read_kwargs)

    if "_index" not in df.columns:
        raise KeyError(f"{tsv_path} is missing a '_index' column")

    df["_index"] = df["_index"].str.replace(r".*anndata_", "", regex=True)
    df = df.set_index("_index")
    return df


def attach_metadata(
    adata: ad.AnnData,
    tsv_path_1: Path,
    tsv_path_2: Path,
    tsv_2_cols: Iterable[int],
) -> ad.AnnData:
    """Left-join two metadata TSVs onto adata.obs by cell barcode."""
    log.info(f"Reading metadata TSV 1: {tsv_path_1}")
    meta1 = _load_metadata(tsv_path_1)

    overlap = adata.obs_names.isin(meta1.index).sum()
    log.info(f"  overlap with adata: {overlap}/{adata.n_obs} cells")
    if overlap == 0:
        log.warning(
            "Zero index overlap with the first metadata file — check the "
            "barcode format in the TSV vs. adata.obs_names."
        )

    adata.obs = adata.obs.join(meta1, how="left", lsuffix="_adata", rsuffix="_meta")

    log.info(f"Reading metadata TSV 2: {tsv_path_2}")
    meta2 = _load_metadata(tsv_path_2, usecols=tsv_2_cols)
    adata.obs = adata.obs.join(meta2, how="left", lsuffix="_adata", rsuffix="_meta")

    return adata


def clean_obs_strings(adata: ad.AnnData) -> ad.AnnData:
    """Cast object/category columns to plain strings and replace null-ish
    tokens with 'Unknown'. Prevents downstream categorical write errors."""
    if "Batch" in adata.obs.columns:
        adata.obs["Batch"] = adata.obs["Batch"].fillna("Unknown").astype(str)

    target_cols = adata.obs.select_dtypes(include=["object", "category"]).columns
    for col in target_cols:
        adata.obs[col] = adata.obs[col].astype(str)
        adata.obs[col] = adata.obs[col].replace(
            ["nan", "NaN", "None", "NONE", ""], "Unknown"
        )
        adata.obs[col] = adata.obs[col].fillna("Unknown")

    return adata


def resolve_batch_column(adata: ad.AnnData, candidates: Iterable[str]) -> Optional[str]:
    """Return the first candidate name that actually exists in adata.obs.

    Also stores the resolved name in adata.uns['batch_col'] so downstream
    modules can look it up without re-guessing.
    """
    for name in candidates:
        if name in adata.obs.columns:
            adata.uns["batch_col"] = name
            log.info(f"Using batch column: {name}")
            return name
    log.warning(
        f"None of {list(candidates)} found in adata.obs. "
        "Batch-dependent steps will be skipped."
    )
    adata.uns["batch_col"] = None
    return None
