"""
Loading and concatenation utilities.
"""
from __future__ import annotations

import gc
import logging
from pathlib import Path
from typing import Dict

import anndata as ad
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def load_h5ad_dict(txt_file_path: Path) -> Dict[str, ad.AnnData]:
    """
    Read a plain-text file of absolute .h5ad paths (one per line) and return a
    {sample_name: AnnData} dict. The sample name is the file basename with
    the '_anndata.h5ad' suffix stripped.
    """
    with open(txt_file_path, "r") as fh:
        h5ad_files = [Path(line.strip()) for line in fh if line.strip()]

    adata_dict: Dict[str, ad.AnnData] = {}
    for path in h5ad_files:
        sample_name = path.name.replace("_anndata.h5ad", "")
        log.info(f"Loading {sample_name}...")
        adata_dict[sample_name] = sc.read_h5ad(path)

    log.info(f"Loaded {len(adata_dict)} samples.")
    return adata_dict


def concat_samples(adata_dict: Dict[str, ad.AnnData]) -> ad.AnnData:
    """
    Reindex each sample onto the union of features, then outer-concatenate.

    Uses join="outer" (not "inner") so the union-reindex above is not undone.
    Missing features become zeros in .X.

    Cell barcodes are assumed to be unique across batches; ad.concat is called
    WITHOUT index_unique so the raw barcodes are preserved (this is what the
    downstream metadata joins expect).
    """
    log.info("Building global master feature list to preserve Ensembl IDs...")
    all_features = pd.Index([])
    for sample_adata in adata_dict.values():
        all_features = all_features.union(sample_adata.var_names)
    log.info(f"Total unique features across all batches: {len(all_features)}")

    log.info("Reindexing individual samples to the union feature layout...")
    for sample_name, sample_adata in adata_dict.items():
        adata_dict[sample_name] = sample_adata[:, all_features].copy()

    log.info("Concatenating samples (join='outer', no barcode suffix)...")
    adata = ad.concat(
        adata_dict,
        join="outer",
        label="sample_id",
    )

    if not adata.obs_names.is_unique:
        n_dup = adata.n_obs - adata.obs_names.nunique()
        raise ValueError(
            f"Barcode uniqueness assumption violated: {n_dup} duplicate cell "
            "barcodes across batches. Re-enable index_unique in concat or fix "
            "upstream barcode generation."
        )

    # Free the per-sample copies as soon as we can — big memory saver.
    adata_dict.clear()
    gc.collect()

    return adata
