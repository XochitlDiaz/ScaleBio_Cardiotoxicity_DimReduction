"""
Random cell subsampling.

Creates a downsampled AnnData containing n_cells // denominator cells picked
uniformly at random (without replacement, seeded). Used to produce a lightweight
version of the full dataset for quick exploration in notebooks.
"""
from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np

log = logging.getLogger(__name__)


def make_random_subset(
    adata: ad.AnnData,
    denominator: int,
    seed: int,
) -> ad.AnnData:
    """Return a new AnnData with n_cells // denominator randomly-chosen cells.

    Uses np.random.default_rng(seed) with replace=False, so the subset size is
    exact and the same cells are picked on every re-run with the same seed.
    """
    total = adata.n_obs
    target = total // denominator
    log.info(f"Random subset: {total} -> {target} cells (seed={seed}).")

    rng = np.random.default_rng(seed=seed)
    idx = rng.choice(total, size=target, replace=False)

    subset = adata[idx].copy()
    subset.obs.index.name = None
    return subset


def save_mini_adata(
    adata: ad.AnnData,
    save_path: Path,
    denominator: int,
    seed: int,
    filename: str = "mini_full_adata.h5ad",
) -> Path:
    """Downsample and write to disk. Returns the output path."""
    save_path = Path(save_path)
    save_path.mkdir(parents=True, exist_ok=True)

    subset = make_random_subset(adata, denominator=denominator, seed=seed)

    out = save_path / filename
    subset.write_h5ad(out)
    log.info(f"Mini AnnData saved: {out} ({subset.n_obs} cells)")
    return out