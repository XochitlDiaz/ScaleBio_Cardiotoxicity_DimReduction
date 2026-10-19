"""
Central configuration for the scRNA-seq dim-reduction pipeline.

Edit paths / thresholds here; everything else imports from this module.
"""
from pathlib import Path

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
#H5AD_LIST_PATH = Path(
#    "/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/"
#    "dim_reduction/code/ScaleBio_h5ad_file_paths.txt"
#)

H5AD_LIST_PATH = Path(
    "/g/data/wy83/cardiotoxicity/analysis/anndata/data/Toxcicity/"
    "TestData/ScaleBio_h5ad_file_paths_minitest.txt"
)

METADATA_TSV_1 = Path(
    "/g/data/wy83/cardiotoxicity/analysis/anndata/data/Toxcicity/ScaleBio/"
    "Cardiotoxicity_per_cell_metadata.tsv"
)

METADATA_TSV_2 = Path(
    "/g/data/wy83/cardiotoxicity/analysis/anndata/BioTuring/TBD_metadata.tsv"
)
METADATA_TSV_2_COLS = [0, 1, 22]  # columns to read from TSV_2

SAVE_PATH = Path(
    "/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/"
    "dim_reduction/results_v2"
)
FIGDIR = SAVE_PATH / "figures"
CSV_EXPORT_PATH = SAVE_PATH / "full_csv_exports"
PER_BATCH_DIR = SAVE_PATH / "per_sample_adata"
DEG_DIR = SAVE_PATH / "deg"

# Optional Ensembl -> gene-symbol mapping (TSV with two columns: ensembl_id, symbol)
# Leave as None if adata.var already has a symbol column (see VAR_SYMBOL_COL_CANDIDATES)
ENSEMBL_TO_SYMBOL_TSV: Path | None = Path(
    "/g/data/wy83/cardiotoxicity/data/Cardiotoxicity_transcriptomics/"
    "GIMR-Neavin-SO10000428-QuantumScale/processed/samples/"
    "Sample-1.merged.filtered.matrix/features.tsv.gz"
)
VAR_SYMBOL_COL_CANDIDATES = ("feature_name", "gene_symbols", "gene_symbol",
                             "symbol", "gene_name", "Symbol")

# -----------------------------------------------------------------------------
# Compute
# -----------------------------------------------------------------------------
N_JOBS = 6
RANDOM_SEED = 42

# -----------------------------------------------------------------------------
# Filtering thresholds
# -----------------------------------------------------------------------------
BAD_VILLAGES = ["Village_21", "Village_27"]
MIN_GENES_PER_CELL = 600
MIN_COUNTS_PER_CELL = 3500
MIN_CELLS_PER_GENE = 20
MAX_PCT_MT = 6.0

# -----------------------------------------------------------------------------
# Column names
# -----------------------------------------------------------------------------
# Which obs column carries the batch label to split on when exporting per-batch.
# After the metadata join this may end up as "Batch" or "Batch_meta" depending
# on whether the original adata.obs already had "Batch". The metadata module
# resolves this and stores the final name in adata.uns["batch_col"].
BATCH_COL_CANDIDATES = ("Batch_meta", "Batch")

# Leiden clustering key stored in adata.obs after running the leiden step.
LEIDEN_RESOLUTIONS = [0.1, 0.15, 0.18, 0.2]
LEIDEN_PRIMARY_RESOLUTION = 0.15
LEIDEN_KEY_PREFIX = "leiden_r"
LEIDEN_KEY = f"{LEIDEN_KEY_PREFIX}{LEIDEN_PRIMARY_RESOLUTION}"  # "leiden_r0.15"

# -----------------------------------------------------------------------------
# Dim reduction
# -----------------------------------------------------------------------------
N_PCS_FOR_NEIGHBORS = 11
N_NEIGHBORS = 30

# -----------------------------------------------------------------------------
# Subsetting
# -----------------------------------------------------------------------------
SUBSET_DENOMINATOR = 16  # keep n_cells // 16 cells at random

# -----------------------------------------------------------------------------
# DEG
# -----------------------------------------------------------------------------
DEG_METHOD = "wilcoxon"           # 'wilcoxon' | 't-test' | 't-test_overestim_var' | 'logreg'
DEG_KEY = "rank_genes_groups"     # where results live in adata.uns
DEG_TOP_N_FOR_SUMMARY = 25        # top-N table + per-cluster ranked plot
DEG_TOP_N_FOR_DOTPLOT = 5         # top-N shown in the DEG dotplot / heatmap

# -----------------------------------------------------------------------------
# UMAP covariates to colour by (each becomes one saved figure)
# -----------------------------------------------------------------------------
UMAP_COVARIATES = [
    "Batch_meta",              # will fall back gracefully if the col doesn't exist
    "sample_id",
    "ct_cardiovascular_disease",
    "village_id",
    "Donor_ID",
    "pct_counts_mt",
    "n_genes_by_counts",
    "total_counts",
]

# -----------------------------------------------------------------------------
# Mitochondrial genes (Ensembl IDs, since var_names are Ensembl)
# -----------------------------------------------------------------------------
MT_GENES_TO_ENSEMBL = {
    "MT-ND1":  "ENSG00000198888",
    "MT-ND2":  "ENSG00000198763",
    "MT-CO1":  "ENSG00000198804",
    "MT-CO2":  "ENSG00000198712",
    "MT-ATP8": "ENSG00000228253",
    "MT-ATP6": "ENSG00000198899",
    "MT-CO3":  "ENSG00000198938",
    "MT-ND3":  "ENSG00000198840",
    "MT-ND4L": "ENSG00000212907",
    "MT-ND4":  "ENSG00000198886",
    "MT-ND5":  "ENSG00000198786",
    "MT-ND6":  "ENSG00000198695",
    "MT-CYB":  "ENSG00000198727",
}

# -----------------------------------------------------------------------------
# Cell cycle scoring (Tirosh et al. 2016 / Regev lab gene lists; symbols)
# -----------------------------------------------------------------------------
# Layer used for scoring — cell-cycle scoring expects log-normalized values.
CELL_CYCLE_LAYER = "log1p_norm"

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# -----------------------------------------------------------------------------
# Marker gene groups (symbols; converted to Ensembl at plot time)
# -----------------------------------------------------------------------------
MARKER_GROUPS = {
    "Endothelium": [
        "PCDH12", "CDH5", "CD93", "PECAM1", "VIP", "DIPK2B", "MYCT1", "CD34",
        "GIMAP4", "PIK3CG", "RASGRP3", "ROBO4", "TAL1", "BTNL9", "AC114757.1",
    ],
    "Glial": [
        "SOX10", "CDH19", "S100B", "PRSS56", "MPZ", "INSC", "PLP1", "CHRND",
        "FOXD3", "GPR17", "ATP10B", "ERBB3", "ROPN1B", "AC099489.1", "MOXD1",
    ],
    "Mesenchymal": [
        "AC093849.1", "TCF21", "AP003555.2", "AC103702.2", "HAND2", "HAND2-AS1",
        "HOXA13", "RXFP2", "HAND1", "HOXA10", "TARID", "FOXF1", "CCN4",
        "CCDC68", "AC103702.1",
    ],
    "Fibroblasts": [
        "CRYM", "LINC02188", "SLC9A9-AS1", "FOXD1", "AC009154.1", "PRRX1", "MKX",
        "AC116563.1", "LAMA4", "Z99289.2", "FAP", "AC018450.1", "ARHGAP20",
        "CPED1", "ABCA9",
    ],
    "Neurons": [
        "SLC17A6", "DRGX", "NEUROD1", "LINC02487", "STMN2", "UNCX", "AC010247.2",
        "AC091980.2", "AL031658.2", "PHOX2B-AS1", "NEUROD4", "XKR7", "TH",
        "ADCYAP1", "LINC00682",
    ],
    "Neuron Progenitors": [
        "WNT1", "WNT3A", "WNT7A", "RFX4", "AL359091.1", "WNT8B", "ILDR2",
        "WNT10B", "EFCC1", "C6orf118", "RMST", "KCNK10", "PAX7", "OLIG3", "WNT7B",
    ],
    "Epithelium": [
        "GRHL2", "ESRP1", "PRSS22", "TMPRSS4", "ARHGEF38", "TP63", "CDH1",
        "TACSTD2", "ZNF750", "VGLL1", "MACC1", "MUC15", "COL17A1", "CLDN4",
        "MAB21L4",
    ],
    "Cardiac muscle": [
        "MYH6", "TTN", "LRRC10", "MYBPC3", "XIRP1", "SMYD1", "MIR133A1HG",
        "AC020909.2", "TNNT2", "AC022916.2", "LMOD2", "BANCR", "UNC45B",
        "AL136171.2", "PLN",
    ],
}

