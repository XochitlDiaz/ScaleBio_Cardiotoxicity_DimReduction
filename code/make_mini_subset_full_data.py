import scanpy as sc
import scanpy.external as sce
import anndata as ad
from pathlib import Path
import pandas as pd
import os
import gc
import matplotlib.pyplot as plt
import numpy as np

sc.settings.n_jobs = 16

# 1.1. Define the path to your text file with an absolute path to .h5ad files (one per line)
txt_file_path = "/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/code/ScaleBio_h5ad_file_paths.txt" # Update this path!

# 1.2. Define output file 
save_path="/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/results/test_v2"

with open(txt_file_path, 'r') as file:
    # .strip() removes hidden newline characters (\n) or spaces at the end of each line
    h5ad_files = [Path(line.strip()) for line in file if line.strip()]


# Read each file and store them in a dictionary
adata_dict = {}
for file_path in h5ad_files:
    # Extract a sample name from the filename (e.g., "Sample-1")
    # Adjust the replacement string if your naming convention is different
    sample_name = file_path.name.replace('_anndata.h5ad', '')
    
    print(f"Loading {sample_name}...")
    adata_dict[sample_name] = sc.read_h5ad(file_path)

# Secure Ensembl IDs & Perform Safe Outer Concatenation
print("Building global master feature list to preserve Ensembl IDs...")
all_features = pd.Index([])
for sample_adata in adata_dict.values():
    all_features = all_features.union(sample_adata.var_names)
print(f"Total unique features across all batches: {len(all_features)}")

print("Reindexing individual samples to uniform master layout...")
for sample_name, sample_adata in adata_dict.items():
    # Slicing maps missing features to 0 counts in .X and NaN in .var metadata
    adata_dict[sample_name] = sample_adata[:, all_features].copy()


# Concatenate all AnnData objects into one
adata = ad.concat(
    adata_dict, 
    join="inner",         
    label="sample_id"    
)

adata.obs_names = adata.obs_names.str.rsplit('-Sample').str[0]

del adata_dict
gc.collect()

BioT_metadata = pd.read_csv("/g/data/wy83/cardiotoxicity/analysis/anndata/data/Toxcicity/ScaleBio/Cardiotoxicity_per_cell_metadata.tsv",sep='\t')
BioT_metadata["_index"] = BioT_metadata["_index"].str.replace(r'.*anndata_', '', regex=True)
BioT_metadata= BioT_metadata.set_index("_index")
adata.obs = adata.obs.join(BioT_metadata, how='left', lsuffix='_adata', rsuffix='_meta')

BioT_metadata2 = pd.read_csv('/g/data/wy83/cardiotoxicity/analysis/anndata/BioTuring/TBD_metadata.tsv', usecols=[0, 1, 22], sep='\t')
BioT_metadata2["_index"] = BioT_metadata2["_index"].str.replace(r'.*anndata_', '', regex=True)
BioT_metadata2 = BioT_metadata2.set_index("_index")
adata.obs = adata.obs.join(BioT_metadata2, how='left', lsuffix='_adata', rsuffix='_meta')

if 'Batch' in adata.obs.columns:
    # 1. Fill missing unmapped values with a string placeholder
    adata.obs['Batch'] = adata.obs['Batch'].fillna('Unknown')
    # 2. Cast the entire column explicitly to string text
    adata.obs['Batch'] = adata.obs['Batch'].astype(str)

target_cols = adata.obs.select_dtypes(include=['object', 'category']).columns

for col in target_cols:
    # 1. Cast to a plain string first (this eliminates the strict categorical restrictions)
    adata.obs[col] = adata.obs[col].astype(str)
    
    # 2. Safely replace any text representations of missing data with "Unknown"
    # (Since it's now a string, 'nan' will be a literal string value from the original NaNs)
    adata.obs[col] = adata.obs[col].replace(['nan', 'NaN', 'None', 'NONE', ''], 'Unknown')
    
    # 3. Explicit check to catch any leftover underlying null positions
    adata.obs[col] = adata.obs[col].fillna('Unknown')

### Cell Filtering ###

print(f"Number of cells: {len(adata)}")

# Doublet filtering
adata = adata[~(adata.obs['Donor_ID'] == 'doublet')]

print(f"Filtered cells after doublet filtering: {len(adata)}")

# Not expected donors filtering
# 1. Define the villages you want to target for exclusion
bad_villages = ["Village_21", "Village_27"]

# 2. Check if a row has a bad village in either column
has_bad = (
    adata.obs['village_id'].isin(bad_villages) | 
    adata.obs['village_id_2'].isin(bad_villages)
)

# 3. Helper function to check if a column contains a valid "saving" village
def has_valid_other(col_series):
    # A valid "other" is not a bad village, not NaN, and not an empty/null string
    return (
        ~col_series.isin(bad_villages) & 
        col_series.notna() & 
        ~col_series.isin(["", "nan", "NaN", "NONE"])
    )

# 4. Check if the row has a saving village in either column
has_other = has_valid_other(adata.obs['village_id']) | has_valid_other(adata.obs['village_id_2'])

# 5. The drop condition: they have a bad village AND lack a saving village.
# We want to KEEP the opposite of the drop condition.
keep_mask = ~(has_bad & ~has_other)

# 6. Apply the filter to create your new anndata object
adata = adata[keep_mask]

# Print a summary of what happened
print(f"Filtered cells after unexpected donors: {len(adata)}")


# Remove cells with low gene diversity
sc.pp.filter_cells(adata, min_genes=2000)
print(f"Remaining cells after filtering for low gene diversity: {adata.n_obs}")

# Remove genes expressed in few cells
print(f"Total features before filtering (genes): {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=20)
print(f"Remaining genes after filtering for low cell expression: {adata.n_vars}")

# REmove high mt percentage cells
adata = adata[adata.obs["pct_counts_mt"] <= 6, :].copy()
print(f"Remaining cells after filtering for cell with < 5% MT counts: {adata.n_obs}")


### Define Mitochondiral genes ###
MT_genes_to_ensemble_ids= {
    'MT-ND1'  : "ENSG00000198888",
    'MT-ND2'  : "ENSG00000198763",
    'MT-CO1'  : "ENSG00000198804",
    'MT-CO2'  : "ENSG00000198712",
    'MT-ATP8' : "ENSG00000228253",
    'MT-ATP6' : "ENSG00000198899",
    'MT-CO3'  : "ENSG00000198938",
    'MT-ND3'  : "ENSG00000198840",
    'MT-ND4L' : "ENSG00000212907",
    'MT-ND4'  : "ENSG00000198886",
    'MT-ND5'  : "ENSG00000198786",
    'MT-ND6'  : "ENSG00000198695",
    'MT-CYB'  : "ENSG00000198727"
}
adata.var["mt"] = adata.var_names.isin(MT_genes_to_ensemble_ids.values())



### Calculate QC metrics ###
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=["mt"], 
    inplace=True, 
    percent_top=[20], 
    log1p=True
)


### Normalize and log-transform the data ###
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)


# Find highly variable genes, perform PCA, and scale the data
sc.pp.highly_variable_genes(adata, layer="log1p_norm")
hvg_df = adata.var[adata.var['highly_variable']]

print("Number of highly variable genes: " + str(len(hvg_df)))

#Compute PCs
sc.pp.pca(adata, svd_solver="arpack",use_highly_variable=True, layer="log1p_norm")

hvg_cols = [c for c in adata.var.columns if 'highly_variable' in c or c in ['means', 'dispersions', 'dispersions_norm']]

print(f"All full data CSV files successfully saved to: {csv_export_path}")

### MAKE THE SUBSET ###
total_cells = adata.n_obs
print(f"Original dataset size: {total_cells} cells.")

target_size = total_cells // 16
print(f"Target subset size: {target_size} cells.")

# 3. Perform a truly random selection
# Setting a seed means if you re-run this, you get the exact same random pool
rng = np.random.default_rng(seed=42)

# Pick random row integers without replacement
random_row_indices = rng.choice(total_cells, size=target_size, replace=False)

# 4. Slice the AnnData object and explicitly call .copy() to allocate new memory
adata_subset = adata[random_row_indices].copy()
print(f"Subset successfully created: {adata_subset.n_obs} cells remaining.")
adata_subset.obs.index.name = None

# 5. Save the new downsampled object to disk
print("Saving subsetted AnnData object to disk...")
adata_subset.write_h5ad(f'{save_path}/mini_full_data_subset.h5ad')
print(f"Successfully saved to: {save_path}")
