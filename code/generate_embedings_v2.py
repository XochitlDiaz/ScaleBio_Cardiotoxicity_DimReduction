import scanpy as sc
import scanpy.external as sce
import anndata as ad
from pathlib import Path
import pandas as pd
import os
import gc
import matplotlib.pyplot as plt

sc.settings.n_jobs = 16

# 1.1. Define the path to your text file with an absolute path to .h5ad files (one per line)
txt_file_path = "/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/code/ScaleBio_h5ad_file_paths.txt" # Update this path!

# 1.2. Define output file 
save_path="/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/results/test_v2"

# 1.3. Define number of PCs to use for the UMAP
pc_counts = [6]

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


# Concatenate all AnnData objects into one
adata = ad.concat(
    adata_dict, 
    join="outer",         
    label="sample_id",    
    index_unique="-"      
)

del adata_dict
gc.collect()

print(f"Combined AnnData shape: {adata.shape}")




### Load metadata ###
BioT_metadata = pd.read_csv("/g/data/wy83/cardiotoxicity/analysis/anndata/data/Toxcicity/ScaleBio/Cardiotoxicity_per_cell_metadata.tsv",sep='\t')
BioT_metadata["_index"] = BioT_metadata["_index"].str.replace(r'.*anndata_', '', regex=True)
BioT_metadata= BioT_metadata.set_index("_index")
adata.obs = adata.obs.join(BioT_metadata, how='left', lsuffix='_adata', rsuffix='_meta')

BioT_metadata2 = pd.read_csv('/g/data/wy83/cardiotoxicity/analysis/anndata/BioTuring/TBD_metadata.tsv', usecols=[0, 1, 22], sep='\t')
BioT_metadata2["_index"] = BioT_metadata2["_index"].str.replace(r'.*anndata_', '', regex=True)
BioT_metadata2 = BioT_metadata2.set_index("_index")
adata.obs = adata.obs.join(BioT_metadata2, how='left', lsuffix='_adata', rsuffix='_meta')


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
sc.pp.filter_cells(adata, min_genes=500)
print(f"Remaining cells after filtering for low gene diversity: {adata.n_obs}")

# Remove genes expressed in few cells
print(f"Total features before filtering (genes): {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=20)
print(f"Remaining genes after filtering for low cell expression: {adata.n_vars}")



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

sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.gcf().savefig(
    os.path.join(save_path, "pca_variance_ratio.png"),
    bbox_inches="tight",
    dpi=300,
)
plt.close()


# Save PC coordinates to CSV
pca_df = pd.DataFrame(
    adata.obsm['X_pca'], 
    index=adata.obs_names,
    columns=[f"PC_{i+1}" for i in range(adata.obsm['X_pca'].shape[1])]
)


pca_csv_file = f"{save_path}/pca_coords.csv"
pca_df.to_csv(pca_csv_file)
print(f"Saved: {pca_csv_file}")

# Save gene loadings for the PCs
is_hvg = adata.var['highly_variable']
hvg_loadings = adata.varm['PCs'][is_hvg, :]

loadings_df = pd.DataFrame(
    hvg_loadings, 
    index=adata.var_names[is_hvg], 
    columns=[f"PC_{i+1}" for i in range(adata.varm['PCs'].shape[1])]
)

loadings_df.to_csv(os.path.join(save_path, "pca_gene_loadings.csv"))
print(f"Saved: {os.path.join(save_path, 'pca_gene_loadings.csv')}")


# Run harmony integration
print("Running Harmony integration...")
sce.pp.harmony_integrate(
    adata, 
    #key=["Batch", "village_id", "disease_status"]
    key="Batch"
)

# Save the Harmony adjusted coordinates to disk
harmony_df = pd.DataFrame(
    adata.obsm['X_pca_harmony'], 
    index=adata.obs_names,

    columns=[f"Harmony_{i+1}" for i in range(adata.obsm['X_pca_harmony'].shape[1])]
)

harmony_csv_file = f"{save_path}/harmony_coords.csv"
harmony_df.to_csv(harmony_csv_file)
print(f"Saved: {harmony_csv_file}")


for row_idx, n_pcs in enumerate(pc_counts):

    # Compute neighbors and UMAP per PC setting
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Save UMAP coordinates to CSV
    umap_df = pd.DataFrame(
        adata.obsm['X_umap'],
        index=adata.obs_names,
        columns=['UMAP_1', 'UMAP_2']
    )
    csv_file = f"{save_path}/umap_coords_{n_pcs}_PCs_harmony_Batch.csv"
    umap_df.to_csv(csv_file)
    print(f"Saved: {csv_file}")

adata.write("g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/results/clean_full_toxicity_adata.h5ad")