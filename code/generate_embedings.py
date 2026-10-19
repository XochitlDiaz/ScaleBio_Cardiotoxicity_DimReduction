import scanpy as sc
import scanpy.external as sce
import anndata as ad
from pathlib import Path
import pandas as pd
import os
import gc


sc.settings.n_jobs = 16

# 1.1. Define the path to your text file with an absolute path to .h5ad files (one per line)
txt_file_path = "/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/code/ScaleBio_h5ad_file_paths_test.txt" # Update this path!

# 1.2. Define output file 
save_path="/g/data/wy83/cardiotoxicity/analysis/anndata/full_data_analysis/dim_reduction/results/test"

# 1.3. Define number of PCs to use for the UMAP
pc_counts = [11]

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


# Load features files
features_file = "/g/data/wy83/cardiotoxicity/data/Cardiotoxicity_transcriptomics/GIMR-Neavin-SO10000428-QuantumScale/processed/samples/Sample-1.merged.filtered.matrix/features.tsv.gz"
features = pd.read_csv(features_file, sep="\t", header=None, names=["ensembl_id", "gene_symbol", "feature_type"])
ensembl_to_gene = dict(zip(features["ensembl_id"], features["gene_symbol"]))
adata.var["ensembl_id"] = adata.var_names
adata.var_names = adata.var_names.to_series().map(ensembl_to_gene).fillna(adata.var_names.to_series())

# Find duplicated gene symbols and append the Ensembl ID to make them unique
# 1. Map the Ensembl IDs to Gene Symbols
adata.var_names = adata.var_names.to_series().map(ensembl_to_gene).fillna(adata.var_names.to_series())

# 2. Make Genes unique
adata.var_names_make_unique()

# 3. Find what the raw symbols were supposed to be
raw_symbols = adata.var["ensembl_id"].map(ensembl_to_gene).fillna(adata.var["ensembl_id"])
is_duplicated = raw_symbols.duplicated(keep=False)

# 4. Build the tracker using .values so pandas doesn't panic over the new unique indices
duplicate_tracker = pd.DataFrame({
    "Ensembl_ID": adata.var["ensembl_id"].values,
    "Raw_Symbol": raw_symbols.values,
    "Final_Modified_Symbol": adata.var_names.values  # <--- This will now contain the -1 modifications
})

# 5. Filter and print
collided_genes = duplicate_tracker[is_duplicated.values].sort_values("Raw_Symbol")

BioT_metadata = pd.read_csv("/g/data/wy83/cardiotoxicity/analysis/anndata/BioTuring/TBD_metadata.tsv",sep='\t')
BioT_metadata["_index"] = BioT_metadata["_index"].str.replace(r'.*anndata_', '', regex=True)
BioT_metadata= BioT_metadata.set_index("_index")
adata.obs = adata.obs.join(BioT_metadata, how='left', lsuffix='_adata', rsuffix='_meta')


metadata = pd.read_csv("/g/data/wy83/cardiotoxicity/data/Cardiotoxicity_transcriptomics/BioTuring_upload/BBrowserX_formatted_metadata.tsv",sep='\t')
metadata["_index"] = metadata["_index"].str.replace(r'.*anndata_', '', regex=True)
metadata= metadata.set_index("_index")
adata.obs = adata.obs.join(metadata, how='left', lsuffix='_o', rsuffix='_m')

donor_assignment = pd.read_csv("/g/data/wy83/cardiotoxicity/data/Cardiotoxicity_transcriptomics/BioTuring_upload/merged_donor_assignments/master_donor_dictionary.tsv",sep='\t',header=None)
first_col = donor_assignment.columns[1]
donor_assignment.loc[donor_assignment[first_col].str.contains('+', regex=False, na=False), first_col] = 'doublet'
donor_assignment = donor_assignment.rename(columns={0: '_index', 1: 'Donor'})
donor_assignment = donor_assignment.set_index('_index')
adata.obs = adata.obs.join(donor_assignment, how='left')

adata = adata[adata.obs['Donor'] != 'doublet'].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

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
loadings_df = pd.DataFrame(
    adata.varm['PCs'], 
    index=adata.var_names[adata.var['highly_variable']], 
    columns=[f"PC_{i+1}" for i in range(adata.varm['PCs'].shape[1])]
)
loadings_df.to_csv(f"{save_path}/pca_gene_loadings.csv")


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

