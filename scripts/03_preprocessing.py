import scanpy as sc

# Load filtered dataset from "02_quality_control.py"
adata = sc.read_h5ad("microglia_AD_qc.h5ad")

# Normalize cell counts to median
sc.pp.normalize_total(adata, target_sum=1e4) # Counts per 10k

# sc.pl.highest_expr_genes(adata, n_top=20) # Check whether highly expressed genes are worth excluding

# Logarithmize data
sc.pp.log1p(adata)

# Save normalized data before feature selection (filtering for highly variable/most informative genes)
adata.raw = adata.copy()

# Identify top 2000 highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(adata.var["highly_variable"].value_counts()) # Confirm selection

# Visualize HVGs
sc.pl.highly_variable_genes(adata)

# Subset to HVGs
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"Dataset after HVG selection: {adata.n_obs} cells × {adata.n_vars} genes")

# Scale data in preparation for PCA and clustering
sc.pp.scale(adata, max_value=10) # max_value = 10 to limit effect of extreme outliers ahead of PCA

# Save preprocessed dataset
output_path = "microglia_AD_preprocessed.h5ad"
adata.write(output_path)