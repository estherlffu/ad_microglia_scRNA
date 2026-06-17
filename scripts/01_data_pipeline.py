import scanpy as sc
import numpy as np

# Load data (download dataset in README if not already downloaded)
adata = sc.read_h5ad("data.h5ad", backed='r')
obs = adata.obs

# Preprocess data (filtered for microglial cells and AD status)
micro = obs[obs["cell_type"] == "microglial cell"]
micro = micro[micro["AD_status"].isin(["Yes", "No"])]

# Stratify samples by AD
ad_ids = micro.index[micro["AD_status"] == "Yes"]
ctrl_ids = micro.index[micro["AD_status"] == "No"]

n = 3000

ad_sample = np.random.choice(ad_ids, n, replace=False)
ctrl_sample = np.random.choice(ctrl_ids, n, replace=False)

cell_ids = np.concatenate([ad_sample, ctrl_sample])

# Subset original data and write to new file
adata_subset = adata[cell_ids]
adata_subset = sc.AnnData(
    X=adata_subset.X[:],
    obs=adata_subset.obs.copy(),
    var=adata_subset.var.copy()
)

# Save dataset
output_path = "microglia_subset.h5ad"
adata_subset.write(output_path)

# Quick quality checks
print("\nAD status counts:")
print(adata_subset.obs["AD_status"].value_counts())

print("\nCell type counts:")
print(adata_subset.obs["cell_type"].value_counts())

print("\nPreview:")
print(adata_subset.obs.head(2))