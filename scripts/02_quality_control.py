'''
Perform quality control on the microglial AD subset from 01_data_pipeline.py. 

This script performs the following tasks:
- Assesses cell quality metrics (number of genes, number of counts, mitochondrial transcript percentage)
- Removes genes that are expressed in fewer than 10 cells
- Saves the filtered dataset to a new file
'''

from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc

# Load data (run 01_data_pipeline.py first to generate input file)
adata = sc.read_h5ad("microglia_subset.h5ad") 

# Create output directory if it doesn't exist
Path("results/qc").mkdir(parents=True, exist_ok=True)

# Explore gene and cell metrics
summary = adata.obs.groupby("AD_status")[["n_genes", "n_counts"]].describe()
print("\n Cell quality summary by AD status:")
print(summary)

sc.pl.violin(adata, ["n_genes", "n_counts"], groupby="AD_status", multi_panel=True)
plt.savefig("results/qc/genes_counts.png")
plt.close()
# Distributions of detected genes per cell and total counts per cell are similar between AD and control groups. Downstream analysis/comparison will not be driven by technical differences in cell quality or sequencing depth.

# Mitochondrial transcript percentage check
adata.var["mt"] = adata.var["gene_name"].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

#adata.obs["pct_counts_mt"].describe()
mt_summary =adata.obs.groupby("AD_status")[["pct_counts_mt"]].describe()
print("\n Mitochondrial transcript percentage summary:")
print(mt_summary)
#sc.pl.violin(adata, ["pct_counts_mt"], groupby="AD_status")

plt.figure(figsize=(6, 4))
adata.obs["pct_counts_mt"].hist(bins=100)

plt.xlabel("Mitochondrial transcript percentage")
plt.ylabel("# of cells")
plt.tight_layout()
plt.savefig("results/qc/mitochondrial_histogram.png")
plt.close()
# Distributions of mitochondrial transcript percentage between AD and control groups are comparable, with a mean of 0.198 and 0.189 respectively, suggesting that RNA is stable in both groups.

# Mitochondrial percentage max > 1.0, therefore check head
# adata.obs["pct_counts_mt"].sort_values(ascending=False).head() # Probably rounding errors (% > 1.0) in three cells

# Gene filtering
# Mitcohondrial percentages are low across cells (median = 0.13, mean = 0.19), so filtering by mitochondria percentage is not necessary.

# Filter cells with low gene counts and genes with low detection

#adata.obs[["n_genes", "n_counts", "pct_counts_mt"]].describe()

# Lowest gene counts
#adata.obs.sort_values("n_genes").head() # Lowest n_genes = 1001, seems like dataset is already filtered.

#adata.var["n_cells"].describe() # Reflects original atlas rather than subset
adata.var["n_cells_by_counts"].describe()

print("Genes before filtering:", adata.n_vars)
print("Genes expressed in 0 cells:", (adata.var["n_cells_by_counts"] == 0).sum())
print("Genes expressed in < 10 cells:", (adata.var["n_cells_by_counts"] < 10).sum()) # genes detected in ~ 0.2% of cells

sc.pp.filter_genes(adata, min_cells=10)

print("Genes after filtering:", adata.n_vars) # Filtered out 9835 genes

# Clean and save filtered dataset
adata.obs.isnull().sum()

adata.obs_names.is_unique # returns True

print("\n Final dataset dimensions:")
print(adata) # Check dimensions (6000 x 24341)

# Save as new file
output_path = "microglia_AD_qc.h5ad"
adata.write(output_path)


print("Quality Control Summary")
print("-------------------------------")
print(f"Cells (not filtered): {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"AD cells: {(adata.obs['AD_status'] == 'Yes').sum()}")
print(f"Control cells: {(adata.obs['AD_status'] == 'No').sum()}")
print(f"Unique donors: {adata.obs['donor_id'].nunique()}")
print(f"Median # of genes: {adata.obs['n_genes'].median():.0f}")
print(f"Median # of counts: {adata.obs['n_counts'].median():.0f}")
print(f"Median mitochondrial %: {adata.obs['pct_counts_mt'].median():.2f}%")
