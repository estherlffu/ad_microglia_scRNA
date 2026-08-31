import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load normalized data from "03_preprocessing.py"
adata = sc.read_h5ad("data/microglia_AD_normalized.h5ad") # Need non-scaled normalized data for log2 computation

# Load clustered data from "04_dim_reduction_clustering.py"
clustered = sc.read_h5ad("data/microglia_AD_clustered.h5ad")

# Marker gene identification
# adata.X.min(), adata.X.max() # Scaled data

adata.obs["leiden"] = clustered.obs["leiden"]
adata.obsm["X_umap"] = clustered.obsm["X_umap"]

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# Visualize marker genes
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, gene_symbols="gene_name")

# Check that genes are variable across clusters
adata.var_names = adata.var["gene_name"]
adata.var_names_make_unique()

sc.pl.violin(
    adata,
    ["APOE", "TREM2", "SPP1", "GPNMB"],
    groupby="leiden", rotation=90
)

adata.obs_names.equals(clustered.obs_names)

markers = sc.get.rank_genes_groups_df(adata, group=None)

markers["gene_symbol"] = markers["names"].map(adata.var["gene_name"])

markers.to_csv("/content/drive/MyDrive/Projects/Single-cell RNA-seq AD Project/results/marker_genes_all_clusters.csv", index=False)

# Check known microglia markers
microglia_markers = ["P2RY12", "TMEM119", "CX3CR1", "CSF1R", "AIF1","TYROBP", "TREM2", "APOE", "LPL", "SPP1", "GPNMB", "SYT1", "NRXN1", "CNTNAP2", "CX3CR1"]
sc.pl.dotplot(adata, microglia_markers, groupby="leiden", gene_symbols="gene_name")
sc.pl.umap(adata, color=["P2RY12", "TREM2", "APOE", "SPP1"], gene_symbols="gene_name")

# Choose clusters based on annotations from exploring top 30 genes per cluster
keep_clusters = ["1", "2", "3", "6", "7", "12", "14", "15", "17"]
adata_clean = adata[adata.obs["leiden"].isin(keep_clusters)].copy()

# Re-cluster
sc.pp.highly_variable_genes(adata_clean, n_top_genes=2000)
sc.tl.pca(adata_clean)
sc.pp.neighbors(adata_clean, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_clean, resolution=0.5)
sc.tl.umap(adata_clean)

# Visualize re-clustering
sc.pl.umap(adata_clean, color="leiden", legend_loc="on data", title="Leiden re-cluster")
microglia_markers = ["P2RY12", "TMEM119", "CX3CR1", "CSF1R", "AIF1","TYROBP", "TREM2", "APOE", "LPL", "SPP1", "GPNMB", "SYT1", "NRXN1", "CNTNAP2", "CX3CR1"]
sc.pl.dotplot(adata_clean, microglia_markers, groupby="leiden", standard_scale="var")
#sc.pl.dotplot(adata_clean, microglia_markers, groupby="leiden")

sc.tl.rank_genes_groups(adata_clean, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata_clean, n_genes=25, sharey=False)
markers = sc.get.rank_genes_groups_df(adata_clean, group=None)
markers.head()

adata_clean.obs["leiden"].value_counts()

# markers[markers["group"] == "8"].head(20)[["names", "scores", "logfoldchanges"]]

# Check AD status by cluster
pd.crosstab(adata_clean.obs["leiden"], adata_clean.obs["AD_status"], normalize="index")

cluster_AD = pd.crosstab(adata_clean.obs["leiden"], adata_clean.obs["AD_status"], normalize="index")

cluster_AD.plot(kind="bar", stacked=True)
plt.ylabel("Fraction of cells")
plt.xlabel("Leiden cluster")
plt.title("AD composition across microglial clusters")
plt.show()

sc.tl.rank_genes_groups(adata_clean, groupby="AD_status", method="wilcoxon")
sc.pl.rank_genes_groups(adata_clean, n_genes=25, sharey=False)

result = adata_clean.uns["rank_genes_groups"]

groups = result["names"].dtype.names

for group in groups:
    df = pd.DataFrame({
        "gene": result["names"][group],
        "score": result["scores"][group],
        "logfoldchange": result["logfoldchanges"][group],
        "pval": result["pvals"][group],
        "pval_adj": result["pvals_adj"][group]
    })

    print("GROUP:", group)

sc.pl.dotplot(adata_clean, microglia_markers, groupby="AD_status", standard_scale="var")

adata_clean.obs["annotated_cell_type"] = adata_clean.obs["leiden"].map({
    "0": "Homeostatic microglia",
    "1": "Intermediate microglia",
    "2": "Activated/stressed microglia",
    "3": "Phagocytic/lipid-associated microglia",
    "4": "DAM/SPP1 inflammatory",
    "5": "Activated/stressed microglia",
    "6": "Interferon microglia",
    "7": "Excluded donor-specific",
    "8": "Excluded donor-specific"
})

adata_clean.obs.to_csv("data/microglia_cluster_annotations.csv")

adata_clean.write("data/microglia_clusters.h5ad")