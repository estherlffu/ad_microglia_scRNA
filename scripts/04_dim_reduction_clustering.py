import scanpy as sc

# Load preprocessed dataset from "03_preprocessing.py"
adata = sc.read_h5ad("data/microglia_AD_preprocessed.h5ad")

# PCA
sc.tl.pca(adata)
#sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False)

# Plot PCA to check
sc.pl.pca(adata, color=["AD_status", "AD_status", "pct_counts_mt", "pct_counts_mt"], dimensions=[(0, 1), (2, 3), (0, 1), (2, 5)], ncols=2, size=2,)

# Build kNN for UMAP and Leiden clustering
sc.pp.neighbors(adata, n_pcs=30)

# UMAP: 2d kNN mapping
sc.tl.umap(adata)
sc.pl.umap(adata, color="AD_status", size=3)

# Leiden clustering: network community detection
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2) # tune resolution, standards default 0.5 to 1
sc.pl.umap(adata, color="leiden", title=f"Leiden clustering (resolution = 1.0)")

#sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)
#sc.pl.umap(adata, color="leiden", title=f"Leiden clustering (resolution = 0.5)")

#sc.tl.leiden(adata, resolution=1.5, flavor="igraph", n_iterations=2)
#sc.pl.umap(adata, color="leiden", title=f"Leiden clustering (resolution = 1.5)")

#sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
#sc.pl.umap(adata, color=["leiden"])

# Interpret Leiden clustering
adata.obs["leiden"].value_counts().sort_index()

# Save clustered analyses
adata.write("data/microglia_AD_clustered.h5ad")