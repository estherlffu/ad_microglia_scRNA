# ad_microglia_scRNA

DATA: RADC_Cohort

Download dataset here: "https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad"

The RADC_Cohort dataset is a "population-scale cross-disorder atlas of the human [dorsolateral] prefrontal cortex at single-cell resolution". The full dataset consists of 693,682 cells.

The data pipeline filters the full dataset by disease (Alzheimer disease and normal) and cell type (microglia). Due to RAM limitation, a stratified random sample of 3000 cells per group is used to produce the subset. Our final dataset consists of 6000 total cells with 34,176 genes, with 147 unique donor IDs.