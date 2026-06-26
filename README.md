# ad_microglia_scRNA

# Microglial Transcriptomics Analysis in Alzheimer's Disease

## Overview
This project examines the question: "How does Alzheimer's disease alter gene expression in microglia?" by analyzing single-cell RNA-sequence data from human dorsolateral prefrontal cortex. It will explore transcriptional differences in microglia immune cells between Alzheimer's disease and control samples.

## Dataset
The RADC_Cohort dataset by CELLxGENE Census is a "population-scale cross-disorder atlas of the human [dorsolateral] prefrontal cortex at single-cell resolution". The full dataset consists of 693,682 cells.

The data pipeline filters the full dataset by disease (Alzheimer's disease and normal) and cell type (microglia). Due to RAM limitations, a stratified random sample of 3000 cells per group is used to produce the subset. Our final dataset consists of 6000 total cells with 34,176 genes, and 147 unique donor IDs.

Download RADC_Cohort dataset here: "https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad"

## Pipeline
1. Data preprocessing
2. Quality Control
3. Preprocessing
4. Dimensionality Reduction
5. Clustering and Annotation
6. Differential Expression
7. Pathway Analysis

## Citations
SCANPY: large-scale single-cell gene expression data analysis
F. Alexander Wolf, Philipp Angerer, Fabian J. Theis
Genome Biology 2018 Feb 06. doi: 10.1186/s13059-017-1382-0.
