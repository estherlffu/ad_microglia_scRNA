# ad_microglia_scRNA

# Microglial Transcriptomics Analysis in Alzheimer's Disease

## Overview
This project examines the question: "How does Alzheimer's disease alter gene expression and mechanistic pathways in microglia?" by analyzing single-cell RNA-sequence data of human dorsolateral prefrontal cortex. It will explore transcriptional differences in microglia immune cells between Alzheimer's disease and control (non-Alzheimer's disease) samples.

## Dataset
The project uses the [RADC Cohort dataset](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f) hosted by CZ CELLxGENE Discover. The complete RADC dataset contains 693,682 cells, from which this project analyzes a balanced subset of 6,000 microglial cells (3,000 Alzheimer's disease and 3,000 control cells). :contentReference[oaicite:1]{index=1}

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

## References
1. Lee D, Koutrouli M, Masse NY, et al.; PsychAD Consortium. Single-cell atlas of transcriptomic vulnerability across multiple neurodegenerative and neuropsychiatric diseases. medRxiv. 2024. https://doi.org/10.1101/2024.10.31.24316513
2.  Mathys H, Davila-Velderrain J, Peng Z, et al. Single-cell transcriptomic analysis of Alzheimer's disease. Nature. 2019;570(7761):332–337. https://doi.org/10.1038/s41586-019-1195-2
3. Wolf FA, Angerer P, Theis FJ. *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology. 2018;19:15. https://doi.org/10.1186/s13059-017-1382-0
