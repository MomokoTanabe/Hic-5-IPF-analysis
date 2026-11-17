# Single-cell RNA-seq Analysis of Human Lung Tissue

## Overview

This repository contains R scripts for the analysis of publicly available single-cell RNA sequencing (scRNA-seq) datasets from human lung tissue.

The analyses focus on the expression of **TGFB1I1 (Hic-5)** in lung cell populations, particularly fibroblasts, and its association with extracellular matrix (ECM)-related gene expression in idiopathic pulmonary fibrosis (IPF).

Two publicly available scRNA-seq datasets were analyzed:

- **GSE136831**: scRNA-seq data from control and IPF lung tissues.
- **GSE122960**: Publicly available scRNA-seq data from human lung tissue.

## Software environment

- R software (v4.2.1)
- Seurat package (v4.4.0)

Other R packages required for the analyses are specified in the corresponding R scripts.

## Datasets and analysis workflows

### 1. GSE136831

- GEO accession: [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831)
- Study groups: Control and IPF

**Analysis workflow:**

1. Creation of individual Seurat objects from the original scRNA-seq data.
2. Quality control and preprocessing.
3. Integration of Control and IPF datasets.
4. Cell-type annotation and visualization.
5. Analysis of TGFB1I1 (Hic-5) expression in lung cell populations.
6. Characterization of fibroblast populations.
7. Analysis of ECM-related gene expression in fibroblasts.

The R scripts include the procedures for Seurat object generation and integration, as well as downstream analyses.

### 2. GSE122960

- GEO accession: [GSE122960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122960)

**Analysis workflow:**

1. Loading a previously processed and annotated Seurat object.
2. Visualization of cell populations using UMAP.
3. Comparison of TGFB1I1 expression among cell types.
4. Visualization of TGFB1I1 expression using violin and feature plots.
5. Extraction and analysis of fibroblast populations.
6. Comparison of gene expression between TGFB1I1-positive and TGFB1I1-negative fibroblasts.
7. Analysis of ECM-related genes, including ACTA2, COL1A1, COL1A2, and CTHRC1.

Cell-type annotations generated using Sc-type in the previously processed Seurat object were used for downstream analyses.

**Required input:**

- A preprocessed Seurat object (`integrated_seurat.rds`).

This Seurat object is not included in this repository and must be obtained separately before running the analysis scripts.

## Input data

The original sequencing datasets are publicly available through the NCBI Gene Expression Omnibus (GEO).

For GSE136831, the analysis scripts include the procedures for generating and integrating Seurat objects from the original data.

For GSE122960, the previously processed Seurat object must be prepared separately, as described above.

## Notes

- The analyses were performed using publicly available scRNA-seq datasets.
- Cell-type annotations and downstream analyses were based on the corresponding Seurat objects.
- The R scripts are provided to document the analyses and facilitate reproducibility.
- Results may vary depending on the versions of R packages and the input data.

## Citation

Please cite the original publications associated with the corresponding GEO accession numbers when using these datasets.

