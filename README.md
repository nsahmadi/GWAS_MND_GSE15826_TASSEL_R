# GWAS_MND_GSE15826_TASSEL_R
A pipeline for GWAS using TASSEL and R on motor neuron disease dataset (GSE15826) from the Affymetrix 500K human SNP array
# GWAS Pipeline for Motor Neuron Disease (GSE15826)

This project demonstrates a full GWAS pipeline using TASSEL and R for 164 samples from the GEO dataset GSE15826. It covers preprocessing, quality control, association mapping, and visualization.

---

## Dataset

- GEO Dataset: [GSE15826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15826)
- Platform: Affymetrix 500K Mapping Array

---

## Pipeline Steps

### 1. Preprocessing in R
- Converted GEO data to HapMap format for TASSEL input
- Generated phenotype file for binary case/control trait

### 2. Association Mapping in TASSEL
- Performed MLM (Mixed Linear Model) analysis using kinship & PCA covariates
- Bonferroni threshold: None significant
- FDR correction: 2 SNPs passed (rs6458646, rs6931955 on Chr6)

### 3. Visualization in R
- Created QQ and Manhattan plots using qqman package

---

## Scripts

- 1_format_genotype.R – Prepares HapMap file from GEO
- 2_format_phenotype.R – Builds phenotype file
- 3_plot_results.R – Generates Manhattan and QQ plots

---

## Sample Plots

Plots are available in the /plots/ folder.

---

## Tools Used

- TASSEL 5
- R (qqman, dplyr, GEOquery)
