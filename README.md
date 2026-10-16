# Statistical Analysis of Gout and Septic Arthritis Gene Expression

## Overview

This project investigates differences in clinical measurements and blood gene-expression patterns between patients with gout, septic arthritis (SA), and healthy controls (HC).

Since gout and septic arthritis can present with similar clinical symptoms, the aim was to explore whether gene-expression data could help distinguish the two conditions.

## Dataset

The analysis uses supplied whole-blood RNA-seq expression data, clinical information, and differential expression results.

- **27 samples:** 9 HC, 9 GOUT, and 9 SA
- **29,814 genes**
- **Clinical variables:** Sex and neutrophil counts

The original dataset labels septic arthritis samples as `SEPSIS`, which is retained in the code.

## Methods

Statistical analyses were performed in R, including:

- Descriptive statistics and clinical group comparisons
- Fisher's exact test and one-way ANOVA
- Principal Component Analysis (PCA)
- Exploration of supplied differential expression results
- Welch's t-tests and Pearson correlation analysis
- Benjamini-Hochberg multiple-testing correction

## Results

Neutrophil counts differed significantly across groups (ANOVA p = 3.25e-05), with the highest mean observed in septic arthritis. Sex distribution showed no significant difference between groups (Fisher's exact test p = 1.0).

In the revised exploratory GOUT vs SA comparison, **11,041 genes** met the BH-adjusted p-value threshold of 0.05, suggesting potentially distinguishable expression patterns despite the clinical similarities.

However, the small sample size and uncertainty around expression normalisation mean that further analysis and validation are needed before drawing diagnostic conclusions.

## Repository Structure

- `scripts/` — R analysis code
- `results/` — PCA and gene-expression visualisations

## Requirements

R with the `ggplot2` and `tibble` packages.

The original datasets are not included in this repository. To reproduce the analysis, provide the five input tables referenced in the script and run:

`DATA_DIR="/path/to/data" Rscript scripts/clinical_gene_expression_analysis.R`