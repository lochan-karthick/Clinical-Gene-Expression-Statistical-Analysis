# Statistical Analysis of Gout and Septic Arthritis Gene Expression

## Overview

This project explores differences in clinical measurements and blood gene-expression patterns between patients with gout, septic arthritis (SA), and healthy controls (HC).

Gout and septic arthritis can present with similar symptoms, making them difficult to distinguish clinically. The aim was to investigate whether gene-expression data could reveal differences between the two conditions.

## Dataset

The analysis uses supplied whole-blood RNA-seq expression data, clinical information, and differential expression results.

- 27 samples: 9 HC, 9 GOUT, and 9 SA
- 29,814 genes
- Clinical variables: sex and neutrophil counts

The dataset labels septic arthritis samples as SEPSIS, which is retained in the R script.

## Methods

Statistical analyses were performed in R, including:

- Descriptive statistics and clinical group comparisons
- Fisher's exact test and one-way ANOVA
- Principal Component Analysis (PCA)
- Exploration of supplied differential expression results
- Welch's t-tests and Pearson correlation analysis
- Benjamini-Hochberg multiple-testing correction

## Results

Neutrophil counts differed significantly across the three groups (ANOVA p = 3.25e-05), while sex distribution did not (Fisher's exact test p = 1.0).

In the revised exploratory GOUT vs SA comparison, 11,041 genes met the BH-adjusted p-value threshold of 0.05.

These results are exploratory and do not establish diagnostic biomarkers. The small sample size and unverified expression normalisation method limit biological interpretation.

## Repository Structure

- `scripts/` - R analysis code
- `results/` - PCA, differential expression and gene-expression visualisations

## Requirements

R with the `ggplot2` and `tibble` packages.

The original datasets are not included. To reproduce the analysis, provide the five input tables described in the script and run:

    DATA_DIR="/path/to/data" Rscript scripts/clinical_gene_expression_analysis.R

