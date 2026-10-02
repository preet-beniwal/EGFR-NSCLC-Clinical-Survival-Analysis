# EGFR NSCLC Clinical Survival Analysis

## Overview
This project is a clinical validation study of EGFR mutations in Non-Small Cell Lung Cancer (NSCLC) using the TCGA-LUAD PanCancer Atlas dataset (N = 505). It serves as the clinical translation phase of a multi-project portfolio, linking computational biology findings to real-world patient outcomes.

## Methodology
1. **Data Acquisition:** Clinical and genomic data were fetched directly via the cBioPortal API.
2. **Cohort:** 505 patients with complete overall survival (OS) data. 70 patients (13.8%) harbored EGFR mutations.
3. **Statistical Analysis:** 
   - Descriptive statistics and Independent-Samples T-Test (Age vs. EGFR Status).
   - Kaplan-Meier Survival Analysis with Log-Rank test.
   - Cox Proportional Hazards Regression (adjusted for Age, Sex, and Stage).

## Results
- **Survival Analysis:** The Log-Rank test yielded a p-value of 0.1530, indicating no statistically significant difference in overall survival between EGFR-mutant and wildtype groups in this cohort.
- **Cox Regression:** EGFR mutation status was not a significant predictor of overall survival (HR = 1.164, 95% CI: 0.758 - 1.788, p = 0.487).
- **Conclusion:** High technical clustering accuracy in single-cell foundation models does not necessarily translate to statistically significant clinical survival differences in every dataset. This highlights the need for disease-specific biological interpretation (as explored in Project 2).

## Repository Structure
- `scripts/`: Python scripts for API fetching, data cleaning, and statistical analysis.
- `spss_output/`: Clean CSV datasets and formatted statistical tables ready for SPSS.
- `figures/`: High-resolution Kaplan-Meier survival curves.

## Requirements
- Python 3.10+
- pandas, requests, lifelines, scipy, matplotlib
