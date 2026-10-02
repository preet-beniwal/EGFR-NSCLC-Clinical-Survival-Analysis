# EGFR NSCLC Clinical Survival Analysis

**A treatment-context-dependent prognostic effect of EGFR mutation class in lung adenocarcinoma.**

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Key Finding

**Pooled EGFR mutation status is not prognostic in lung adenocarcinoma. Stratified by variant class and treatment context, it is.**

In 505 TCGA-LUAD patients (treatment-naive era, 2006–2013), L858R substitutions carried worse overall survival than wildtype after adjusting for age, sex, and stage (HR 2.00, 95% CI 1.04–3.84, p=0.037). In 302 patients from OncoSG (post-TKI era, up to 2020), the same comparison **reversed**: L858R was protective (HR 0.34, p=0.005). The prognostic effect of EGFR mutation class is therefore not intrinsic to the variant — it is conditional on the therapeutic context in which the patient is treated.

---

## Abstract

**Background.** EGFR activating mutations in lung adenocarcinoma (LUAD) are not a monolithic entity. Exon 19 deletions, L858R substitutions, Exon 20 insertions, and uncommon variants have different biology, different responsiveness to tyrosine kinase inhibitors (TKIs), and different clinical trajectories. Pooling them into a single "EGFR-mutant" category obscures clinically meaningful heterogeneity.

**Methods.** We analyzed two independent LUAD cohorts from cBioPortal: TCGA-LUAD (n=505, treatment-naive era) and OncoSG (n=302, post-TKI era). EGFR variants were classified into Exon 19 deletion, L858R, Exon 20 insertion, T790M, uncommon activating, and other. Kaplan-Meier and Cox proportional hazards models were fit for overall survival (OS), with disease-free (DFS), disease-specific (DSS), and progression-free (PFS) survival as secondary endpoints where available. Co-mutation status for TP53, KRAS, STK11, KEAP1, NF1, BRAF, and PIK3CA was included as effect modifiers.

**Results.** In TCGA-LUAD, the pooled EGFR-mutant vs wildtype comparison was null (log-rank p=0.153). Stratified by class, L858R was associated with worse OS after adjustment (HR 2.00, p=0.037), and the Rare class was associated with worse DSS (HR 3.39, p=0.011) and PFS (HR 2.55, p=0.019). The Rare effect was independent of TP53 co-mutation. In OncoSG, the L858R effect reversed direction (HR 0.34–0.47 across models, p<0.05). Within TKI-treated OncoSG patients (n=52), all class-specific effects collapsed.

**Conclusion.** EGFR mutation class carries prognostic information in LUAD, but its direction is treatment-context dependent. Reporting EGFR status as binary is a source of lost signal; class-stratified analysis reveals a therapeutic-eras modifier of prognosis that is invisible when patients are pooled.

---

## Project Structure

    EGFR-NSCLC-Clinical-Survival-Analysis/
    |-- data/
    |   |-- clinical_patient.tsv         # TCGA patient-level attributes
    |   |-- clinical_sample.tsv          # TCGA sample-level stage
    |   |-- mutations_*.tsv              # TCGA mutations per gene (8 files)
    |   |-- egfr_variant_classes.tsv     # TCGA classified EGFR variants
    |   |-- cohort_master.csv            # TCGA analysis-ready cohort
    |   +-- oncosg/                      # OncoSG cohort (parallel structure)
    |-- scripts/
    |   |-- 01_fetch_clinical.py         # cBioPortal API: clinical data
    |   |-- 02_fetch_mutations.py        # cBioPortal API: mutations per gene
    |   |-- 03_classify_egfr.py          # EGFR variant classification
    |   |-- 04_build_cohort.py           # TCGA cohort assembly
    |   |-- 04b_fix_stage.py             # stage string cleanup
    |   |-- 05_survival_analysis.py      # stratified KM + Cox
    |   |-- 05b_forest_plot_classes.py   # HR forest plot across endpoints
    |   |-- 06_comutation.py             # co-mutation landscape
    |   |-- 06b_cox_clean.py             # clean co-mutation Cox
    |   |-- 07_sensitivity.py            # bootstrap + covariate sensitivity
    |   |-- 09_fetch_oncosg.py           # OncoSG cohort fetch
    |   |-- 10_classify_egfr_oncosg.py   # OncoSG variant classification
    |   |-- 11_build_cohort_oncosg.py    # OncoSG cohort assembly
    |   +-- 12_survival_oncosg.py        # OncoSG survival analysis
    |-- results/
    |   |-- cox_stratified_*.tsv         # per-endpoint Cox results (TCGA)
    |   |-- comutation_frequencies.tsv   # co-mutation frequencies
    |   |-- sensitivity_*.tsv            # bootstrap + sensitivity results
    |   +-- cox_oncosg_summary.tsv       # external validation summary
    |-- figures/
    |   |-- km_egfr_class_*.png          # TCGA KM curves (4 endpoints)
    |   |-- km_egfr_tp53.png             # TP53 co-mutation KM
    |   |-- km_egfr_aggressive.png       # STK11/KEAP1 co-mutation KM
    |   |-- km_oncosg_egfr_class_os.png  # OncoSG KM
    |   +-- forest_plot_egfr_classes.png # HR forest plot
    |-- docs/
    |   +-- MANUSCRIPT.md                # full scientific writeup
    |-- LICENSE
    +-- README.md

---

## Reproducibility

All data is fetched directly from the cBioPortal REST API — no manual downloads required.

    conda create -n egfr_survival python=3.10 -y
    conda activate egfr_survival
    conda install -c conda-forge pandas requests lifelines scipy matplotlib tabulate -y

    python scripts/01_fetch_clinical.py
    python scripts/02_fetch_mutations.py
    python scripts/03_classify_egfr.py
    python scripts/04_build_cohort.py
    python scripts/04b_fix_stage.py
    python scripts/05_survival_analysis.py
    python scripts/06_comutation.py
    python scripts/06b_cox_clean.py
    python scripts/07_sensitivity.py
    python scripts/09_fetch_oncosg.py
    python scripts/10_classify_egfr_oncosg.py
    python scripts/11_build_cohort_oncosg.py
    python scripts/12_survival_oncosg.py

See `docs/MANUSCRIPT.md` for the full scientific methods and discussion.

---

## Software

Python 3.10 · pandas · numpy · requests · lifelines · scipy · matplotlib

## Data Sources

- **TCGA-LUAD** (PanCancer Atlas, 2018): `luad_tcga_pan_can_atlas_2018`
- **OncoSG** (Singapore, 2020): `luad_oncosg_2020`
- Both accessed via cBioPortal REST API (`https://www.cbioportal.org/api`)

## Author

**Preet Beniwal** — Independent bioinformatics portfolio project.
GitHub: [@preet-beniwal](https://github.com/preet-beniwal)

## License

MIT
