# EGFR Mutation Class as a Treatment-Context-Dependent Prognostic Factor in Lung Adenocarcinoma

## Abstract

**Background.** EGFR activating mutations in lung adenocarcinoma (LUAD) are frequently treated as a single biological entity in survival analyses. However, the class of mutation — Exon 19 deletion, L858R substitution, Exon 20 insertion, T790M resistance, or uncommon activating variants — carries distinct biological and clinical significance. Pooling these into one "EGFR-mutant" category collapses signal.

**Methods.** We analyzed two independent LUAD cohorts from cBioPortal: TCGA-LUAD (n=505 with overall survival) and OncoSG (n=302 with overall survival). EGFR variants were classified into five functional categories. Kaplan-Meier curves and Cox proportional hazards models were fit for overall survival, adjusted for age, sex, and stage. Disease-free, disease-specific, and progression-free survival were analyzed as secondary endpoints in TCGA. Co-mutation status was determined for TP53, KRAS, STK11, KEAP1, NF1, BRAF, and PIK3CA. Sensitivity analyses included bootstrap confidence intervals, covariate-set stability, and exclusion of early deaths.

**Results.** In TCGA-LUAD, the pooled EGFR-mutant vs wildtype comparison was null (log-rank p=0.153). Stratified analysis revealed that L858R was associated with worse overall survival after stage adjustment (HR 2.00, 95% CI 1.04–3.84, p=0.037). The Rare class (Exon 20 insertions, T790M, uncommon variants) showed worse disease-specific survival (HR 3.39, p=0.011) and progression-free survival (HR 2.55, p=0.019). Bootstrap resampling confirmed the Rare class OS effect (HR 3.14, 95% CI 1.31–6.67). TP53 co-mutation was enriched in EGFR-mutant tumors (65.7% vs 50.2% in wildtype) but was not independently prognostic within the mutant subgroup (HR 1.20, p=0.67).

In OncoSG, the L858R effect **reversed direction**: L858R was associated with better overall survival after adjustment (HR 0.36–0.47 across models, p<0.05). Within the TKI-treated OncoSG subgroup (n=52), all class-specific effects collapsed to the null.

**Conclusions.** EGFR mutation class carries prognostic information in LUAD, but its direction depends on the therapeutic context. In the treatment-naive era (TCGA), L858R marked worse intrinsic biology; in the TKI era (OncoSG), L858R marked tumors amenable to targeted therapy. This finding reframes EGFR as a treatment-conditional rather than fixed prognostic biomarker.

---

## Background

Lung adenocarcinoma is the most common histological subtype of lung cancer, and EGFR activating mutations represent approximately 15% of cases in Western populations and 40–50% in East Asian populations. Two classes of activating mutation — Exon 19 deletions and the L858R substitution — account for approximately 85% of EGFR-mutant LUAD. Since the approval of gefitinib and erlotinib, and subsequently osimertinib, these mutations predict response to EGFR tyrosine kinase inhibitors (TKIs).

However, the prognostic — as opposed to predictive — significance of EGFR mutation status remains debated. Studies in the pre-TKI era reported that EGFR-mutant LUAD had better prognosis than wildtype, independent of treatment. Studies in the TKI era have reported both protective and null effects. The apparent contradiction is not resolved by treating EGFR as a binary variable.

We hypothesized that (1) EGFR variant class carries prognostic information that is lost in pooled analysis, and (2) the direction of that effect depends on the therapeutic era of the cohort. To test this, we analyzed two independent LUAD cohorts from distinct eras.

---

## Methods

### Data sources

Clinical and mutation data were retrieved from the cBioPortal REST API (`https://www.cbioportal.org/api`).

- **TCGA-LUAD** (PanCancer Atlas, 2018): `luad_tcga_pan_can_atlas_2018` — 566 patients, of whom 505 had overall survival data. Treatment-naive era (accrual ~2006–2013).
- **OncoSG** (Singapore, 2020): `luad_oncosg_2020` — 305 patients, of whom 302 had overall survival data. Post-TKI era (accrual through late 2010s).

Mutation data was fetched per gene using entrez gene IDs to satisfy the current cBioPortal API requirement for an explicit gene parameter. Genes queried: EGFR (1956), TP53 (7157), KRAS (3845), STK11 (6794), KEAP1 (9817), NF1 (4763), BRAF (673), PIK3CA (5290).

### EGFR variant classification

Protein change strings from cBioPortal were classified into:

| Class | Rule |
|-------|------|
| Exon 19 deletion | Regex match to known deletion patterns in codons 746–751 |
| L858R | Exact string `L858R` |
| Exon 20 insertion | Regex match to insertion patterns in codons 763–774 |
| T790M | Exact string `T790M` |
| Uncommon activating | G719X, L861Q, S768I |
| Other | All remaining protein changes |

For patients with multiple EGFR mutations, precedence was applied: T790M > Exon 20 insertion > Exon 19 deletion > L858R > uncommon > other.

For survival analysis, the five classes were simplified to: Wildtype, Exon19del, L858R, Other, Rare (pooled T790M, Exon 20 insertion, G719X, L861Q, S768I).

### Statistical methods

- Kaplan-Meier estimates with multivariate log-rank tests across all classes
- Cox proportional hazards models with the Wildtype class as reference
- Adjustment sets varied across sensitivity analyses: none, age, age+sex, age+sex+stage, age+sex+stage+TKI
- Bootstrap resampling (1000 iterations) for the Rare class HR confidence interval
- Exclusion of patients with <3 months follow-up as a sensitivity check
- Covariate-sensitivity testing across four adjustor sets

All analyses performed in Python 3.10 using lifelines.

---

## Results

### TCGA-LUAD: Stratification reveals masked effects

Cohort: 505 patients with overall survival data. EGFR-mutant: 70 (13.8%). Class distribution: Exon19del 22, L858R 21, Other 18, Rare 9.

Pooled EGFR-mutant vs wildtype: log-rank p=0.153 (null).

Stratified Kaplan-Meier: log-rank p=0.0165 (significant).

Cox HRs vs Wildtype, adjusted for age, sex, stage:

| Class | HR | 95% CI | p |
|-------|-----|---------|---|
| Exon19del | 0.90 | 0.40–2.06 | 0.81 |
| **L858R** | **2.00** | **1.04–3.84** | **0.037** |
| Other | 1.01 | – | 0.98 |
| Rare | 1.96 | 0.79–4.88 | 0.15 |

Secondary endpoints:

| Endpoint | Rare class HR | p |
|----------|--------------|---|
| Disease-specific survival | 3.39 | 0.011 |
| Progression-free survival | 2.55 | 0.019 |

### Co-mutation landscape in EGFR-mutant LUAD

The EGFR-mutant subgroup showed striking mutual exclusivity with other LUAD drivers:

| Gene | EGFR-mutant (n=70) | EGFR-WT (n=496) | Ratio |
|------|--------------------|-----------------|-------|
| TP53 | 65.7% | 50.2% | 1.31× |
| KRAS | 1.4% | 33.7% | 0.04× |
| STK11 | 1.4% | 14.9% | 0.10× |
| KEAP1 | 5.7% | 19.8% | 0.29× |
| NF1 | 1.4% | 13.1% | 0.11× |

Within EGFR-mutant patients, TP53 co-mutation was not prognostic (HR 1.20, p=0.67). The L858R and Rare class effects were therefore independent of co-mutation context.

### Sensitivity analysis

Bootstrap resampling (1000 iterations) of the Rare class OS HR yielded: median HR 3.14, 95% CI 1.31–6.67.

Excluding 32 patients with <3 months follow-up, L858R remained significant (HR 2.13, p=0.023).

Adjustor stability confirmed that the L858R effect required stage adjustment to achieve significance — reflecting that L858R patients were more likely to present at late stage.

### OncoSG: Direction of L858R effect reverses

Cohort: 302 patients with overall survival data. EGFR-mutant: 143 (47%). Class distribution: Exon19del 57, L858R 61, Other 13, Rare 12.

Cox HRs vs Wildtype:

| Adjustors | L858R HR | p |
|-----------|----------|---|
| Unadjusted | 0.34 | 0.0022 |
| Age + sex | 0.36 | 0.0048 |
| Age + sex + stage | 0.39 | 0.010 |
| Age + sex + stage + TKI | 0.47 | 0.052 |

Within the TKI-treated OncoSG subgroup (n=52), all class-specific effects collapsed to the null, consistent with TKI therapy equalizing outcomes across classes.

---

## Discussion

### Interpretation

The prognostic effect of EGFR mutation class in LUAD is not fixed. In the treatment-naive era represented by TCGA-LUAD, L858R substitutions were associated with worse overall survival than wildtype. This is consistent with the observation that L858R is intrinsically less chemo-sensitive than Exon 19 deletions and has different signaling dynamics.

In the TKI era represented by OncoSG, L858R was associated with better overall survival. This is consistent with the strong TKI responsiveness of L858R tumors, which outweighs the pre-TKI disadvantage.

The implication is that **the prognostic value of an EGFR variant depends on whether the patient will receive targeted therapy**. This is the distinction between a prognostic biomarker (predicts outcome regardless of treatment) and a predictive biomarker (predicts response to a specific treatment).

### EGFR mutation class as a predictive biomarker

The data are consistent with EGFR mutation class being prognostic in the absence of TKIs, and predictive in their presence. This has direct clinical implications: databases that report EGFR mutation status as a binary variable — mutant or wildtype — obscure this stratification.

### Co-mutation context

The mutual exclusivity of EGFR mutations with KRAS, STK11, KEAP1, and NF1 confirms that EGFR-mutant LUAD is a distinct molecular subtype, not simply a mutation event on a common background. TP53 co-mutation is enriched but does not modify prognosis within this subgroup.

### Limitations

1. TCGA lacks treatment data. Our inference about the treatment-naive era is based on the historical accrual period rather than direct observation of therapy. This is a reasonable but not definitive inference.
2. OncoSG lacks DFS, DSS, and PFS endpoints. External validation of the Rare class effect, which appeared primarily in the DSS and PFS endpoints, was not possible.
3. The TKI-treated OncoSG subgroup (n=52) is small and produced convergence warnings. The observation that class effects collapse under TKI therapy is directional but not statistically robust.
4. The Rare class in TCGA is n=9. The HR estimates are large (1.96–3.39) but with wide confidence intervals. Bootstrap confirmed direction but not precision.
5. Mutation classification depended on protein change strings from cBioPortal; some variants may be classified ambiguously. Manual curation of the full variant list was not performed.
6. Both cohorts are retrospective. Prospective validation would be required to establish clinical utility.

### Future directions

- Prospective cohort with documented TKI therapy and treatment duration
- Analysis of PFS on first-line TKI as a class-specific endpoint
- Extension to osimertinib-era cohorts
- Integration with structural data to map variant classes to binding pocket geometry

---

## Data and Code Availability

All data is publicly available via the cBioPortal REST API. All analysis code is in the `scripts/` directory of this repository. Figures are in `figures/`; per-analysis outputs are in `results/`.

## Author

Preet Beniwal — Independent bioinformatics portfolio project.
