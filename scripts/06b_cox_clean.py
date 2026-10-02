#!/usr/bin/env python3
"""Phase 5.2b — Clean Cox model on EGFR-mutant patients.

Drops covariates with too few events to estimate (STK11 n=1, NF1 n=1).
Reports adjusted HRs for TP53 and KEAP1 only, with a caveat that KEAP1
has only 4 events and wide CIs.
"""
import os
import pandas as pd
from lifelines import CoxPHFitter

COHORT = "data/cohort_master.csv"
OUT    = "results/cox_egfr_mutant_adjusted_clean.tsv"

df = pd.read_csv(COHORT)
mut = df[df["egfr_status"] == 1].dropna(subset=["OS_MONTHS", "OS_event"]).copy()

print(f"EGFR-mutant with OS data: {len(mut)}")
print("Covariate counts:")
for g in ["TP53", "KEAP1", "STK11", "NF1"]:
    print(f"  mut_{g}: {int(mut[f'mut_{g}'].sum())}")

# Build the clean model matrix
cols = ["OS_MONTHS", "OS_event", "AGE", "SEX", "stage_group",
        "mut_TP53", "mut_KEAP1"]
d = mut[cols].dropna().copy()

for col in ["SEX", "stage_group"]:
    dummies = pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=int)
    d = pd.concat([d.drop(columns=[col]), dummies], axis=1)

cph = CoxPHFitter()
cph.fit(d, duration_col="OS_MONTHS", event_col="OS_event")
summary = cph.summary[["exp(coef)", "exp(coef) lower 95%",
                       "exp(coef) upper 95%", "p"]].copy()
summary.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
summary = summary.round(4)
print()
print(summary.to_string())
summary.to_csv(OUT, sep="\t")
print(f"\nwrote {OUT}")
