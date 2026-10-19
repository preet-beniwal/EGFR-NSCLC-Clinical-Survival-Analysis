#!/usr/bin/env python3
"""
Phase 3 Enhancement: Formal interaction test.
Combines TCGA and OncoSG, then tests if the L858R hazard ratio 
is significantly different between the two cohorts.
"""

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import os

TCGA_PATH = "data/cohort_master.csv"
ONCOSG_PATH = "data/oncosg/cohort_master.csv"

def main():
    print("[Interaction Test] Loading cohorts...")
    
    # 1. Load both cohorts
    tcga = pd.read_csv(TCGA_PATH)
    oncosg = pd.read_csv(ONCOSG_PATH)

    # 2. Standardize columns
    # TCGA uses 'stage_group', OncoSG uses 'stage_group'
    # TCGA has 'OS_event', OncoSG has 'OS_event'
    cols_needed = ["OS_MONTHS", "OS_event", "egfr_class_simplified", "AGE", "SEX", "stage_group"]
    
    tcga = tcga[cols_needed].copy()
    oncosg = oncosg[cols_needed].copy()

    tcga["cohort"] = "TCGA"
    oncosg["cohort"] = "OncoSG"
    
    # 3. Combine
    df = pd.concat([tcga, oncosg], ignore_index=True)
    df = df.dropna(subset=["OS_MONTHS", "OS_event", "egfr_class_simplified", "AGE", "SEX", "stage_group"])
    
    # 4. Create binary variables for the model
    # We want to specifically test L858R vs Wildtype (Reference)
    df = df[df["egfr_class_simplified"].isin(["L858R", "Wildtype"])].copy()
    df["is_L858R"] = (df["egfr_class_simplified"] == "L858R").astype(int)
    
    # Cohort binary: 1 for OncoSG, 0 for TCGA
    df["is_OncoSG"] = (df["cohort"] == "OncoSG").astype(int)
    
    # Sex: 1 for Male, 0 for Female
    df["is_male"] = (df["SEX"] == "Male").astype(int)
    
    # Stage: 1 for Late, 0 for Early
    df["is_late"] = (df["stage_group"] == "Late").astype(int)

    print(f"  Combined cohort: {len(df)} patients (TCGA: {len(df[df['cohort']=='TCGA'])}, OncoSG: {len(df[df['cohort']=='OncoSG'])})")
    print(f"  L858R patients: {df['is_L858R'].sum()}")

    # 5. Fit the Cox model with interaction
    # The formula 'is_L858R * is_OncoSG' expands to:
    # is_L858R + is_OncoSG + (is_L858R * is_OncoSG)
    cph = CoxPHFitter()
    cph.fit(df, duration_col="OS_MONTHS", event_col="OS_event", 
            formula="is_L858R * is_OncoSG + AGE + is_male + is_late")
    
    cph.print_summary()
    
    # 6. Extract the interaction p-value
    interaction_term = "is_L858R:is_OncoSG"
    summary = cph.summary
    if interaction_term in summary.index:
        pval = summary.loc[interaction_term, "p"]
        hr = summary.loc[interaction_term, "exp(coef)"]
        print("\n" + "="*50)
        print("INTERACTION TEST RESULT")
        print(f"  Interaction HR: {hr:.3f}")
        print(f"  Interaction p-value: {pval:.4f}")
        if pval < 0.05:
            print("SIGNIFICANT: The L858R effect is statistically different between TCGA and OncoSG.")
        else:
            print("NOT SIGNIFICANT: We cannot prove the reversal is statistically robust.")
        print("="*50)

if __name__ == "__main__":
    main()
