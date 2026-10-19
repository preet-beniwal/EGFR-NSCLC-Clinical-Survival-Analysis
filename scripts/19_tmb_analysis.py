#!/usr/bin/env python3
"""
Analyze TMB by L858R status. Compares TCGA and MSK.
Fixes TCGA column mapping for is_male and is_late.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from scipy import stats
import os

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("[TMB Analysis] Loading cohorts...")
    
    # 1. Load TCGA
    tcga = pd.read_csv("data/cohort_master.csv")
    tcga = tcga[tcga["egfr_class_simplified"].isin(["L858R", "Wildtype"])].copy()
    tcga["TMB"] = pd.to_numeric(tcga["TMB"], errors="coerce")
    tcga["cohort"] = "TCGA"
    
    # Create binary columns for TCGA to match MSK
    tcga["is_L858R"] = (tcga["egfr_class_simplified"] == "L858R").astype(int)
    tcga["is_male"] = (tcga["SEX"] == "Male").astype(int)
    tcga["is_late"] = (tcga["stage_group"] == "Late").astype(int)
    
    # 2. Load MSK
    msk = pd.read_csv("data/msk/cohort_msk.csv")
    msk["TMB"] = pd.to_numeric(msk["TMB"], errors="coerce")
    msk["cohort"] = "MSK"
    
    # Standardize columns for merging
    cols = ["OS_MONTHS", "OS_event", "is_L858R", "TMB", "AGE", "is_male", "is_late", "cohort"]
    
    # Drop missing data
    tcga = tcga[cols].dropna()
    msk = msk[cols].dropna()
    
    df = pd.concat([tcga, msk], ignore_index=True)
    
    # 3. Statistical test: TMB difference
    print("\n--- TMB Comparison (L858R vs Wildtype) ---")
    for c in ["TCGA", "MSK"]:
        sub = df[df["cohort"] == c]
        wt_tmb = sub[sub["is_L858R"] == 0]["TMB"]
        l858r_tmb = sub[sub["is_L858R"] == 1]["TMB"]
        
        stat, pval = stats.mannwhitneyu(wt_tmb, l858r_tmb)
        print(f"  {c}: WT median TMB = {wt_tmb.median():.2f} | L858R median TMB = {l858r_tmb.median():.2f} | p = {pval:.4e}")

    # 4. Cox Model with TMB
    print("\n--- Cox Model including TMB ---")
    cph = CoxPHFitter()
    cph.fit(df, duration_col="OS_MONTHS", event_col="OS_event", 
            formula="is_L858R + TMB + AGE + is_male + is_late + cohort")
    cph.print_summary()

    # 5. Boxplot Visualization
    fig, ax = plt.subplots(figsize=(8, 5))
    
    tcga_wt = df[(df["cohort"] == "TCGA") & (df["is_L858R"] == 0)]["TMB"]
    tcga_mut = df[(df["cohort"] == "TCGA") & (df["is_L858R"] == 1)]["TMB"]
    msk_wt = df[(df["cohort"] == "MSK") & (df["is_L858R"] == 0)]["TMB"]
    msk_mut = df[(df["cohort"] == "MSK") & (df["is_L858R"] == 1)]["TMB"]
    
    positions = [1, 2, 3.5, 4.5]
    labels = ["WT\n(TCGA)", "L858R\n(TCGA)", "WT\n(MSK)", "L858R\n(MSK)"]
    colors = ["#999999", "#d62728", "#999999", "#2ca02c"]
    
    all_vals = [tcga_wt, tcga_mut, msk_wt, msk_mut]
    bp = ax.boxplot(all_vals, positions=positions, patch_artist=True, widths=0.6)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Tumor Mutational Burden (mutations/Mb)", fontsize=12)
    ax.set_title("TMB is Lower in L858R-Mutant LUAD across Cohorts", fontsize=14, fontweight='bold')
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/tmb_l858r_comparison.png", dpi=300)
    print(f"\n  Wrote {FIG_DIR}/tmb_l858r_comparison.png")

if __name__ == "__main__":
    main()
