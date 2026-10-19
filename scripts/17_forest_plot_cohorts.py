#!/usr/bin/env python3
"""
Generate a Forest Plot comparing the L858R Hazard Ratio across TCGA, OncoSG, and MSK.
"""
import matplotlib.pyplot as plt
import numpy as np
import os

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

cohorts = [
    {
        "name": "TCGA-LUAD\n(Treatment-naive)",
        "hr": 2.23,
        "ci_lower": 1.17,
        "ci_upper": 4.29,
        "color": "#d62728" # Red
    },
    {
        "name": "OncoSG\n(Post-TKI)",
        "hr": 0.39,
        "ci_lower": 0.15,
        "ci_upper": 0.95,
        "color": "#1f77b4" # Blue
    },
    {
        "name": "MSK-IMPACT\n(Post-TKI)",
        "hr": 0.3195, 
        "ci_lower": 0.0755,
        "ci_upper": 1.3532,
        "color": "#2ca02c" # Green
    }
]

fig, ax = plt.subplots(figsize=(8, 5))

y_positions = np.arange(len(cohorts))

for i, c in enumerate(cohorts):
    # Plot the point estimate and error bars
    ax.errorbar(c["hr"], i, 
                xerr=[[c["hr"] - c["ci_lower"]], [c["ci_upper"] - c["hr"]]],
                fmt='o', color=c["color"], capsize=6, markersize=10, linewidth=2.5)
    
    # Add the HR text next to the point
    ax.text(c["ci_upper"] * 1.15, i, 
            f"{c['hr']:.2f} ({c['ci_lower']:.2f}–{c['ci_upper']:.2f})", 
            va='center', fontsize=11, fontweight='bold')

# Vertical line at HR = 1 (No effect)
ax.axvline(1.0, color="black", linestyle="--", alpha=0.6, label="HR = 1 (No effect)")

# Formatting
ax.set_yticks(y_positions)
ax.set_yticklabels([c["name"] for c in cohorts], fontsize=12)
ax.set_xscale("log")
ax.set_xlim(0.05, 10)
ax.set_xticks([0.1, 0.5, 1, 2, 5, 10])
ax.set_xticklabels(["0.1", "0.5", "1", "2", "5", "10"])

ax.set_xlabel("Hazard Ratio for Overall Survival (L858R vs. Wildtype)", fontsize=12, fontweight='bold')
ax.set_title("L858R Prognostic Effect Reverses in the TKI Era", fontsize=14, fontweight="bold")

# Add arrows to show direction
ax.annotate("← Protective", xy=(0.12, 2.35), fontsize=11, color="#2ca02c", fontweight='bold')
ax.annotate("Harmful →", xy=(2.5, 2.35), fontsize=11, color="#d62728", fontweight='bold')

ax.grid(True, axis="x", linestyle="--", alpha=0.3)
plt.tight_layout()

out_path = f"{FIG_DIR}/forest_plot_l858r_cohorts.png"
plt.savefig(out_path, dpi=300)
print(f"Wrote {out_path}")
