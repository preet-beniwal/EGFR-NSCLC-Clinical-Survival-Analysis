#!/usr/bin/env python3
"""Phase 5.1b — Forest plot of EGFR class HRs across three endpoints."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ENDPOINTS = [
    ("OS",  "Overall survival"),
    ("DSS", "Disease-specific survival"),
    ("PFS", "Progression-free survival"),
]
CLASSES = ["egfr_Exon19del", "egfr_L858R", "egfr_Other", "egfr_Rare"]

fig, ax = plt.subplots(figsize=(9, 5))
y_labels = []
y_pos = 0

for ep, ep_label in ENDPOINTS:
    path = f"results/cox_stratified_{ep}.tsv"
    if not os.path.exists(path):
        continue
    df = pd.read_csv(path, sep="\t", index_col=0)
    for cls in CLASSES:
        if cls not in df.index:
            continue
        hr = df.loc[cls, "HR"]
        lo = df.loc[cls, "CI_lower"]
        hi = df.loc[cls, "CI_upper"]
        p  = df.loc[cls, "p_value"]

        ax.errorbar(hr, y_pos,
                    xerr=[[hr - lo], [hi - hr]],
                    fmt="o", color="black", capsize=3, markersize=6)
        marker = " *" if p < 0.05 else ""
        y_labels.append(f"{ep} | {cls.replace('egfr_', '')}{marker}")
        y_pos += 1

ax.axvline(1.0, color="red", linestyle="--", alpha=0.5, label="HR = 1 (no effect)")
ax.set_yticks(range(len(y_labels)))
ax.set_yticklabels(y_labels, fontsize=10)
ax.set_xlabel("Hazard ratio vs Wildtype (log scale)", fontsize=11)
ax.set_xscale("log")
ax.set_title("EGFR class effect on survival — Cox HR with 95% CI\n"
             "* indicates p < 0.05", fontsize=12)
ax.grid(True, linestyle="--", alpha=0.3, which="both")
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout()
plt.savefig("figures/forest_plot_egfr_classes.png", dpi=300)
print("wrote figures/forest_plot_egfr_classes.png")

