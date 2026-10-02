#!/usr/bin/env python3
"""
Phase 5.1 — Stratified survival analysis by EGFR mutation class.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

COHORT  = "data/cohort_master.csv"
FIG_DIR = "figures"
RES_DIR = "results"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

CLASS_ORDER = ["Wildtype", "Exon19del", "L858R", "Other", "Rare"]
CLASS_COLORS = {
    "Wildtype":  "#7f7f7f",
    "Exon19del": "#1f77b4",
    "L858R":     "#d62728",
    "Other":     "#2ca02c",
    "Rare":      "#9467bd",
}

ENDPOINTS = [
    ("OS_MONTHS",  "OS_event",  "Overall survival"),
    ("DFS_MONTHS", "DFS_event", "Disease-free survival"),
    ("DSS_MONTHS", "DSS_event", "Disease-specific survival"),
    ("PFS_MONTHS", "PFS_event", "Progression-free survival"),
]


def km_plot(df, months_col, event_col, title, out_path):
    fig, ax = plt.subplots(figsize=(9, 6))
    kmf = KaplanMeierFitter()
    for cls in CLASS_ORDER:
        sub = df[df["egfr_class_simplified"] == cls]
        if len(sub) < 3:
            continue
        kmf.fit(sub[months_col], sub[event_col], label=f"{cls} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=CLASS_COLORS[cls])

    lr = multivariate_logrank_test(
        df[months_col], df["egfr_class_simplified"], df[event_col]
    )

    ax.set_title(f"{title} by EGFR mutation class (TCGA-LUAD)\n"
                 f"multivariate log-rank p = {lr.p_value:.4f}", fontsize=12)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return lr.p_value, len(df)


def cox_stratified(df, months_col, event_col):
    d = df[[months_col, event_col, "AGE", "SEX",
            "stage_group", "egfr_class_simplified"]].dropna().copy()
    d["egfr_class_simplified"] = pd.Categorical(
        d["egfr_class_simplified"], categories=CLASS_ORDER
    )
    dummies = pd.get_dummies(d["egfr_class_simplified"],
                             prefix="egfr", drop_first=True, dtype=int)
    d = pd.concat([d.drop(columns=["egfr_class_simplified"]), dummies], axis=1)

    for col in ["SEX", "stage_group"]:
        dummies = pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=int)
        d = pd.concat([d.drop(columns=[col]), dummies], axis=1)

    cph = CoxPHFitter()
    cph.fit(d, duration_col=months_col, event_col=event_col)
    summary = cph.summary[["exp(coef)", "exp(coef) lower 95%",
                           "exp(coef) upper 95%", "p"]].copy()
    summary.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
    return summary.round(4)


def main():
    print("[05] Stratified survival analysis by EGFR class\n")
    df = pd.read_csv(COHORT)
    print(f"  loaded {len(df)} patients")
    print(f"  EGFR class counts:")
    for c in CLASS_ORDER:
        print(f"    {c:<10} {int((df['egfr_class_simplified']==c).sum()):>4}")

    logrank_rows = []

    for months_col, event_col, label in ENDPOINTS:
        if months_col not in df.columns or event_col not in df.columns:
            print(f"\n  [{label}] missing columns — skipped")
            continue

        d = df.dropna(subset=[months_col, event_col]).copy()
        if len(d) < 30:
            print(f"\n  [{label}] too few patients ({len(d)}) — skipped")
            continue

        print(f"\n  [{label}]  n = {len(d)}")

        fname = f"km_egfr_class_{months_col.split('_')[0].lower()}.png"
        pval, n = km_plot(d, months_col, event_col, label,
                          os.path.join(FIG_DIR, fname))
        print(f"    log-rank p = {pval:.4f}   -> {fname}")

        logrank_rows.append({
            "endpoint": label,
            "n_patients": n,
            "multivariate_logrank_p": round(pval, 4),
        })

        # Skip Cox if any class is too small — prevents convergence failures
        class_counts = d["egfr_class_simplified"].value_counts()
        if (class_counts < 5).any():
            small = class_counts[class_counts < 5].to_dict()
            print(f"    skipping Cox: class too small {small} — KM + log-rank only")
        else:
            try:
                summary = cox_stratified(d, months_col, event_col)
                out = os.path.join(RES_DIR, f"cox_stratified_{months_col.split('_')[0]}.tsv")
                summary.to_csv(out, sep="\t")
                print(f"    Cox HRs (reference = Wildtype):")
                for row_name in summary.index:
                    if row_name.startswith("egfr_"):
                        hr = summary.loc[row_name, "HR"]
                        p  = summary.loc[row_name, "p_value"]
                        print(f"      {row_name:<20} HR={hr:<6}  p={p}")
                print(f"    wrote {out}")
            except Exception as e:
                print(f"    Cox failed: {e}")

    summary_df = pd.DataFrame(logrank_rows)
    summary_df.to_csv(os.path.join(RES_DIR, "km_logrank_summary.tsv"),
                      sep="\t", index=False)
    print(f"\n  wrote results/km_logrank_summary.tsv")
    print("\n[05] Done.")


if __name__ == "__main__":
    main()
