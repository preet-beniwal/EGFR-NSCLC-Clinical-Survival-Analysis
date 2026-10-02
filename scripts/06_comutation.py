#!/usr/bin/env python3
"""
Phase 5.2 - Co-mutation effect modification.

Within EGFR-mutant patients, does co-mutation status in TP53, STK11,
KEAP1, or KRAS modify overall survival?
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

COHORT  = "data/cohort_master.csv"
FIG_DIR = "figures"
RES_DIR = "results"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)


def km_two_groups(df, months_col, event_col, group_col, labels,
                  colors, title, out_path):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    kmf = KaplanMeierFitter()

    group0 = df[df[group_col] == 0]
    group1 = df[df[group_col] == 1]

    kmf.fit(group0[months_col], group0[event_col],
            label=f"{labels[0]} (n={len(group0)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=colors[0])

    kmf.fit(group1[months_col], group1[event_col],
            label=f"{labels[1]} (n={len(group1)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=colors[1])

    lr = logrank_test(
        group0[months_col], group1[months_col],
        event_observed_A=group0[event_col],
        event_observed_B=group1[event_col],
    )

    ax.set_title(f"{title}\nlog-rank p = {lr.p_value:.4f}", fontsize=12)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=10)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return lr.p_value


def main():
    print("[06] Co-mutation effect modification analysis\n")
    df = pd.read_csv(COHORT)
    print(f"  total cohort: {len(df)} patients")

    mut = df[df["egfr_status"] == 1].copy()
    print(f"  EGFR-mutant subgroup: {len(mut)} patients\n")

    co_genes = ["TP53", "KRAS", "STK11", "KEAP1", "NF1", "BRAF", "PIK3CA"]
    print("  Co-mutation frequencies within EGFR-mutant:")
    freq_rows = []
    for g in co_genes:
        col = f"mut_{g}"
        if col not in mut.columns:
            continue
        n = int(mut[col].sum())
        pct = 100 * n / len(mut)
        print(f"    {g:<7} {n:>3} / {len(mut)} ({pct:.1f}%)")
        freq_rows.append({"gene": g, "n_mutated": n,
                          "n_total": len(mut), "pct": round(pct, 1)})

    wt = df[df["egfr_status"] == 0].copy()
    print("\n  Co-mutation frequencies in EGFR-WT:")
    for g in co_genes:
        col = f"mut_{g}"
        if col not in wt.columns:
            continue
        n = int(wt[col].sum())
        pct = 100 * n / len(wt)
        print(f"    {g:<7} {n:>3} / {len(wt)} ({pct:.1f}%)")

    pd.DataFrame(freq_rows).to_csv(
        os.path.join(RES_DIR, "comutation_frequencies.tsv"),
        sep="\t", index=False
    )

    mut_os = mut.dropna(subset=["OS_MONTHS", "OS_event"]).copy()
    print(f"\n  EGFR-mutant with OS data: {len(mut_os)}")

    if "mut_TP53" in mut_os.columns:
        p = km_two_groups(
            mut_os, "OS_MONTHS", "OS_event", "mut_TP53",
            ["EGFR only (no TP53)", "EGFR + TP53 co-mut"],
            ["#2ca02c", "#d62728"],
            "Overall survival within EGFR-mutant patients: TP53 status",
            os.path.join(FIG_DIR, "km_egfr_tp53.png"),
        )
        print(f"\n  KM OS (EGFR-mutant, TP53 +/-): log-rank p = {p:.4f}")
        print(f"    wrote figures/km_egfr_tp53.png")

    mut_os["aggressive_co"] = (
        (mut_os.get("mut_STK11", 0) == 1) | (mut_os.get("mut_KEAP1", 0) == 1)
    ).astype(int)

    p = km_two_groups(
        mut_os, "OS_MONTHS", "OS_event", "aggressive_co",
        ["EGFR only", "EGFR + STK11/KEAP1"],
        ["#2ca02c", "#9467bd"],
        "Overall survival within EGFR-mutant patients: STK11/KEAP1 status",
        os.path.join(FIG_DIR, "km_egfr_aggressive.png"),
    )
    print(f"\n  KM OS (EGFR-mutant, STK11/KEAP1 +/-): log-rank p = {p:.4f}")
    print(f"    wrote figures/km_egfr_aggressive.png")

    print("\n  Cox model on EGFR-mutant patients only:")
    cox_cols = ["OS_MONTHS", "OS_event", "AGE", "SEX",
                "stage_group", "mut_TP53", "mut_STK11", "mut_KEAP1"]
    cox_cols = [c for c in cox_cols if c in mut_os.columns]
    d = mut_os[cox_cols].dropna().copy()

    for col in ["SEX", "stage_group"]:
        dummies = pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=int)
        d = pd.concat([d.drop(columns=[col]), dummies], axis=1)

    try:
        cph = CoxPHFitter()
        cph.fit(d, duration_col="OS_MONTHS", event_col="OS_event")
        summary = cph.summary[["exp(coef)", "exp(coef) lower 95%",
                               "exp(coef) upper 95%", "p"]].copy()
        summary.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
        summary = summary.round(4)
        out = os.path.join(RES_DIR, "cox_egfr_mutant_adjusted.tsv")
        summary.to_csv(out, sep="\t")
        print(summary.to_string())
        print(f"\n  wrote {out}")
    except Exception as e:
        print(f"  Cox failed: {e}")

    print("\n[06] Done.")


if __name__ == "__main__":
    main()
