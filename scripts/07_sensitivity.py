#!/usr/bin/env python3
"""Phase 6 — Sensitivity analysis for the EGFR class effect on OS.

Tests whether the stratified finding (Rare class HR 3-5 vs Wildtype)
is robust to:
    6.1  bootstrap resampling (1000 iter, 95% CI for Rare HR)
    6.2  restricting to Stage III/IV
    6.3  different adjustor sets (unadjusted, +age, +age+sex, +age+sex+stage)
    6.4  excluding short-term deaths (<3 months follow-up)
    6.5  summary table

Outputs:
    results/sensitivity_bootstrap.tsv
    results/sensitivity_stage_restricted.tsv
    results/sensitivity_adjustors.tsv
    results/sensitivity_excl_early.tsv
    results/sensitivity_summary.tsv
"""

import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

COHORT  = "data/cohort_master.csv"
RES_DIR = "results"
os.makedirs(RES_DIR, exist_ok=True)
np.random.seed(42)

CLASS_ORDER = ["Wildtype", "Exon19del", "L858R", "Other", "Rare"]


def build_design(d, adjustors):
    """Return a Cox-ready DataFrame with Wildtype as reference level."""
    cols = ["OS_MONTHS", "OS_event", "egfr_class_simplified"] + adjustors
    cols = [c for c in cols if c in d.columns]
    dd = d[cols].dropna().copy()

    dd["egfr_class_simplified"] = pd.Categorical(
        dd["egfr_class_simplified"], categories=CLASS_ORDER
    )
    dummies = pd.get_dummies(dd["egfr_class_simplified"],
                             prefix="egfr", drop_first=True, dtype=int)
    dd = pd.concat([dd.drop(columns=["egfr_class_simplified"]), dummies], axis=1)

    for col in [c for c in adjustors if c in ("SEX", "stage_group")]:
        dm = pd.get_dummies(dd[col], prefix=col, drop_first=True, dtype=int)
        dd = pd.concat([dd.drop(columns=[col]), dm], axis=1)

    return dd


def fit_cox(dd):
    cph = CoxPHFitter()
    cph.fit(dd, duration_col="OS_MONTHS", event_col="OS_event")
    s = cph.summary[["exp(coef)", "exp(coef) lower 95%",
                     "exp(coef) upper 95%", "p"]].copy()
    s.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
    return s.round(4)


def main():
    print("[07] Sensitivity analysis for EGFR class effect on OS\n")
    df = pd.read_csv(COHORT)
    d_os = df.dropna(subset=["OS_MONTHS", "OS_event"]).copy()
    print(f"  base cohort: {len(d_os)} patients with OS data\n")

    base_adjustors = ["AGE", "SEX", "stage_group"]

    # ---- 6.1 Bootstrap for Rare HR ----
    print("  [6.1] Bootstrap (1000 iterations) for Rare class HR")
    rare_hrs = []
    design = build_design(d_os, base_adjustors)
    n = len(design)
    for i in range(1000):
        sample = design.sample(n=n, replace=True, random_state=i)
        try:
            cph = CoxPHFitter()
            cph.fit(sample, duration_col="OS_MONTHS", event_col="OS_event")
            if "egfr_Rare" in cph.summary.index:
                rare_hrs.append(cph.summary.loc["egfr_Rare", "exp(coef)"])
        except Exception:
            continue
    rare_hrs = np.array(rare_hrs)
    boot_summary = pd.DataFrame([{
        "metric": "Rare class HR (OS, adjusted)",
        "point_estimate": round(float(np.median(rare_hrs)), 3),
        "ci_lower_2.5pct": round(float(np.percentile(rare_hrs, 2.5)), 3),
        "ci_upper_97.5pct": round(float(np.percentile(rare_hrs, 97.5)), 3),
        "n_bootstrap_iterations": len(rare_hrs),
    }])
    print(boot_summary.to_string(index=False))
    boot_summary.to_csv(os.path.join(RES_DIR, "sensitivity_bootstrap.tsv"),
                        sep="\t", index=False)

    # ---- 6.2 Stage III/IV only ----
    print("\n  [6.2] Restricted to Stage III/IV")
    late = d_os[d_os["stage_group"] == "Late"].copy()
    print(f"        n = {len(late)}")
    if len(late) >= 20:
        try:
            s = fit_cox(build_design(late, ["AGE", "SEX"]))
            print(s.loc[[i for i in s.index if i.startswith("egfr_")]].to_string())
            s.to_csv(os.path.join(RES_DIR, "sensitivity_stage_restricted.tsv"),
                     sep="\t")
        except Exception as e:
            print(f"        Cox failed: {e}")

    # ---- 6.3 Adjustor sets ----
    print("\n  [6.3] Different adjustor sets")
    rows = []
    for label, adj in [
        ("unadjusted",            []),
        ("age",                   ["AGE"]),
        ("age + sex",             ["AGE", "SEX"]),
        ("age + sex + stage",     ["AGE", "SEX", "stage_group"]),
    ]:
        try:
            s = fit_cox(build_design(d_os, adj))
            for cls in ["egfr_Rare", "egfr_L858R", "egfr_Exon19del"]:
                if cls in s.index:
                    rows.append({
                        "adjustors": label,
                        "class": cls.replace("egfr_", ""),
                        "HR": s.loc[cls, "HR"],
                        "CI_lower": s.loc[cls, "CI_lower"],
                        "CI_upper": s.loc[cls, "CI_upper"],
                        "p_value": s.loc[cls, "p_value"],
                    })
        except Exception as e:
            print(f"        {label}: failed — {e}")
    adj_df = pd.DataFrame(rows)
    print(adj_df.to_string(index=False))
    adj_df.to_csv(os.path.join(RES_DIR, "sensitivity_adjustors.tsv"),
                  sep="\t", index=False)

    # ---- 6.4 Exclude early deaths ----
    print("\n  [6.4] Excluding patients with OS_MONTHS < 3")
    excl = d_os[d_os["OS_MONTHS"] >= 3].copy()
    print(f"        n = {len(excl)} (removed {len(d_os)-len(excl)})")
    try:
        s = fit_cox(build_design(excl, base_adjustors))
        print(s.loc[[i for i in s.index if i.startswith("egfr_")]].to_string())
        s.to_csv(os.path.join(RES_DIR, "sensitivity_excl_early.tsv"), sep="\t")
    except Exception as e:
        print(f"        Cox failed: {e}")

    # ---- 6.5 Summary table ----
    print("\n  [6.5] Summary")
    print("        See individual TSVs in results/ for details.")
    print("\n[07] Done.")


if __name__ == "__main__":
    main()
