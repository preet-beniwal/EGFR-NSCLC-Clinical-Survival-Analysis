#!/usr/bin/env python3
"""Phase 7.6 — Stratified survival analysis in OncoSG (external validation)."""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

COHORT  = "data/oncosg/cohort_master.csv"
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
    ax.set_title(f"{title} by EGFR mutation class (OncoSG)\n"
                 f"multivariate log-rank p = {lr.p_value:.4f}", fontsize=12)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return lr.p_value


def cox_model(df, adjustors):
    cols = ["OS_MONTHS", "OS_event", "egfr_class_simplified"] + adjustors
    cols = [c for c in cols if c in df.columns]
    d = df[cols].dropna().copy()
    d["egfr_class_simplified"] = pd.Categorical(
        d["egfr_class_simplified"], categories=CLASS_ORDER
    )
    dummies = pd.get_dummies(d["egfr_class_simplified"],
                             prefix="egfr", drop_first=True, dtype=int)
    d = pd.concat([d.drop(columns=["egfr_class_simplified"]), dummies], axis=1)
    for col in [c for c in adjustors if c in ("SEX", "stage_group",
                                                "TKI_TREATMENT", "SMOKING_STATUS")]:
        dm = pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=int)
        d = pd.concat([d.drop(columns=[col]), dm], axis=1)
    cph = CoxPHFitter()
    cph.fit(d, duration_col="OS_MONTHS", event_col="OS_event")
    s = cph.summary[["exp(coef)", "exp(coef) lower 95%",
                     "exp(coef) upper 95%", "p"]].copy()
    s.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
    return s.round(4)


def main():
    print("[12] OncoSG survival analysis (external validation)\n")
    df = pd.read_csv(COHORT)
    print(f"  total cohort: {len(df)}")

    d_os = df.dropna(subset=["OS_MONTHS", "OS_event"]).copy()
    print(f"  with OS data: {len(d_os)}")
    print(f"  EGFR class counts:")
    for c in CLASS_ORDER:
        n = int((d_os["egfr_class_simplified"] == c).sum())
        print(f"    {c:<10} {n:>4}")

    fname = "km_oncosg_egfr_class_os.png"
    pval = km_plot(d_os, "OS_MONTHS", "OS_event",
                   "Overall survival", os.path.join(FIG_DIR, fname))
    print(f"\n  log-rank p = {pval:.4f}  -> {fname}")

    print("\n  Cox HRs (reference = Wildtype):")
    models = [
        ("unadjusted",            []),
        ("+ age + sex",           ["AGE", "SEX"]),
        ("+ age + sex + stage",   ["AGE", "SEX", "stage_group"]),
        ("+ age + sex + stage + TKI", ["AGE", "SEX", "stage_group", "TKI_TREATMENT"]),
    ]
    rows = []
    for label, adj in models:
        try:
            s = cox_model(d_os, adj)
            print(f"\n  [{label}]")
            for cls in ["egfr_Exon19del", "egfr_L858R", "egfr_Rare"]:
                if cls in s.index:
                    hr = s.loc[cls, "HR"]
                    p  = s.loc[cls, "p_value"]
                    print(f"    {cls:<20} HR={hr:<6}  p={p}")
                    rows.append({
                        "cohort": "OncoSG", "model": label,
                        "class": cls.replace("egfr_", ""),
                        "HR": hr, "CI_lower": s.loc[cls, "CI_lower"],
                        "CI_upper": s.loc[cls, "CI_upper"], "p_value": p,
                    })
        except Exception as e:
            print(f"  [{label}] failed: {e}")

    # TKI-treated subgroup with class-size guard
    if "TKI_TREATMENT" in d_os.columns:
        print("\n  TKI-treated subgroup only:")
        tki = d_os[d_os["TKI_TREATMENT"].astype(str).str.upper().isin(
            ["YES", "Y", "1", "TRUE"])].copy()
        print(f"    n = {len(tki)} TKI-treated patients")
        print(f"    class distribution in TKI subgroup:")
        for c in CLASS_ORDER:
            n = int((tki["egfr_class_simplified"] == c).sum())
            print(f"      {c:<10} {n:>4}")
        small = [c for c in CLASS_ORDER
                 if c != "Wildtype" and (tki["egfr_class_simplified"] == c).sum() < 5]
        if small:
            print(f"    skipping Cox: classes too small {small}")
        else:
            try:
                s = cox_model(tki, ["AGE", "SEX", "stage_group"])
                print("    Cox HRs:")
                for cls in ["egfr_Exon19del", "egfr_L858R", "egfr_Rare"]:
                    if cls in s.index:
                        print(f"      {cls:<20} HR={s.loc[cls, 'HR']}  p={s.loc[cls, 'p_value']}")
            except Exception as e:
                print(f"    Cox failed: {e}")

    pd.DataFrame(rows).to_csv(
        os.path.join(RES_DIR, "cox_oncosg_summary.tsv"),
        sep="\t", index=False
    )
    print(f"\n  wrote results/cox_oncosg_summary.tsv")
    print("\n[12] Done.")


if __name__ == "__main__":
    main()
