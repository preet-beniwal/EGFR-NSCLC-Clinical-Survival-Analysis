#!/usr/bin/env python3
"""
Phase 7.5 — Build the OncoSG analysis cohort.

Differences from TCGA cohort build (script 04):
    - Stage column is 'STAGE' not 'AJCC_PATHOLOGIC_TUMOR_STAGE'
    - Adds TKI_TREATMENT, CHEMOTHERAPY, SMOKING_STATUS as covariates
    - 4 survival endpoints: OS, DFS, DSS, PFS (same names as TCGA)
    - Co-mutation genes are the same 7

Output: data/oncosg/cohort_master.csv
"""

import os
import re
import pandas as pd

CLINICAL = "data/oncosg/clinical_patient.tsv"
EGFR     = "data/oncosg/egfr_variant_classes.tsv"
OUT      = "data/oncosg/cohort_master.csv"

CO_GENES = ["TP53", "KRAS", "STK11", "KEAP1", "NF1", "BRAF", "PIK3CA"]


def parse_survival_status(series):
    def conv(x):
        if pd.isna(x):
            return None
        s = str(x)
        if s.startswith("1"):
            return 1
        if s.startswith("0"):
            return 0
        return None
    return series.apply(conv)


def clean_stage(x):
    if pd.isna(x):
        return None
    return re.sub(r"^STAGE\s+", "", str(x).upper()).strip()


def stage_group(s):
    if pd.isna(s) or s is None:
        return "Unknown"
    s = str(s)
    if s.startswith("III") or s.startswith("IV"):
        return "Late"
    if s.startswith("I") or s.startswith("II"):
        return "Early"
    return "Unknown"


def main():
    print("[11] Building OncoSG cohort...\n")

    clin = pd.read_csv(CLINICAL, sep="\t")
    print(f"  clinical: {len(clin)} patients, {len(clin.columns)} cols")
    print(f"  available columns: {sorted(clin.columns.tolist())}\n")

    # Clinical columns to keep (adjust for OncoSG naming)
    keep = [
        "patientId", "AGE", "SEX", "STAGE",
        "SMOKING_STATUS", "SMOKING_PACK_YEARS",
        "TKI_TREATMENT", "CHEMOTHERAPY",
        "ETHNICITY", "HISTOLOGICAL_GRADE",
        "OS_MONTHS", "OS_STATUS",
        "DFS_MONTHS", "DFS_STATUS",
        "DSS_MONTHS", "DSS_STATUS",
        "PFS_MONTHS", "PFS_STATUS",
    ]
    keep = [c for c in keep if c in clin.columns]
    clin = clin[keep].copy()

    # numeric coercion
    for c in ["AGE", "OS_MONTHS", "DFS_MONTHS", "DSS_MONTHS", "PFS_MONTHS",
              "SMOKING_PACK_YEARS"]:
        if c in clin.columns:
            clin[c] = pd.to_numeric(clin[c], errors="coerce")

    for c in ["OS_STATUS", "DFS_STATUS", "DSS_STATUS", "PFS_STATUS"]:
        if c in clin.columns:
            clin[c.replace("_STATUS", "_event")] = parse_survival_status(clin[c])

    # stage cleaning
    if "STAGE" in clin.columns:
        clin["STAGE_CLEAN"] = clin["STAGE"].apply(clean_stage)
        clin["stage_group"] = clin["STAGE_CLEAN"].apply(stage_group)
        print("  stage distribution:")
        for s, n in clin["stage_group"].value_counts().items():
            print(f"    {s:<10} {n:>4}")

    # ---- EGFR class merge ----
    egfr = pd.read_csv(EGFR, sep="\t")
    print(f"\n  EGFR classes: {len(egfr)} mutant patients")

    rare_set = {"L861Q", "G719X", "T790M", "Exon20ins", "S768I"}
    egfr["egfr_class_simplified"] = egfr["egfr_primary_class"].apply(
        lambda c: "Rare" if c in rare_set else c
    )

    clin = clin.merge(
        egfr[["patientId", "egfr_class_simplified", "egfr_primary_class"]],
        on="patientId", how="left",
    )
    clin["egfr_class_simplified"] = clin["egfr_class_simplified"].fillna("Wildtype")
    clin["egfr_status"] = (clin["egfr_class_simplified"] != "Wildtype").astype(int)

    print("\n  EGFR class distribution after merge:")
    for cls, n in clin["egfr_class_simplified"].value_counts().items():
        print(f"    {cls:<12} {n:>4}")

    # ---- co-mutations ----
    print("\n  Co-mutation flags:")
    for gene in CO_GENES:
        path = f"data/oncosg/mutations_{gene}.tsv"
        if not os.path.exists(path):
            print(f"    {gene}: missing")
            continue
        m = pd.read_csv(path, sep="\t")
        mutated = set(m["patientId"].unique())
        clin[f"mut_{gene}"] = clin["patientId"].isin(mutated).astype(int)
        print(f"    {gene:<7} mutated in {clin[f'mut_{gene}'].sum():>4} / {len(clin)}")

    # ---- save ----
    clin.to_csv(OUT, index=False)
    print(f"\n  wrote {OUT}")
    print(f"  final shape: {clin.shape}")

    print("\n  Survival coverage (non-null months):")
    for c in ["OS_MONTHS", "DFS_MONTHS", "DSS_MONTHS", "PFS_MONTHS"]:
        if c in clin.columns:
            n = clin[c].notna().sum()
            print(f"    {c:<12} {n:>4} patients with data")


if __name__ == "__main__":
    main()
