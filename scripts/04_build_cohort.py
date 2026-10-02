#!/usr/bin/env python3
"""
Phase 4.4 — Build the master analysis cohort.

Merges three data sources into one row-per-patient table:
    1. Clinical attributes (from data/clinical_patient.tsv)
    2. Sample-level stage (from data/clinical_sample.tsv)
    3. EGFR variant class (from data/egfr_variant_classes.tsv)
    4. Co-mutation flags for TP53, KRAS, STK11, KEAP1, NF1, BRAF, PIK3CA

Output: data/cohort_master.csv — one row per patient, ready for survival
analysis in Python (lifelines), R, or SPSS.

Column selection and cleaning rationale:
    - OS_MONTHS, DFS_MONTHS, DSS_MONTHS, PFS_MONTHS are the four
      survival endpoints. We keep all four for sensitivity analysis.
    - Survival status columns are strings like '0:LIVING', '1:DECEASED'.
      Converted to numeric 0/1.
    - Stage is on the SAMPLE level, joined to patient by patientId via
      sampleId prefix match.
    - EGFR class precedence: rare > specific activating > other.
    - Co-mutation flags are binary: 1 if the patient has any mutation in
      that gene, 0 otherwise.
"""

import os
import re
import pandas as pd

CLINICAL = "data/clinical_patient.tsv"
SAMPLES  = "data/clinical_sample.tsv"
EGFR     = "data/egfr_variant_classes.tsv"
OUT      = "data/cohort_master.csv"

CO_GENES = ["TP53", "KRAS", "STK11", "KEAP1", "NF1", "BRAF", "PIK3CA"]


def parse_survival_status(series):
    """'0:LIVING' -> 0, '1:DECEASED' -> 1, NaN -> NaN."""
    def conv(x):
        if pd.isna(x):
            return None
        s = str(x)
        if s.startswith("1:") or s.startswith("1"):
            return 1
        if s.startswith("0:") or s.startswith("0"):
            return 0
        return None
    return series.apply(conv)


def main():
    print("[04] Building master cohort...\n")

    # ---- clinical ----
    clin = pd.read_csv(CLINICAL, sep="\t")
    print(f"  clinical:       {len(clin)} patients, {len(clin.columns)} columns")

    # keep only the columns we will actually use
    keep = [
        "patientId", "AGE", "SEX",
        "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "HISTORY_NE0ADJUVANT_TRTYN", "RADIATION_THERAPY",
        "OS_MONTHS", "OS_STATUS",
        "DFS_MONTHS", "DFS_STATUS",
        "DSS_MONTHS", "DSS_STATUS",
        "PFS_MONTHS", "PFS_STATUS",
        "TMB_NONSYNONYMOUS", "FRACTION_GENOME_ALTERED",
    ]
    keep = [c for c in keep if c in clin.columns]
    clin = clin[keep].copy()

    # coerce numeric
    for c in ["AGE", "OS_MONTHS", "DFS_MONTHS", "DSS_MONTHS", "PFS_MONTHS",
              "TMB_NONSYNONYMOUS", "FRACTION_GENOME_ALTERED"]:
        if c in clin.columns:
            clin[c] = pd.to_numeric(clin[c], errors="coerce")

    # convert survival status
    for c in ["OS_STATUS", "DFS_STATUS", "DSS_STATUS", "PFS_STATUS"]:
        if c in clin.columns:
            clin[c.replace("_STATUS", "_event")] = parse_survival_status(clin[c])

    # normalize stage strings: "STAGE IA" -> "IA"
    def clean_stage(x):
        if pd.isna(x):
            return None
        return re.sub(r"^STAGE\s+", "", str(x).upper()).strip()

        clin["AJCC_PATHOLOGIC_TUMOR_STAGE"] = clin["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(clean_stage)

    # stage grouping: I/II versus III/IV (common survival stratification)
    def stage_group(s):
        if pd.isna(s):
            return "Unknown"
        s = str(s)
        if s.startswith("I") and not s.startswith("III") and not s.startswith("IV"):
            return "Early"
        if s.startswith("II") and not s.startswith("III"):
            return "Early"
        if s.startswith("III") or s.startswith("IV"):
            return "Late"
        return "Unknown"
    clin["stage_group"] = clin["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(stage_group)

    # ---- EGFR class ----
    egfr = pd.read_csv(EGFR, sep="\t")
    print(f"  EGFR classes:   {len(egfr)} mutant patients")

    # collapse rare classes into 'Rare'
    rare_set = {"L861Q", "G719X", "T790M", "Exon20ins"}
    egfr["egfr_class_simplified"] = egfr["egfr_primary_class"].apply(
        lambda c: "Rare" if c in rare_set else c
    )

    clin = clin.merge(
        egfr[["patientId", "egfr_class_simplified", "egfr_primary_class"]],
        on="patientId", how="left",
    )
    # every patient without an EGFR mutation record is wildtype
    clin["egfr_class_simplified"] = clin["egfr_class_simplified"].fillna("Wildtype")
    clin["egfr_status"] = (clin["egfr_class_simplified"] != "Wildtype").astype(int)

    print("\n  EGFR class distribution after merge:")
    for cls, n in clin["egfr_class_simplified"].value_counts().items():
        print(f"    {cls:<12} {n:>4}")

    # ---- co-mutation flags ----
    print("\n  Co-mutation flags:")
    for gene in CO_GENES:
        path = f"data/mutations_{gene}.tsv"
        if not os.path.exists(path):
            print(f"    {gene}: missing file")
            continue
        m = pd.read_csv(path, sep="\t")
        mutated = set(m["patientId"].unique())
        clin[f"mut_{gene}"] = clin["patientId"].isin(mutated).astype(int)
        print(f"    {gene:<7} mutated in {clin[f'mut_{gene}'].sum():>4} / {len(clin)}")

    # ---- save ----
    clin.to_csv(OUT, index=False)
    print(f"\n  wrote {OUT}")
    print(f"  final shape: {clin.shape}")

    # quick sanity: survival coverage
    print("\n  Survival coverage (non-null endpoints):")
    for c in ["OS_MONTHS", "DFS_MONTHS", "DSS_MONTHS", "PFS_MONTHS"]:
        if c in clin.columns:
            n = clin[c].notna().sum()
            print(f"    {c:<12} {n:>4} patients with data")


if __name__ == "__main__":
    main()
