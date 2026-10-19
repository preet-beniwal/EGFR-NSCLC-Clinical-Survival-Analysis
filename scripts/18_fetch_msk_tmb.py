#!/usr/bin/env python3
"""
Fetch SAMPLE-level TMB data for MSK and merge into cohort_msk.csv.
"""
import requests
import pandas as pd

STUDY = "luad_mskcc_2020"
BASE = "https://www.cbioportal.org/api"

print("Fetching MSK SAMPLE-level clinical data for TMB...")
# CRITICAL FIX: clinicalDataType is now "SAMPLE" instead of "PATIENT"
r = requests.get(f"{BASE}/studies/{STUDY}/clinical-data", 
                 params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"}, timeout=60)
df = pd.DataFrame(r.json())

# Pivot to wide format (index by patientId if available, else sampleId)
index_col = "patientId" if "patientId" in df.columns else "sampleId"
df = df.pivot(index=index_col, columns="clinicalAttributeId", values="value").reset_index()

# Find TMB column
tmb_col = next((c for c in df.columns if "TMB" in c.upper()), None)
if not tmb_col:
    print(f"Available sample columns: {list(df.columns)}")
    raise ValueError("No TMB column found in MSK samples.")

print(f"Found TMB column: {tmb_col}")

df = df[[index_col, tmb_col]].rename(columns={tmb_col: "TMB"})
df["TMB"] = pd.to_numeric(df["TMB"], errors="coerce")

# Merge into existing cohort
cohort = pd.read_csv("data/msk/cohort_msk.csv")
cohort = cohort.merge(df, left_on="patientId", right_on=index_col, how="left", suffixes=("", "_y"))

# Drop duplicate column if it was merged
if f"{index_col}_y" in cohort.columns:
    cohort = cohort.drop(columns=[f"{index_col}_y"])

cohort.to_csv("data/msk/cohort_msk.csv", index=False)
print(f"Merged TMB data. Valid TMB values: {cohort['TMB'].notna().sum()} out of {len(cohort)} patients.")
