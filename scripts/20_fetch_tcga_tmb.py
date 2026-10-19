#!/usr/bin/env python3
"""
Fetch SAMPLE-level TMB for TCGA-LUAD and merge into cohort_master.csv.
"""
import requests
import pandas as pd

STUDY = "luad_tcga_pan_can_atlas_2018"
BASE = "https://www.cbioportal.org/api"

print("Fetching TCGA SAMPLE-level clinical data...")
r = requests.get(f"{BASE}/studies/{STUDY}/clinical-data", 
                 params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"}, timeout=60)
df = pd.DataFrame(r.json())

# Pivot to wide format
df = df.pivot(index="sampleId", columns="clinicalAttributeId", values="value").reset_index()

# Find the TMB column
tmb_col = next((c for c in df.columns if "TMB" in c.upper()), None)
if not tmb_col:
    print(f"Available sample columns: {list(df.columns)}")
    raise ValueError("No TMB column found.")

print(f"Found TMB column: {tmb_col}")

# Extract patientId from sampleId (e.g., TCGA-05-4244-01 -> TCGA-05-4244)
df["patientId"] = df["sampleId"].str[:12]

df = df[["patientId", tmb_col]].rename(columns={tmb_col: "TMB"})
df["TMB"] = pd.to_numeric(df["TMB"], errors="coerce")

# Merge into cohort_master.csv
cohort = pd.read_csv("data/cohort_master.csv")
cohort = cohort.merge(df, on="patientId", how="left")
cohort.to_csv("data/cohort_master.csv", index=False)
print(f"Merged TMB. Valid TMB values: {cohort['TMB'].notna().sum()} out of {len(cohort)} patients.")
