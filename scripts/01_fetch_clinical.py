#!/usr/bin/env python3
"""
Phase 4.1 — Fetch clinical data from cBioPortal REST API.

Two endpoints:
  PATIENT-level: age, sex, survival months, survival status, smoking
  SAMPLE-level:  stage (AJCC pathologic tumor stage lives on the sample)

Both are fetched and pivoted from long to wide format. Missing values
are preserved as NaN; filtering happens downstream in 04_build_cohort.py.
"""

import os
import sys
import requests
import pandas as pd

STUDY = "luad_tcga_pan_can_atlas_2018"
BASE = "https://www.cbioportal.org/api"
OUT_DIR = "data"
os.makedirs(OUT_DIR, exist_ok=True)


def fetch(endpoint, params=None):
    """GET request with explicit status check. Fails loudly if the API
    returns non-200 instead of silently producing an empty DataFrame."""
    url = f"{BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=60)
    if r.status_code != 200:
        sys.exit(f"API error {r.status_code} on {url}\n{r.text[:400]}")
    print(f"  fetched {len(r.content):>8} bytes  {url}")
    return r.json()


def pivot_long_to_wide(long_df, index_col):
    """cBioPortal returns clinical data as one row per attribute. Pivot
    to wide so each patient is one row with columns for each attribute."""
    return long_df.pivot(
        index=index_col,
        columns="clinicalAttributeId",
        values="value",
    ).reset_index()


def main():
    print("[01] Fetching clinical data from cBioPortal...")

    # Patient-level attributes
    print("  patient-level:")
    pat_long = fetch(
        f"studies/{STUDY}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
    )
    pat_wide = pivot_long_to_wide(pd.DataFrame(pat_long), "patientId")
    print(f"  -> {len(pat_wide)} patients, {len(pat_wide.columns)} attributes")

    # Sample-level attributes (stage)
    print("  sample-level:")
    samp_long = fetch(
        f"studies/{STUDY}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    )
    samp_wide = pivot_long_to_wide(pd.DataFrame(samp_long), "sampleId")
    print(f"  -> {len(samp_wide)} samples, {len(samp_wide.columns)} attributes")

    # Save
    pat_path = os.path.join(OUT_DIR, "clinical_patient.tsv")
    samp_path = os.path.join(OUT_DIR, "clinical_sample.tsv")
    pat_wide.to_csv(pat_path, sep="\t", index=False)
    samp_wide.to_csv(samp_path, sep="\t", index=False)
    print(f"  wrote {pat_path}")
    print(f"  wrote {samp_path}")

    # Show the columns available so we can confirm the fields we need
    print("\n  patient columns:", sorted(pat_wide.columns.tolist()))
    print("\n  sample columns:", sorted(samp_wide.columns.tolist()))


if __name__ == "__main__":
    main()
