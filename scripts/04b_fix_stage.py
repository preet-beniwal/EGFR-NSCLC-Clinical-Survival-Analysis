#!/usr/bin/env python3
"""Phase 4.5 — One-off fix: clean stage strings and recompute stage_group.

The original cohort build applied stage_group to raw 'STAGE IB' strings
instead of cleaned 'IB' strings, producing all-Unknown. This script
loads the CSV, re-cleans the stage, recomputes stage_group, and writes
it back.
"""
import re
import pandas as pd

COHORT = "data/cohort_master.csv"

df = pd.read_csv(COHORT)
print(f"loaded {len(df)} patients")

# clean the raw stage strings
def clean_stage(x):
    if pd.isna(x):
        return None
    return re.sub(r"^STAGE\s+", "", str(x).upper()).strip()

df["AJCC_PATHOLOGIC_TUMOR_STAGE"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(clean_stage)
print("\ncleaned stage values:")
print(df["AJCC_PATHOLOGIC_TUMOR_STAGE"].value_counts(dropna=False))

# recompute stage_group on cleaned values
def stage_group(s):
    if pd.isna(s) or s is None:
        return "Unknown"
    s = str(s)
    if s.startswith("III") or s.startswith("IV"):
        return "Late"
    if s.startswith("I"):  # reaches here only if not III or IV
        return "Early"
    if s.startswith("II"):  # II or IIA/IIB but not III
        return "Early"
    return "Unknown"

df["stage_group"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(stage_group)
print("\nnew stage_group distribution:")
print(df["stage_group"].value_counts(dropna=False))

df.to_csv(COHORT, index=False)
print(f"\nwrote {COHORT}")
