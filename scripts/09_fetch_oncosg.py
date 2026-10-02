#!/usr/bin/env python3
"""
Phase 7.3 — Fetch OncoSG LUAD cohort (luad_oncosg_2020) for external validation.

This is the independent cohort for reproducing the TCGA finding that
EGFR mutation class predicts survival. Different from TCGA in:
    - Geographic origin: Singapore (Asian cohort, higher EGFR mut rate)
    - Sequencing: targeted panel
    - Treatment era: post-TKI
    - Has SMOKING_STATUS, which TCGA lacks

Fetches clinical + mutations for EGFR + 7 co-mutation genes.
Output directory: data/oncosg/
"""

import os
import sys
import requests
import pandas as pd

STUDY = "luad_oncosg_2020"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_all"
BASE = "https://www.cbioportal.org/api"
OUT_DIR = "data/oncosg"
os.makedirs(OUT_DIR, exist_ok=True)

GENES = {
    "EGFR":   1956, "TP53":   7157, "KRAS":   3845, "STK11":  6794,
    "KEAP1":  9817, "NF1":    4763, "BRAF":   673,  "PIK3CA": 5290,
}


def fetch(endpoint, params=None, label="request"):
    url = f"{BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=90)
    if r.status_code != 200:
        sys.exit(f"API error {r.status_code} on {label}\n{r.text[:400]}")
    print(f"  {label:<25} {len(r.content):>8} bytes")
    return r.json()


def pivot_long_to_wide(long_df, index_col):
    return long_df.pivot(
        index=index_col, columns="clinicalAttributeId", values="value"
    ).reset_index()


def main():
    print(f"[09] Fetching OncoSG cohort ({STUDY})\n")

    print("  clinical (patient-level):")
    pat = fetch(f"studies/{STUDY}/clinical-data",
                {"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
                "patient-level")
    pat_wide = pivot_long_to_wide(pd.DataFrame(pat), "patientId")
    print(f"  -> {len(pat_wide)} patients, {len(pat_wide.columns)} cols")
    print(f"  columns: {sorted(pat_wide.columns.tolist())}")
    pat_wide.to_csv(os.path.join(OUT_DIR, "clinical_patient.tsv"),
                    sep="\t", index=False)

    print("\n  clinical (sample-level):")
    samp = fetch(f"studies/{STUDY}/clinical-data",
                 {"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
                 "sample-level")
    samp_wide = pivot_long_to_wide(pd.DataFrame(samp), "sampleId")
    print(f"  -> {len(samp_wide)} samples, {len(samp_wide.columns)} cols")
    print(f"  columns: {sorted(samp_wide.columns.tolist())}")
    samp_wide.to_csv(os.path.join(OUT_DIR, "clinical_sample.tsv"),
                     sep="\t", index=False)

    print("\n  mutations:")
    for gene, entrez in GENES.items():
        data = fetch(f"molecular-profiles/{PROFILE}/mutations",
                     {"sampleListId": SAMPLE_LIST,
                      "entrezGeneId": entrez,
                      "projection": "SUMMARY"},
                     gene)
        if not data:
            print(f"    WARNING: no mutations for {gene}")
            continue
        df = pd.DataFrame(data)
        keep = [c for c in ["sampleId", "patientId", "proteinChange",
                            "mutationType", "chromosome", "startPosition",
                            "referenceAllele", "variantAllele"]
                if c in df.columns]
        df = df[keep].drop_duplicates()
        out = os.path.join(OUT_DIR, f"mutations_{gene}.tsv")
        df.to_csv(out, sep="\t", index=False)
        print(f"    wrote {out}  ({len(df)} rows)")

    print("\n[09] Done.")


if __name__ == "__main__":
    main()
