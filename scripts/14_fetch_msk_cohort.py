#!/usr/bin/env python3
"""
Fetch and build the MSK LUAD (luad_mskcc_2020) cohort.
"""
import os
import requests
import pandas as pd

STUDY = "luad_mskcc_2020"
BASE = "https://www.cbioportal.org/api"
OUT_DIR = "data/msk"
os.makedirs(OUT_DIR, exist_ok=True)

def fetch(endpoint, params=None):
    url = f"{BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json()

def main():
    print("[1] Fetching clinical data...")
    clin = fetch(f"studies/{STUDY}/clinical-data", 
                 {"clinicalDataType": "PATIENT", "projection": "SUMMARY"})
    df_clin = pd.DataFrame(clin)
    
    # Pivot long to wide
    df_clin = df_clin.pivot(index="patientId", columns="clinicalAttributeId", values="value").reset_index()
    
    # --- Dynamic Column Resolution ---
    # cBioPortal schemas vary; we must map the correct MSK columns safely.
    stage_col = "AJCC_PATHOLOGIC_TUMOR_STAGE" if "AJCC_PATHOLOGIC_TUMOR_STAGE" in df_clin.columns else "PATHOLOGIC_STAGE"
    age_col = "AGE_AT_SEQ_REPORT" if "AGE_AT_SEQ_REPORT" in df_clin.columns else "AGE_AT_SURGERY"
    sex_col = "SEX" if "SEX" in df_clin.columns else "SEX_AT_BIRTH"
    
    if stage_col not in df_clin.columns:
        print(f"WARNING: Stage column not found. Available columns: {df_clin.columns.tolist()}")
    
    # Keep relevant columns
    keep = ["patientId", "OS_MONTHS", "OS_STATUS", age_col, sex_col]
    if stage_col in df_clin.columns:
        keep.append(stage_col)
        
    df_clin = df_clin[keep]
    print(f"  -> {len(df_clin)} patients")

    print("[2] Fetching EGFR mutations...")
    muts = fetch(f"molecular-profiles/{STUDY}_mutations/mutations",
                 {"sampleListId": f"{STUDY}_all", "entrezGeneId": 1956, "projection": "SUMMARY"})
    df_mut = pd.DataFrame(muts)
    print(f"  -> {len(df_mut)} EGFR mutation records")

    print("[3] Classifying variants...")
    def classify(pc):
        if pd.isna(pc): return "Wildtype"
        pc = str(pc)
        if pc == "L858R": return "L858R"
        if pc == "T790M": return "Rare"
        if pc.startswith("E746") or pc.startswith("L747") or pc.startswith("A750") or pc.startswith("T751"):
            return "Exon19del"
        if pc in ("G719A", "G719C", "G719S", "L861Q", "S768I"): return "Rare"
        if "ins" in pc: return "Rare"
        return "Other"

    df_mut["egfr_class"] = df_mut["proteinChange"].apply(classify)
    
    # Collapse to patient-level (take the first/primary class)
    patient_classes = df_mut.groupby("patientId")["egfr_class"].apply(
        lambda x: "Rare" if "Rare" in list(x) else list(x)[0]
    ).reset_index()

    print("[4] Building merged cohort...")
    df = df_clin.merge(patient_classes, on="patientId", how="left")
    df["egfr_class"] = df["egfr_class"].fillna("Wildtype")
    
    # Filter for L858R vs Wildtype for the interaction test
    df = df[df["egfr_class"].isin(["L858R", "Wildtype"])].copy()
    
    # Clean up numeric types
    df["OS_MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
    df["OS_event"] = df["OS_STATUS"].apply(lambda x: 1 if str(x).startswith("1") else (0 if str(x).startswith("0") else None))
    df["is_L858R"] = (df["egfr_class"] == "L858R").astype(int)
    
    # Simplify stage using the dynamically resolved stage_col
    def stage_group(s):
        if pd.isna(s): return "Unknown"
        s = str(s).upper()
        if s.startswith("III") or s.startswith("IV"): return "Late"
        if s.startswith("I") or s.startswith("II"): return "Early"
        return "Unknown"
        
    if stage_col in df.columns:
        df["stage_group"] = df[stage_col].apply(stage_group)
    else:
        df["stage_group"] = "Unknown"
        
    df["is_late"] = (df["stage_group"] == "Late").astype(int)
    df["is_male"] = (df[sex_col] == "Male").astype(int)
    df["AGE"] = pd.to_numeric(df[age_col], errors="coerce")

    df.to_csv(f"{OUT_DIR}/cohort_msk.csv", index=False)
    print(f"\n  Saved {len(df)} patients to data/msk/cohort_msk.csv")
    print(f"  L858R: {df['is_L858R'].sum()} | Wildtype: {(df['is_L858R']==0).sum()}")

if __name__ == "__main__":
    main()
