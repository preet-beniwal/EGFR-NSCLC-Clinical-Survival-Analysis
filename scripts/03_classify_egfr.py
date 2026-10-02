#!/usr/bin/env python3
"""
Phase 4.3 — Classify EGFR mutations into functional classes.

The clinical significance of EGFR in LUAD depends on WHICH mutation a
patient carries, not just whether EGFR is mutated. Pooling Ex19del,
L858R, Ex20ins, T790M, and everything else into one "mutant" group
washes out any signal — this is the core hypothesis of the rebuild.

Classification rules (based on protein change strings):
    Ex19del     : contains 'del' AND targets codons E746-A750 region
    L858R       : exact 'L858R'
    Ex20ins     : contains 'ins' AND targets codons A763-P772 region
    T790M       : exact 'T790M'
    G719X       : 'G719A', 'G719C', 'G719D', 'G719S'
    L861Q       : 'L861Q'
    S768I       : 'S768I'
    Other       : any other protein change

A patient may carry more than one EGFR mutation (e.g., L858R + T790M).
We record all classes and pick a primary class based on precedence:
    T790M > Ex20ins > Ex19del > L858R > uncommon > other

Output: data/egfr_variant_classes.tsv (one row per patient)
"""

import os
import re
import pandas as pd

IN_FILE  = "data/mutations_EGFR.tsv"
OUT_FILE = "data/egfr_variant_classes.tsv"


# ---- classification rules ----

# Exon 19 deletions cluster in codons 746-751. Regex catches E746_A750del,
# L747_P753del, E746_T751delinsA, etc.
EX19DEL_RE = re.compile(r"^(E746|L747|A750|T751|E749|P753)[_A-Za-z0-9]*(del|delins)")

# Exon 20 insertions cluster in codons 763-774. Regex catches
# A763_Y764insFQEA, D770_N771insNPG, N771_P772insH, V769_D770insASV, etc.
EX20INS_RE = re.compile(r"^(A763|Y764|V765|V769|D770|N771|P772|H773|V774)[_A-Za-z0-9]*ins")


def classify(pc):
    """Return the functional class for a single protein change string."""
    if pd.isna(pc):
        return "Unknown"
    pc = str(pc).strip()

    if pc == "L858R":
        return "L858R"
    if pc == "T790M":
        return "T790M"
    if EX19DEL_RE.match(pc):
        return "Exon19del"
    if EX20INS_RE.match(pc):
        return "Exon20ins"
    if pc in ("G719A", "G719C", "G719D", "G719S"):
        return "G719X"
    if pc == "L861Q":
        return "L861Q"
    if pc == "S768I":
        return "S768I"
    return "Other"


# precedence for choosing a patient's primary class if they carry
# more than one EGFR mutation
PRECEDENCE = ["T790M", "Exon20ins", "Exon19del", "L858R",
              "G719X", "L861Q", "S768I", "Other", "Unknown"]


def pick_primary(classes):
    for cls in PRECEDENCE:
        if cls in classes:
            return cls
    return "Unknown"


def main():
    print("[03] Classifying EGFR mutations...\n")
    df = pd.read_csv(IN_FILE, sep="\t")
    print(f"  loaded {len(df)} mutation records")
    print(f"  unique patients: {df['patientId'].nunique()}")

    df["egfr_class"] = df["proteinChange"].apply(classify)

    # Print the class distribution at the mutation level
    print("\n  Mutation-level class distribution:")
    for cls, n in df["egfr_class"].value_counts().items():
        print(f"    {cls:<12} {n:>3}")

    # Collapse to one row per patient, aggregating classes
    patient_classes = (
        df.groupby("patientId")["egfr_class"]
        .apply(lambda s: sorted(set(s)))
        .reset_index()
        .rename(columns={"egfr_class": "egfr_classes"})
    )
    patient_classes["egfr_primary_class"] = patient_classes["egfr_classes"].apply(pick_primary)
    patient_classes["egfr_classes"] = patient_classes["egfr_classes"].apply(
        lambda lst: ",".join(lst)
    )

    print("\n  Patient-level primary class distribution:")
    for cls, n in patient_classes["egfr_primary_class"].value_counts().items():
        print(f"    {cls:<12} {n:>3}")

    patient_classes.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"\n  wrote {OUT_FILE}")
    print(f"  total EGFR-mutant patients: {len(patient_classes)}")


if __name__ == "__main__":
    main()
