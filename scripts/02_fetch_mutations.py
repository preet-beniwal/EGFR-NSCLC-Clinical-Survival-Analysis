#!/usr/bin/env python3
"""
Phase 4.2 — Fetch mutation data from cBioPortal REST API.

The old downloader failed because cBioPortal now requires an explicit
entrezGeneId on every mutation query. The previous "fetch all mutations"
endpoint was deprecated and returns 400. This script queries one gene
at a time and fails loudly if the API misbehaves.

Genes fetched:
    EGFR (target)  +  7 co-mutation genes (TP53, KRAS, STK11, KEAP1,
    NF1, BRAF, PIK3CA)

Output: data/mutations_<GENE>.tsv  (one TSV per gene)
"""

import os
import sys
import requests
import pandas as pd

STUDY = "luad_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_all"
BASE = "https://www.cbioportal.org/api"
OUT_DIR = "data"
os.makedirs(OUT_DIR, exist_ok=True)

# gene symbol -> entrez gene ID (from NCBI)
GENES = {
    "EGFR":   1956,
    "TP53":   7157,
    "KRAS":   3845,
    "STK11":  6794,
    "KEAP1":  9817,
    "NF1":    4763,
    "BRAF":   673,
    "PIK3CA": 5290,
}


def fetch_gene_mutations(gene, entrez_id):
    """Fetch all mutations for one gene across the cohort."""
    url = f"{BASE}/molecular-profiles/{PROFILE}/mutations"
    params = {
        "sampleListId": SAMPLE_LIST,
        "entrezGeneId": entrez_id,
        "projection": "SUMMARY",
    }
    r = requests.get(url, params=params, timeout=120)
    if r.status_code != 200:
        sys.exit(f"API error {r.status_code} for {gene}\n{r.text[:400]}")
    data = r.json()
    print(f"  {gene:<7} -> {len(data):>5} mutation records  ({len(r.content):>8} bytes)")
    return data


def main():
    print("[02] Fetching mutation data from cBioPortal...\n")
    for gene, entrez in GENES.items():
        records = fetch_gene_mutations(gene, entrez)
        if not records:
            print(f"    WARNING: no mutations returned for {gene}")
            continue
        df = pd.DataFrame(records)

        # Keep only the columns we need downstream
        keep = [c for c in [
            "sampleId", "patientId", "proteinChange", "mutationType",
            "chromosome", "startPosition", "endPosition",
            "referenceAllele", "variantAllele",
            "mutationStatus", "validationStatus", "ncbiBuild",
        ] if c in df.columns]
        df = df[keep].drop_duplicates()

        out = os.path.join(OUT_DIR, f"mutations_{gene}.tsv")
        df.to_csv(out, sep="\t", index=False)
        print(f"           wrote {out}  ({len(df)} unique rows, {len(keep)} cols)")

    print("\n[02] Done.")


if __name__ == "__main__":
    main()
