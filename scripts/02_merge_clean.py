import pandas as pd
import requests
import os

os.makedirs('spss_output', exist_ok=True)

print("Fetching clinical data via cBioPortal API...")
clinical_url = "https://www.cbioportal.org/api/studies/luad_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY"
clinical_response = requests.get(clinical_url)
clinical_data = pd.DataFrame(clinical_response.json())
clinical_pivot = clinical_data.pivot(index='patientId', columns='clinicalAttributeId', values='value').reset_index()

print("Fetching EGFR mutation data via cBioPortal API...")
# Fixed: Changed 'entrezGeneIds' to 'entrezGeneId'
mutation_url = "https://www.cbioportal.org/api/molecular-profiles/luad_tcga_pan_can_atlas_2018_mutations/mutations?sampleListId=luad_tcga_pan_can_atlas_2018_all&projection=ID&entrezGeneId=1956"
mutation_response = requests.get(mutation_url)

try:
    mutation_json = mutation_response.json()
    if isinstance(mutation_json, dict): 
        print("API error:", mutation_json)
        mutation_json = [] 
except Exception as e:
    print("JSON parse error:", e)
    mutation_json = []

mutations = pd.DataFrame(mutation_json)

print("Processing EGFR status...")
if not mutations.empty and 'patientId' in mutations.columns:
    egfr_patients = mutations['patientId'].unique()
    print(f"Found {len(egfr_patients)} patients with EGFR mutations.")
else:
    print("Warning: No EGFR mutations found or API error.")
    egfr_patients = []
    
clinical_pivot['EGFR_Status'] = clinical_pivot['patientId'].apply(lambda x: 1 if x in egfr_patients else 0)

print("Cleaning survival data...")
clinical_pivot['OS_MONTHS'] = pd.to_numeric(clinical_pivot['OS_MONTHS'], errors='coerce')
clinical_pivot['OS_STATUS_NUM'] = clinical_pivot['OS_STATUS'].apply(lambda x: 1 if '1:DECEASED' in str(x) else 0)

spss_data = clinical_pivot[['patientId', 'AGE', 'SEX', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'OS_MONTHS', 'OS_STATUS_NUM', 'EGFR_Status']].dropna(subset=['OS_MONTHS', 'OS_STATUS_NUM'])
spss_data.rename(columns={'patientId': 'PATIENT_ID'}, inplace=True)

spss_data.to_csv('spss_output/EGFR_NSCLC_SPSS_Ready.csv', index=False)
print(f"Cohort N = {len(spss_data)}")
print("Data exported successfully to spss_output/EGFR_NSCLC_SPSS_Ready.csv")
