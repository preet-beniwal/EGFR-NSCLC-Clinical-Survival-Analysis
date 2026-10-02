import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
import numpy as np

print("Loading clean data...")
df = pd.read_csv('spss_output/EGFR_NSCLC_SPSS_Ready.csv')

# Clean Stage
df['STAGE'] = df['AJCC_PATHOLOGIC_TUMOR_STAGE'].str.replace('STAGE ', '', regex=False)
df['STAGE'] = df['STAGE'].fillna('Unknown')

# --- Fix the T-Test ---
print("\n--- 2. Independent T-Test (Age vs EGFR Status) ---")
mut_age = df[df['EGFR_Status'] == 1]['AGE'].dropna()
wt_age = df[df['EGFR_Status'] == 0]['AGE'].dropna()
t_stat, p_val = stats.ttest_ind(mut_age, wt_age, equal_var=False)
print(f"T-statistic: {t_stat:.4f}, P-value: {p_val:.4f}")

# --- Generate Clean Table ---
print("\n--- 4. Cox Proportional Hazards Regression ---")
cox_df = df.dropna(subset=['AGE', 'SEX', 'STAGE', 'OS_MONTHS', 'OS_STATUS_NUM', 'EGFR_Status']).copy()
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='OS_MONTHS', event_col='OS_STATUS_NUM', formula="AGE + SEX + STAGE + EGFR_Status")

# Create a clean summary table
summary = cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
summary.columns = ['Hazard Ratio (HR)', '95% CI Lower', '95% CI Upper', 'p-value']
summary = summary.round(3)
print(summary.to_markdown())

# Save to CSV for your records
summary.to_csv('spss_output/Cox_Regression_Table.csv')
print("\nSaved clean Cox Regression table to spss_output/Cox_Regression_Table.csv")
