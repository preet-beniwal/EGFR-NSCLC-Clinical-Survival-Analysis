import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats
import os

os.makedirs('figures', exist_ok=True)

print("Loading clean data...")
df = pd.read_csv('spss_output/EGFR_NSCLC_SPSS_Ready.csv')

# Clean Stage variable for Cox Regression (e.g., "STAGE I" -> "I")
df['STAGE'] = df['AJCC_PATHOLOGIC_TUMOR_STAGE'].str.replace('STAGE ', '', regex=False)
df['STAGE'] = df['STAGE'].fillna('Unknown')

print("\n--- 1. Descriptive Statistics ---")
print(f"Total Cohort N = {len(df)}")
print(f"EGFR Mutants: {df['EGFR_Status'].sum()} | Wildtype: {len(df) - df['EGFR_Status'].sum()}")
print("\nMean Age by EGFR Status:")
print(df.groupby('EGFR_Status')['AGE'].mean())

print("\n--- 2. Independent T-Test (Age vs EGFR Status) ---")
mut_age = df[df['EGFR_Status'] == 1]['AGE']
wt_age = df[df['EGFR_Status'] == 0]['AGE']
t_stat, p_val = stats.ttest_ind(mut_age, wt_age, equal_var=False)
print(f"T-statistic: {t_stat:.4f}, P-value: {p_val:.4f}")

print("\n--- 3. Kaplan-Meier Survival Analysis ---")
kmf = KaplanMeierFitter()

# Plotting
plt.figure(figsize=(10, 6))

# Mutant group
mask_mut = (df['EGFR_Status'] == 1)
kmf.fit(df[mask_mut]['OS_MONTHS'], df[mask_mut]['OS_STATUS_NUM'], label='EGFR Mutant')
ax = kmf.plot_survival_function(ci_show=False)

# Wildtype group
mask_wt = (df['EGFR_Status'] == 0)
kmf.fit(df[mask_wt]['OS_MONTHS'], df[mask_wt]['OS_STATUS_NUM'], label='EGFR Wildtype')
kmf.plot_survival_function(ax=ax, ci_show=False)

plt.title('Overall Survival by EGFR Mutation Status (TCGA-LUAD)')
plt.xlabel('Time (Months)')
plt.ylabel('Survival Probability')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('figures/kaplan_meier_egfr.png', dpi=300)
print("Saved Kaplan-Meier plot to figures/kaplan_meier_egfr.png")

# Log-Rank Test
results = logrank_test(
    df[mask_mut]['OS_MONTHS'], df[mask_wt]['OS_MONTHS'],
    df[mask_mut]['OS_STATUS_NUM'], df[mask_wt]['OS_STATUS_NUM']
)
print(f"Log-Rank Test P-value: {results.p_value:.4f}")

print("\n--- 4. Cox Proportional Hazards Regression ---")
# Prepare data for Cox model (Drop rows with missing stage)
cox_df = df.dropna(subset=['AGE', 'SEX', 'STAGE', 'OS_MONTHS', 'OS_STATUS_NUM', 'EGFR_Status']).copy()

cph = CoxPHFitter()
# Fit model using formula to handle categorical variables (SEX, STAGE)
cph.fit(cox_df, duration_col='OS_MONTHS', event_col='OS_STATUS_NUM', 
        formula="AGE + SEX + STAGE + EGFR_Status")

print(cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']])
print("\nInterpretation: 'exp(coef)' is the Hazard Ratio (HR).")
print("HR > 1 means higher risk of death; HR < 1 means lower risk.")
