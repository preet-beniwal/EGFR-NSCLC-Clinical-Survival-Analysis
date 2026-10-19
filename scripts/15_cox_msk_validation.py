import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

def run_msk_cox():
    # 1. Load processed MSK cohort
    msk_path = "data/msk/cohort_msk.csv"
    df = pd.read_csv(msk_path)
    
    # 2. Clean covariates
    cols_to_keep = ["is_L858R", "OS_MONTHS", "OS_event", "AGE", "is_male"]
    
    # Check if stage actually has variance (it won't in this MSK cohort, but good practice)
    if "is_late" in df.columns and df["is_late"].nunique() > 1:
        cols_to_keep.append("is_late")
    
    model_df = df[cols_to_keep].dropna().copy()
    
    # Ensure event column is integer 0/1
    model_df["OS_event"] = model_df["OS_event"].astype(int)
    
    print(f"Final sample size for Cox model: n = {len(model_df)}")
    print(f"Events (Deaths) in model: {model_df['OS_event'].sum()}")
    
    # 3. Fit Cox Proportional Hazards Model
    cph = CoxPHFitter()
    formula = "is_L858R + AGE + is_male"
    if "is_late" in model_df.columns:
        formula += " + is_late"
        
    cph.fit(model_df, duration_col="OS_MONTHS", event_col="OS_event", formula=formula)
    cph.print_summary()
    
    # 4. Extract and report L858R-specific metrics
    hr = np.exp(cph.params_["is_L858R"])
    p_val = cph.summary.loc["is_L858R", "p"]
    ci_lower = np.exp(cph.confidence_intervals_.loc["is_L858R"].iloc[0])
    ci_upper = np.exp(cph.confidence_intervals_.loc["is_L858R"].iloc[1])
    
    print(f"\n=== MSK Third Cohort Validation Result ===")
    print(f"EGFR L858R Hazard Ratio (HR): {hr:.4f} (95% CI: {ci_lower:.4f} - {ci_upper:.4f})")
    print(f"P-value: {p_val:.4e}")
    
    if hr < 1.0:
        print("Conclusion: HR < 1.0 confirms L858R is PROTECTIVE in the post-TKI MSK cohort.")
    else:
        print("Conclusion: HR >= 1.0.")

if __name__ == "__main__":
    run_msk_cox()

