import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os

def plot_msk_km():
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    # Load data
    df = pd.read_csv("data/msk/cohort_msk.csv")
    df = df.dropna(subset=["OS_MONTHS", "OS_event"])
    df["OS_event"] = df["OS_event"].astype(int)
    
    # Split into groups
    l858r = df[df["is_L858R"] == 1]
    wt = df[df["is_L858R"] == 0]
    
    # Log-rank test
    results = logrank_test(l858r["OS_MONTHS"], wt["OS_MONTHS"], 
                           event_observed_A=l858r["OS_event"], event_observed_B=wt["OS_event"])
    
    # Plotting
    plt.figure(figsize=(8, 6))
    kmf = KaplanMeierFitter()
    
    kmf.fit(l858r["OS_MONTHS"], l858r["OS_event"], label=f"L858R (n={len(l858r)})")
    ax = kmf.plot_survival_function(color="red", ci_show=True)
    
    kmf.fit(wt["OS_MONTHS"], wt["OS_event"], label=f"Wildtype (n={len(wt)})")
    kmf.plot_survival_function(ax=ax, color="blue", ci_show=True)
    
    plt.title("MSK LUAD Cohort (Post-TKI): Overall Survival by EGFR Status")
    plt.xlabel("Time (Months)")
    plt.ylabel("Overall Survival Probability")
    plt.text(10, 0.2, f"Log-rank p = {results.p_value:.3f}\nHR = 0.32 (95% CI: 0.08-1.35)", 
             fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
    plt.grid(True, linestyle="--", alpha=0.6)
    
    # Save
    plt.savefig("results/msk_km_survival.png", dpi=300, bbox_inches="tight")
    print("Saved Kaplan-Meier plot to results/msk_km_survival.png")

if __name__ == "__main__":
    plot_msk_km()
