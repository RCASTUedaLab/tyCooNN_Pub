import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.special import gammaln
import os
from statsmodels.stats.multitest import multipletests

# --- SETTINGS ---
INPUT_FILE = '/home/ueda/project/tycoon_pub_result/bias_corrected_counts.tsv'
OUTPUT_DIR = '/home/ueda/project/tycoon_pub_result/output_plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Dictionary for manual label adjustment
MANUAL_OFFSETS = {}


# --- CORE FUNCTIONS (English comments as requested) ---

def nb_log_likelihood(y, mu, alpha):
    """
    Direct calculation of Negative Binomial Log-Likelihood.
    Formula: Var = mu + alpha * mu^2
    This uses the Gamma-Poisson parameterization.
    """
    if mu <= 0:
        return -1e10  # Penalize non-positive means

    # size parameter r = 1/alpha
    r = 1.0 / alpha
    prob = r / (r + mu)

    # NB Log-likelihood formula using Log-Gamma for numerical stability
    ll = (gammaln(y + r) - gammaln(y + 1) - gammaln(r) +
          r * np.log(prob) + y * np.log(1 - prob + 1e-15))
    return np.sum(ll)


def estimate_common_alpha(df, all_sample_cols):
    """
    Estimate a global dispersion parameter (alpha) from the mean-variance trend.
    Calculated as the median of (Var - Mean) / Mean^2 across all features.
    """
    means = df[all_sample_cols].mean(axis=1)
    variances = df[all_sample_cols].var(axis=1)

    # Filter out features with low counts or no overdispersion
    mask = (means > 1) & (variances > means)
    if not mask.any():
        return 0.1  # Default fallback value

    m = means[mask]
    v = variances[mask]
    alphas = (v - m) / (m ** 2)
    return np.median(alphas)


def run_differential_analysis(df, ctrl_cols, expt_cols, title, common_alpha):
    """
    Perform Differential Expression Analysis using manual Likelihood Ratio Test.
    Compares H1 (Separate means for groups) vs H0 (Single mean for all samples).
    """
    results = []
    print(f"Analyzing: {title}...")

    for _, row in df.iterrows():
        t_name = row['tRNA']
        ctrl_vals = row[ctrl_cols].values.astype(float)
        expt_vals = row[expt_cols].values.astype(float)
        all_vals = np.concatenate([ctrl_vals, expt_vals])

        # 1. Full Model Likelihood (Separate means for Control and Experiment)
        mu_ctrl = np.mean(ctrl_vals)
        mu_expt = np.mean(expt_vals)
        ll_full = (nb_log_likelihood(ctrl_vals, mu_ctrl, common_alpha) +
                   nb_log_likelihood(expt_vals, mu_expt, common_alpha))

        # 2. Null Model Likelihood (Single common mean for all samples)
        mu_all = np.mean(all_vals)
        ll_null = nb_log_likelihood(all_vals, mu_all, common_alpha)

        # 3. Likelihood Ratio Test (LRT)
        # Test statistic = 2 * (LL_full - LL_null) follows Chi-square distribution
        lr_stat = 2 * (ll_full - ll_null)
        lr_stat = max(0, lr_stat)  # Guard against minor floating point precision errors
        p_val = stats.chi2.sf(lr_stat, df=1)  # df=1 for the single parameter difference

        # Calculate Log2 Fold Change
        log2FC = np.log2((mu_expt + 1e-9) / (mu_ctrl + 1e-9))

        if np.isnan(p_val): p_val = 1.0
        neg_log_p = -np.log10(p_val) if p_val > 0 else 0

        results.append({
            'tRNA': t_name,
            'log2FC': log2FC,
            'p_value': p_val,
            'neg_log_p': neg_log_p
        })

    res_df = pd.DataFrame(results)

    # --- Status Categorization ---
    num_tests = len(res_df)
    P_LIMIT = 0.01
    FC_LIMIT = 0.5
    res_df['status'] = 'Stable'
    res_df.loc[(res_df['p_value'] < P_LIMIT) & (res_df['log2FC'] > FC_LIMIT), 'status'] = 'Up'
    res_df.loc[(res_df['p_value'] < P_LIMIT) & (res_df['log2FC'] < -FC_LIMIT), 'status'] = 'Down'

    # Save Stats CSV
    safe_title = title.replace(' ', '_')
    res_df.to_csv(os.path.join(OUTPUT_DIR, f"{safe_title}_stats.csv"), index=False)

    # --- PLOTTING ---
    plt.figure(figsize=(10, 8))
    palette = {'Stable': '#95a5a6', 'Up': '#e74c3c', 'Down': '#3498db'}
    sns.scatterplot(data=res_df, x='log2FC', y='neg_log_p', hue='status',
                    palette=palette, alpha=0.7, s=90, legend=False)

    # Center the X-axis
    max_fc = res_df['log2FC'].abs().max()
    plt.xlim(-(max_fc * 1.3), max_fc * 1.3)

    # Reference Lines
    plt.axhline(-np.log10(P_LIMIT), color='black', lw=1.2, ls=':')
    plt.axvline(FC_LIMIT, color='gray', linestyle='--', lw=1.0)
    plt.axvline(-FC_LIMIT, color='gray', linestyle='--', lw=1.0)

    # Annotation logic
    for i, row in res_df.iterrows():
        if row['status'] != 'Stable':
            off_x = 8 if row['log2FC'] > 0 else -8
            off_y = 8 if i % 2 == 0 else -12
            plt.annotate(
                row['tRNA'], xy=(row['log2FC'], row['neg_log_p']),
                xytext=(off_x, off_y), textcoords='offset points',
                fontsize=9, fontweight='bold',
                ha='left' if off_x > 0 else 'right', va='bottom',
                arrowprops=dict(arrowstyle='-', color='black', lw=0.5, alpha=0.4)
            )

    plt.title(f"{title}\n(Manual NB-LRT, Alpha={common_alpha:.4f})", fontsize=16)
    plt.xlabel('log2 Fold Change', fontsize=12)
    plt.ylabel('-log10 p-value', fontsize=12)
    sns.despine()

    plt.savefig(os.path.join(OUTPUT_DIR, f"{safe_title}_volcano.png"), dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # Load data
    df = pd.read_csv(INPUT_FILE, sep='\t')

    # Identify numeric columns for global dispersion estimation
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    global_alpha = estimate_common_alpha(df, numeric_cols)
    print(f"Global Dispersion (Alpha) estimated at: {global_alpha:.6f}")

    control_samples = ["WT1", "WT2", "WT3"]
    comparisons = {
        "LOG vs Stationary": ["LOG1", "LOG2", "LOG3"],
        "LB vs M9": ["M9.1", "M9.2", "M9.3"],
        "MG vs BW": ["MNMG"],
        "WT vs thiI(KO)": ["THI1", "THI2", "THI3", "THI4"],
        "WT vs trmB(KO)": ["TRMB"],
        "WT vs truA(KO)": ["TRUA"],
        "WT vs truB(KO)": ["TRUB"],
        "WT vs TapT(KO)": ["TAPT1", "TAPT2"]
    }

    for title, expt_samples in comparisons.items():
        avail = [c for c in expt_samples if c in df.columns]
        if avail:
            run_differential_analysis(df, control_samples, avail, title, global_alpha)


if __name__ == "__main__":
    main()