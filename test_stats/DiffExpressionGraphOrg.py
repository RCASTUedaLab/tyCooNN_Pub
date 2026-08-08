import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# --- SETTINGS ---
INPUT_FILE = '/home/ueda/project/tycoon_pub_result/bias_corrected_counts.tsv'
OUTPUT_DIR = '/home/ueda/project/tycoon_pub_result/output_plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Using English for program comments as requested.

def run_differential_analysis(df, ctrl_cols, expt_cols, title):
    """
    Welch's T-test on log2-transformed proportional data.
    Outputs: Statistics CSV and SVG volcano plot.
    """
    results = []
    print(f"--- Processing: {title} ---")

    # Global log-SD for N=1 fallback
    all_cols = ctrl_cols + expt_cols
    global_sd_log = np.log2(df[all_cols] + 0.01).std(axis=1).mean()

    for _, row in df.iterrows():
        t_name = row['tRNA']
        ctrl_vals = row[ctrl_cols].values.astype(float)
        expt_vals = row[expt_cols].values.astype(float)
        n_ctrl, n_expt = len(ctrl_vals), len(expt_vals)

        # 1. Log2 transformation with 0.01 offset
        ctrl_log = np.log2(np.maximum(ctrl_vals, 0) + 0.01)
        expt_log = np.log2(np.maximum(expt_vals, 0) + 0.01)

        # 2. Welch's T-test
        if n_ctrl > 1 and n_expt > 1:
            _, p_val = stats.ttest_ind(expt_log, ctrl_log, equal_var=False)
        else:
            mean_c, mean_e = np.mean(ctrl_log), np.mean(expt_log)
            se = np.sqrt((global_sd_log**2 / n_ctrl) + (global_sd_log**2 / n_expt))
            z_score = (mean_e - mean_c) / (se + 1e-9)
            p_val = stats.norm.sf(abs(z_score)) * 2

        # 3. Log2 Fold Change
        mean_ctrl_raw = np.mean(ctrl_vals)
        mean_expt_raw = np.mean(expt_vals)
        log2FC = np.log2((mean_expt_raw + 1e-9) / (mean_ctrl_raw + 1e-9))

        if np.isnan(p_val) or p_val == 0: p_val = 1.0
        neg_log_p = -np.log10(p_val)

        results.append({
            'tRNA': t_name, 'log2FC': log2FC,
            'p_value': p_val, 'neg_log_p': neg_log_p
        })

    res_df = pd.DataFrame(results)

    P_LIMIT, FC_LIMIT = 0.05, 0.5
    res_df['status'] = 'Stable'
    res_df.loc[(res_df['p_value'] < P_LIMIT) & (res_df['log2FC'] > FC_LIMIT), 'status'] = 'Up'
    res_df.loc[(res_df['p_value'] < P_LIMIT) & (res_df['log2FC'] < -FC_LIMIT), 'status'] = 'Down'

    # --- SAVE CSV ---
    safe_title = title.replace(' ', '_')
    csv_path = os.path.join(OUTPUT_DIR, f"{safe_title}_stats.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"  [CSV Saved] {csv_path}")

    # --- PLOTTING ---
    plt.figure(figsize=(10, 8))
    palette = {'Stable': '#95a5a6', 'Up': '#e74c3c', 'Down': '#3498db'}
    sns.scatterplot(data=res_df, x='log2FC', y='neg_log_p', hue='status',
                    palette=palette, alpha=0.8, s=60, legend=False,
                    edgecolor='white', linewidth=0.5)

    max_fc = res_df['log2FC'].abs().max()
    plt.xlim(-(max_fc * 1.35), max_fc * 1.35)
    plt.axhline(-np.log10(P_LIMIT), color='red', lw=0.8, ls='--', alpha=0.4)
    plt.axvline(FC_LIMIT, color='gray', linestyle=':', lw=0.8)
    plt.axvline(-FC_LIMIT, color='gray', linestyle=':', lw=0.8)

    # Staggered labeling
    sig_hits = res_df[res_df['status'] != 'Stable'].sort_values('neg_log_p', ascending=False)
    for i, (_, row) in enumerate(sig_hits.iterrows()):
        off_y = 10 if i % 2 == 0 else -18
        off_x = 10 if row['log2FC'] > 0 else -10
        plt.annotate(
            row['tRNA'], xy=(row['log2FC'], row['neg_log_p']),
            xytext=(off_x, off_y), textcoords='offset points',
            fontsize=8, fontweight='bold', ha='left' if off_x > 0 else 'right',
            arrowprops=dict(arrowstyle='->', color='black', lw=0.5, alpha=0.3)
        )

    plt.title(title, fontsize=18, fontweight='bold')
    plt.xlabel('log2 Fold Change', fontsize=12)
    plt.ylabel('-log10 p-value', fontsize=12)
    sns.despine()

    # --- SAVE SVG ---
    svg_path = os.path.join(OUTPUT_DIR, f"{safe_title}_volcano.svg")
    plt.savefig(svg_path, format='svg', bbox_inches='tight')
    print(f"  [SVG Saved] {svg_path}")
    plt.close()

def main():
    df = pd.read_csv(INPUT_FILE, sep='\t')
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
            run_differential_analysis(df, control_samples, avail, title)

if __name__ == "__main__":
    main()