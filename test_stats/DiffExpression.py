import pandas as pd
import numpy as np
import os

# --- SETTINGS ---
# Column order will follow this exact list
DIRS = [
    "WT1", "WT2", "WT3", "LOG1", "LOG2", "LOG3", "LB1", "LB2", "LB3",
    "M9.1", "M9.2", "M9.3", "MNMG", "THI1", "THI2", "THI3", "THI4",
    "TRMB", "TRUB", "TRUA", "TAPT1", "TAPT2"
]

FILENAME = 'count.csv'
COEF_FILE = '/home/ueda/project/tycoon_pub_result/coef_dat.csv'
RAW_OUTPUT = "/home/ueda/project/tycoon_pub_result/raw_counts_combined.tsv"
NORM_OUTPUT = "/home/ueda/project/tycoon_pub_result/bias_corrected_counts.tsv"
BaseD = "/mnt/share/bhaskar/rcc_data_extend/infer_50/"


def main():
    # 1. Load data from directories
    all_samples_data = []
    for dname in DIRS:
        filepath = os.path.join(BaseD + dname, FILENAME)
        if not os.path.exists(filepath):
            print(f"Skipping: {filepath} (File not found)")
            continue

        raw_counts = pd.read_csv(filepath, nrows=1)
        rcc_names = [col.replace('_rcc', '') for col in raw_counts.columns]

        res_df = pd.DataFrame({'rcc_key': rcc_names, 'count': raw_counts.iloc[0].values})
        res_df['sample_id'] = dname
        all_samples_data.append(res_df)

    if not all_samples_data:
        print("Error: No data was loaded.")
        return

    # 2. Create Raw Matrix
    full_stacked_df = pd.concat(all_samples_data, axis=0)
    raw_matrix = full_stacked_df.pivot_table(
        index='rcc_key',
        columns='sample_id',
        values='count',
        aggfunc='sum'
    ).reset_index()

    # Reorder columns to match DIRS (handle potential missing directories)
    actual_dirs = [d for d in DIRS if d in raw_matrix.columns]
    raw_matrix = raw_matrix[['rcc_key'] + actual_dirs]
    raw_matrix[actual_dirs] = raw_matrix[actual_dirs].astype(int)

    # Sort rows case-insensitively
    raw_matrix = raw_matrix.sort_values(
        by='rcc_key',
        key=lambda x: x.str.lower()
    ).reset_index(drop=True)

    # Save Raw Counts
    raw_matrix.to_csv(RAW_OUTPUT, sep='\t', index=False)
    print(f"Raw count matrix saved to: {RAW_OUTPUT}")

    # 3. Bias Correction and Column Normalization (to 100%)
    norm_coef = pd.read_csv(COEF_FILE, sep='\t')
    merged = pd.merge(raw_matrix, norm_coef, left_on='rcc_key', right_on='tRNA')

    counts_array = merged[actual_dirs].values.astype(float)
    coef_vector = merged['coef'].values

    # Step A: Bias correction
    corrected_values = counts_array / coef_vector[:, np.newaxis]

    # Step B: Total sum scaling (to 100%)
    column_sums = corrected_values.sum(axis=0)
    normalized_values = (corrected_values / column_sums[np.newaxis, :]) * 100

    # 4. Save Normalized Data
    normalized_df = pd.DataFrame(normalized_values, columns=actual_dirs)
    normalized_df.insert(0, 'tRNA', merged['tRNA'])

    # Sort rows case-insensitively
    normalized_df = normalized_df.sort_values(
        by='tRNA',
        key=lambda x: x.str.lower()
    ).reset_index(drop=True)

    # Ensure the columns follow the DIRS order one last time
    final_cols = ['tRNA'] + actual_dirs
    normalized_df = normalized_df[final_cols]

    # Save to file
    normalized_df.to_csv(NORM_OUTPUT, sep='\t', index=False, float_format='%.2f')
    print(f"Normalized matrix (Ordered by DIRS) saved to: {NORM_OUTPUT}")


if __name__ == "__main__":
    main()