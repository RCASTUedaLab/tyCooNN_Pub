import os
import glob
import pandas as pd


def summarize_csv_metrics(base_dir):
    results = []
    # Search for all CSV files in the directory
    target_path = os.path.join(base_dir, "*.csv")
    files = glob.glob(target_path)

    for file_path in files:
        file_name = os.path.basename(file_path)
        # Identify the target tRNA from the filename (e.g., ala1.csv -> ala1)
        target_tRNA = file_name.replace(".csv", "").lower()

        try:
            # Read CSV: The first column is used as index
            df = pd.read_csv(file_path, index_col=0)

            if not df.empty:
                # Get ONLY the first data row (the 2nd row in the file including header)
                row = df.iloc[0].copy()

                # Normalize column names to lowercase for matching
                row.index = row.index.str.lower()

                # Identify correct mapping count
                if target_tRNA in row.index:
                    correct_reads = row[target_tRNA]
                else:
                    correct_reads = 0

                # Calculate total mapped reads
                total_mapped = row.sum()

                # Metrics Calculation
                # Read Retention Rate = total mapped / 1000
                retention_rate = total_mapped / 1000

                # Accuracy = correct / total mapped
                accuracy = correct_reads / total_mapped if total_mapped > 0 else 0

                results.append({
                    "File": file_name,
                    "Target_tRNA": target_tRNA,
                    "Correct": correct_reads,
                    "Total_Mapped": total_mapped,
                    "Read_Retention_Rate": retention_rate,
                    "Accuracy": accuracy
                })
        except Exception as e:
            print(f"Error processing {file_name}: {e}")

    return pd.DataFrame(results)


# --- Execution ---
directories = [
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/fast5_validation/eval",
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/fast5_validation/eval_by_ivt"
]

for d in directories:
    print(f"\n--- Results for: {d} ---")
    if os.path.exists(d):
        summary_df = summarize_csv_metrics(d)
        if not summary_df.empty:
            print(summary_df.to_string(index=False))
        else:
            print("No CSV files found or data is empty.")
    else:
        print(f"Directory not found: {d}")