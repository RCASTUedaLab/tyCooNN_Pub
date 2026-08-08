import csv
import glob
import os
import pandas as pd


def export_all_rcc_data_v5(input_dir, matrix_out, metrics_out):
    # Find all CSV files
    csv_files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))

    all_counts = {}  # For Raw Count Matrix
    metrics_list = []  # For Accuracy/Throughput summary

    for file_path in csv_files:
        with open(file_path, 'r', newline='') as f:
            reader = csv.reader(f)

            # Row 1: Get ALL columns for the header
            header = next(reader, None)
            # Row 2: Ignore
            data_row = next(reader, None)
            # # Row 3: Get ALL columns for the data row
            # data_row = next(reader, None)

            if not data_row or not header:
                continue

            # Process labels: Clean up suffixes like "_ivt"
            # We keep all columns now (no [3:] slicing)
            ivt_labels = [h.replace("_ivt", "") for h in header]

            try:
                # Convert all data columns to float
                # If the first column is a string ID, we handle it via try-except
                values_for_matrix = []
                numeric_only_values = []

                for v in data_row:
                    try:
                        f_val = float(v)
                        values_for_matrix.append(f_val)
                        numeric_only_values.append(f_val)
                    except ValueError:
                        # Keep original string if it's not a number (like a label)
                        values_for_matrix.append(v)

                sample_name = os.path.basename(file_path).replace(".csv", "")

                # 1. Build dictionary for Raw Count Matrix using all columns
                all_counts[sample_name] = dict(zip(ivt_labels, values_for_matrix))

                # 2. Calculate Metrics using only the numeric parts
                if numeric_only_values:
                    row_max = max(numeric_only_values)
                    row_sum = sum(numeric_only_values)

                    metrics_list.append({
                        "Sample": sample_name,
                        "Max_Count": row_max,
                        "Total_Sum": row_sum,
                        "Accuracy": row_max / row_sum if row_sum > 0 else 0,
                        "Throughput": row_sum / 1000
                    })

            except Exception as e:
                print(f"Error in {file_path}: {e}")

    # --- 1. Export Raw Count Matrix ---
    if all_counts:
        matrix_df = pd.DataFrame.from_dict(all_counts, orient='index')

        # Sort rows and columns case-insensitively
        matrix_df = matrix_df.sort_index(axis=0, key=lambda x: x.str.lower())
        matrix_df = matrix_df.sort_index(axis=1, key=lambda x: x.str.lower())

        matrix_df.to_csv(matrix_out, index_label="RCC_Sample")
        print(f"Raw count matrix (Row 3, All Columns) saved to: {matrix_out}")

    # --- 2. Export Metrics Summary ---
    if metrics_list:
        metrics_df = pd.DataFrame(metrics_list)
        metrics_df = metrics_df.sort_values(by="Sample", key=lambda x: x.str.lower())
        metrics_df.to_csv(metrics_out, index=False)
        print(f"Metrics summary (Row 3, All Columns) saved to: {metrics_out}")


if __name__ == "__main__":
    indir = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/fast5_validation/eval_by_ivt/"
    raw_csv = "/home/ueda/raw_count_matrix.csv"
    summary_csv = "/home/ueda/metrics_summary.csv"

    export_all_rcc_data_v5(indir, raw_csv, summary_csv)

    indir = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/ivt/fast5_validation/eval_by_rcc/"
    raw_csv = "/home/ueda/raw_count_matrix2.csv"
    summary_csv = "/home/ueda/metrics_summary2.csv"

    export_all_rcc_data_v5(indir, raw_csv, summary_csv)