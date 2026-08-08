import csv
import os


def build_mapping_from_file(matrix_filepath):
    """
    Parses the tRNA alignment reference CSV/TSV to map
    alignment indices to (base, sequence_index).
    """
    mapping_db = {}

    if not os.path.exists(matrix_filepath):
        print(f"Error: Reference file '{matrix_filepath}' not found.")
        return None

    with open(matrix_filepath, 'r', encoding='utf-8') as f:
        # Automatically detect if it's comma or tab separated
        dialect = csv.Sniffer().sniff(f.read(2048))
        f.seek(0)
        reader = csv.reader(f, dialect)

        for row in reader:
            if not row:
                continue

            # Col 0: tRNA name, Col 1+: aligned bases
            trna_name = row[0].strip()
            alignment_columns = row[1:]

            mapping = {}
            sequence_index = 1  # 1-based index for actual residues

            for align_idx, base in enumerate(alignment_columns, start=1):
                base = base.strip()
                if base:
                    # If position is not empty, record base and its real position
                    mapping[align_idx] = (base, sequence_index)
                    sequence_index += 1
                else:
                    # Mark gap positions
                    mapping[align_idx] = (None, None)

            mapping_db[trna_name] = mapping

    return mapping_db


def convert_query_file(input_filepath, mapping_db, output_filepath):
    """
    Reads the query list (Name, AlignPos), converts positions,
    and writes results to a new file.
    """
    if not mapping_db:
        return

    results = []
    with open(input_filepath, 'r', encoding='utf-8') as f:
        # Handle simple space, tab, or comma separated inputs
        for line in f:
            parts = line.replace(',', ' ').split()
            if len(parts) < 2:
                continue

            name = parts[0].strip()
            try:
                align_pos = int(parts[1].strip())
            except ValueError:
                continue

            if name in mapping_db:
                base, seq_pos = mapping_db[name].get(align_pos, (None, "Out of range"))
                results.append([name, align_pos, base if base else "GAP", seq_pos])
            else:
                results.append([name, align_pos, "N/A", "tRNA Not Found"])

    # Write output to CSV
    with open(output_filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["tRNA", "Alignment_Pos", "Base", "Sequence_Pos"])
        writer.writerows(results)

    print(f"Conversion complete. Results saved to: {output_filepath}")


# --- CONFIGURATION ---
# Path provided for your reference matrix
matrix_ref = "/mnt/share/ueda/trna_data/ecoliTrna_paper_25.csv"

# Path to the text/csv file containing the list from your image
# (e.g., "Met 37", "Arg2 69" per line)
query_input = "/home/ueda/modpositions.tsv"

# Output destination
output_csv = "/home/ueda/converted_tRNA_positions.csv"

# --- EXECUTION ---
# 1. Load the heavy matrix once into memory
tRNA_db = build_mapping_from_file(matrix_ref)

# 2. Process the list and save results
if tRNA_db:
    convert_query_file(query_input, tRNA_db, output_csv)