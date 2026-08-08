import os
import glob
import pandas as pd


def summarize_tRNA_counts(base_dir):
    results = []
    # Compare tools: bwa and minimap2
    tools = ['bwa', 'minimap2']

    for tool in tools:
        target_path = os.path.join(base_dir, tool, "*.count")
        files = glob.glob(target_path)

        for file_path in files:
            # Get tRNA key from filename (e.g., ala1.count -> ala1)
            file_name = os.path.basename(file_path)
            tRNA_key = file_name.replace(".count", "").lower()

            correct_reads = 0
            incorrect_reads = 0
            unmapped_reads = 0

            with open(file_path, 'r') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if not cols or len(cols) < 4:
                        continue

                    ref_name = cols[0]
                    mapped_count = int(cols[2])
                    unmapped_count = int(cols[3])

                    if ref_name == '*':
                        # Handle unmapped reads from samtools idxstats format
                        unmapped_reads += unmapped_count
                    else:
                        # Check if mapped to the correct tRNA target
                        if tRNA_key in ref_name.lower():
                            correct_reads += mapped_count
                        else:
                            incorrect_reads += mapped_count

            # Initial total calculation
            mapped_total = correct_reads + incorrect_reads
            raw_total = mapped_total + unmapped_reads

            # Adjust for multi-mapping if total exceeds 1000
            # Assuming the original read count was 1000
            if raw_total > 1000:
                excess = raw_total - 1000
                # Subtract excess from Correct reads to normalize total to 1000
                correct_reads = max(0, correct_reads - excess)
                mapped_total = correct_reads + incorrect_reads
                raw_total = 1000

            # Calculate metrics
            # Read Retention Rate: mapped / total
            retention_rate = mapped_total / raw_total if raw_total > 0 else 0
            # Accuracy: correct / mapped
            accuracy = correct_reads / mapped_total if mapped_total > 0 else 0

            results.append({
                "Tool": tool,
                "File": file_name,
                "Target_tRNA": tRNA_key,
                "Correct": correct_reads,
                "Incorrect": incorrect_reads,
                "Unmapped": unmapped_reads,
                "Total": raw_total,
                "Read_Retention_Rate": retention_rate,
                "Accuracy": accuracy
            })

    return pd.DataFrame(results)


# Execution
base_directory = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/fastq_validation"
df = summarize_tRNA_counts(base_directory)

# Output results
print(df.to_string(index=False))