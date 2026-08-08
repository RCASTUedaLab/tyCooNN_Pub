import os
import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file
from tqdm import tqdm


def find_fast5_files_and_read_ids(base_dir):
    """
    Recursively finds all .fast5 files in a directory and maps read_id to file path.
    """
    print("Searching for .fast5 files and extracting read_ids...")
    fast5_map = {}
    for root, _, files in os.walk(base_dir):
        for file in tqdm(files):
            if file.endswith('.fast5'):
                fast5_path = os.path.join(root, file)
                try:
                    with get_fast5_file(fast5_path, mode='r') as f5:
                        for read in f5.get_reads():
                            read_id = read.read_id
                            if read_id:
                                fast5_map[read_id] = fast5_path
                except Exception as e:
                    print(f"Error reading {fast5_path}: {e}")
    return fast5_map


def find_pq_to_fast5_mapping(outpq_dir, fast5_map):
    """
    Maps .pq files to their corresponding .fast5 files based on read_ids.
    """
    print("Mapping .pq files to .fast5 files...")
    pq_mapping = {}
    for pq_file in tqdm(os.listdir(outpq_dir)):
        if pq_file.endswith('.pq'):
            pq_path = os.path.join(outpq_dir, pq_file)
            try:
                df = pd.read_parquet(pq_path, columns=['read_id'])
                read_ids_in_pq = set(df['read_id'].unique())

                corresponding_fast5s = set()
                for read_id in read_ids_in_pq:
                    if read_id in fast5_map:
                        corresponding_fast5s.add(fast5_map[read_id])

                pq_mapping[pq_path] = sorted(list(corresponding_fast5s))
            except Exception as e:
                print(f"Error reading {pq_path}: {e}")
    return pq_mapping


def write_mapping_to_file(mapping, output_file):
    """
    Writes the final mapping to a text file.
    """
    print(f"Writing results to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        for pq_path, fast5_list in mapping.items():
            pq_filename = os.path.basename(pq_path)
            pq_filename = pq_filename.replace(".pq","")
            f.write(f"tRNA: {pq_filename}\n")
            if fast5_list:
                for fast5_path in fast5_list:
                    fast5_filename = os.path.basename(fast5_path)
                    f.write(f"  - {fast5_filename}\n")
            else:
                f.write("  - No corresponding fast5 files found.\n")
    print("Done!")


if __name__ == "__main__":
    # --- Configuration ---
    # Be sure to set the correct paths.
    # Note: If INP_DIR contains many files, this step might take a long time.
    # To run this, you need to set one of the pairs of INP_DIR and OUTPQ_DIR.
    # For example, let's use the last one:
    OUTPQ_DIR = "/mnt/share/bhaskar/pq_binary"
    INP_DIR = "/data/suzukilab/seqdata/basecall"
    OUTPUT_FILE = "/mnt/share/bhaskar/pq_binary/pq_to_fast5_mapping_binary.txt"


    # --- End of Configuration ---

    # 1. Map all fast5 files to their read_ids
    fast5_id_to_path_map = find_fast5_files_and_read_ids(INP_DIR)
    #
    # # 2. Map pq files to their corresponding fast5 files
    pq_to_fast5_map = find_pq_to_fast5_mapping(OUTPQ_DIR, fast5_id_to_path_map)
    #
    # # 3. Write the result to a file
    write_mapping_to_file(pq_to_fast5_map, OUTPUT_FILE)

    OUTPQ_DIR = "/mnt/share/bhaskar/pq_db/ivt_8k/"
    INP_DIR = "/data/suzukilab/seqdata/basecall/split/"
    OUTPUT_FILE = "/mnt/share/bhaskar/pq_binary/pq_to_fast5_mapping_ivt.txt"
    # --- End of Configuration ---

    # 1. Map all fast5 files to their read_ids
    # fast5_id_to_path_map = find_fast5_files_and_read_ids(INP_DIR)
    #
    # # 2. Map pq files to their corresponding fast5 files
    # pq_to_fast5_map = find_pq_to_fast5_mapping(OUTPQ_DIR, fast5_id_to_path_map)
    #
    # # 3. Write the result to a file
    # write_mapping_to_file(pq_to_fast5_map, OUTPUT_FILE)

    OUTPQ_DIR = "/mnt/share/bhaskar/pq_db/rcc_12k/ec/"
    INP_DIR = "/mnt/share/suzukilab/rcc/"
    OUTPUT_FILE = "/mnt/share/bhaskar/pq_binary/pq_to_fast5_mapping_rcc12.txt"
    # --- End of Configuration ---
    #
    # # 1. Map all fast5 files to their read_ids
    # fast5_id_to_path_map = find_fast5_files_and_read_ids(INP_DIR)
    #
    # # 2. Map pq files to their corresponding fast5 files
    # pq_to_fast5_map = find_pq_to_fast5_mapping(OUTPQ_DIR, fast5_id_to_path_map)
    #
    # # 3. Write the result to a file
    # write_mapping_to_file(pq_to_fast5_map, OUTPUT_FILE)




