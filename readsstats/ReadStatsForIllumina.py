import os
import re
import pandas as pd

base_dir = "/mnt/share/bhaskar/fromhome/illumina/short_read/result"

records = []
# t@CAIT
for root, _, files in os.walk(base_dir):
    for fn in files:
        if not fn.endswith(".txt"):
            continue
        path = os.path.join(root, fn)
        before = None
        after  = None
        with open(path, encoding="utf-8") as f:
            lines = f.readlines()
        # es
        for i, line in enumerate(lines):
            l = line.strip().lower()
            # Read1 before filtering As~ total reads T
            if "read1 before filtering" in l or "Read1 before filtering" in l:
                for nl in lines[i+1:]:
                    m = re.search(r"total reads[:F]\s*([0-9]+)", nl, re.IGNORECASE)
                    if m:
                        before = m.group(1)
                        print(fn,m.group(1))
                        break
            # Read1 after filtering AlT
            if "read1 after filtering" in l or "Read1 after filtering" in l:
                for nl in lines[i+1:]:
                    m = re.search(r"total reads[:F]\s*([0-9]+)", nl, re.IGNORECASE)
                    if m:
                        after = m.group(1)
                        print(fn,m.group(1))
                        break
            # o[v
            if before is not None and after is not None:
                break

        records.append({
            "file": fn,
            "before_reads": before,
            "after_reads":  after
        })


# DataFrame \
df = pd.DataFrame(records)
print(df)
# Kv CSV o
df.to_csv("read_counts_summary.csv", index=False)
