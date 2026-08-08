#!/usr/bin/env python3
import os
import glob
import pandas as pd


def summarize_flags(base_dir: str, output_csv="/home/ueda/flag_summary.csv"):
    pattern = os.path.join(base_dir, "*", "flag.csv")
    files = glob.glob(pattern)

    results = []

    for f in files:
        folder = os.path.basename(os.path.dirname(f))

        # flag.csv
        df = pd.read_csv(
            f, delim_whitespace=True, header=None,
            names=["idx1", "idx2", "flag", "count"]
        )

        # flag Wv
        agg = df.groupby("flag")["count"].sum().to_dict()

        row = {"Folder": folder}

        # ---- { ----
        row["Total"] = sum(agg.values())
        row["Pass"] = agg.get("pass", 0)
        row["Mean Phred Quality"] = agg.get("meanqlow", 0)
        row["Max. Signal Length"] = agg.get("siglen", 0)

        # ---- v] ----
        row["Min. Duration Rate low"] = agg.get("readlenlow", 0)
        row["Min. Duration Rate high"] = agg.get("readlenhigh", 0)

        row["Delta Normalised Signal low"] = agg.get("deltalow", 0)
        row["Delta Normalised Signal high"] = agg.get("delthigh", 0)

        # Sequence Length ` readlenhigh `
        row["Sequence Length"] = agg.get("readlenhigh", 0)

        results.append(row)

    out = pd.DataFrame(results)

    # w
    ordered = [
        "Folder",
        "Total",
        "Pass",
        "Mean Phred Quality",
        "Max. Signal Length",
        "Min. Duration Rate low",
        "Min. Duration Rate high",
        "Delta Normalised Signal low",
        "Delta Normalised Signal high",
        "Sequence Length",
    ]

    out = out[ordered]
    out.to_csv(os.path.join(base_dir, output_csv), index=False)

    print(f"\n? Saved  {os.path.join(base_dir, output_csv)}")
    print(out.head())


if __name__ == "__main__":
    summarize_flags("/mnt/share/bhaskar/rcc_data_extend/infer_50")
