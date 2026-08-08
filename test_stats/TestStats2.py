#!/usr/bin/env python3
import os
import glob
import pandas as pd


def summarize_flags(base_dir: str, output_csv="flag_summary2.csv"):
    # �w���������f�B���N�g������ flag.csv ������
    pattern = os.path.join(base_dir, "*", "flag.csv")
    files = glob.glob(pattern)

    results = []
    total_pass_count = 0
    limit = 12000

    for f in files:
        folder = os.path.basename(os.path.dirname(f))

        # flag.csv ����������
        # delim_whitespace=True �� pandas �����V�o�[�W�������� sep='\s+' ��������������
        try:
            df = pd.read_csv(
                f, sep=r'\s+', header=None,
                names=["idx1", "idx2", "flag", "count"]
            )
        except Exception as e:
            print(f"Skipping {f} due to error: {e}")
            continue

        # flag �������J�E���g���W�v
        agg = df.groupby("flag")["count"].sum().to_dict()

        current_pass = agg.get("pass", 0)

        row = {
            "Folder": folder,
            "Total": sum(agg.values()),
            "Pass": current_pass,
            "Mean Phred Quality": agg.get("meanqlow", 0),
            "Max. Signal Length": agg.get("siglen", 0),
            "Min. Duration Rate low": agg.get("readlenlow", 0),
            "Min. Duration Rate high": agg.get("readlenhigh", 0),
            "Delta Normalised Signal low": agg.get("deltalow", 0),
            "Delta Normalised Signal high": agg.get("delthigh", 0),
            "Sequence Length": agg.get("readlenhigh", 0),
        }

        results.append(row)

        # Pass ���������v�Z
        total_pass_count += current_pass

        # 12000�����������I��
        if total_pass_count >= limit:
            print(f"Reached {total_pass_count} 'Pass' counts (limit: {limit}). Stopping early.")
            break

    if not results:
        print("No data found.")
        return

    out = pd.DataFrame(results)

    # �J�����������w��
    ordered = [
        "Folder", "Total", "Pass", "Mean Phred Quality", "Max. Signal Length",
        "Min. Duration Rate low", "Min. Duration Rate high",
        "Delta Normalised Signal low", "Delta Normalised Signal high",
        "Sequence Length",
    ]

    out = out[ordered]

    # �o���p�X�������ibase_dir������ flag_summary2.csv�j
    save_path = os.path.join(base_dir, output_csv)
    out.to_csv(save_path, index=False)

    print(f"\nSaved to: {save_path}")
    print(f"Total accumulated Pass: {total_pass_count}")
    print(out.head())


if __name__ == "__main__":
    # �f�B���N�g���� /mnt/share/bhaskar/rcc_data �����X
    summarize_flags("/mnt/share/bhaskar/rcc_data")