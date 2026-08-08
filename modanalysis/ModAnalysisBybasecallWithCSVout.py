import os
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.transforms as mtransforms
import matplotlib.patches as patches
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from Bio import SeqIO
import pysam
from scipy.stats import chi2_contingency
from sklearn.metrics import roc_curve, auc


# --- Helper Functions ---

def getratio(samfile, id, seq):
    """Calculate match/mismatch counts for each position in the sequence."""
    ret = np.zeros((len(seq), 2))
    end = len(seq)
    for pileupcolumn in samfile.pileup(id, 1, end, min_mapping_quality=0, min_base_quality=0, stepper='nofilter'):
        ref_base = seq[pileupcolumn.pos]
        depth = pileupcolumn.n
        matchcnt = 0
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                query_position = pileupread.query_position
                base = pileupread.alignment.query_sequence[query_position]
                if ref_base == base:
                    matchcnt += 1
            else:
                depth -= 1

        ret[pileupcolumn.pos][0] = matchcnt
        ret[pileupcolumn.pos][1] = depth - matchcnt
    return ret


def getDiffLog(a, b):
    """Calculate log2 ratio of mismatch rates between two samples."""
    ret = []
    for n in range(len(a)):
        if n >= len(b): continue
        aa, bb = a[n], b[n]
        t1, t2 = aa[0] + aa[1], bb[0] + bb[1]
        lg = 0
        if t1 > 10 and t2 > 10:
            r1, r2 = aa[1] / t1, bb[1] / t2
            if r2 > 0 and max(r1, r2) > 0.2:
                diffr = r1 / r2
                if diffr > 0:
                    lg = np.log2(diffr)
        ret.append(lg)
    return ret


def getRangeMap(matrixref):
    """Load reference map from CSV."""
    result_dict = {}
    with open(matrixref, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            result_dict[row[0]] = row[1:]
    return result_dict


def getModIndex(rangemap, ids):
    """Identify modified positions and their types from the range map."""
    result_list, modlist = [], []
    for rowidx, id in enumerate(ids):
        rows = rangemap[id]
        for colidx, row in enumerate(rows):
            if str(row) not in ["A", "U", "G", "C", "nan", "NAN", ""] and len(row) > 0:
                result_list.append((rowidx, colidx))
                modlist.append(str(row))
    return result_list, modlist


def insertNaN(lst, row):
    """Insert NaN values into the list based on empty positions in the reference row."""
    results = []
    count, non_empty_count = 0, 0
    for item in row:
        if item == "":
            count += 1
            if count == 1: start_index = non_empty_count
        else:
            non_empty_count += 1
            if count > 0:
                results.append((start_index, count))
                count = 0
    if count > 0: results.append((start_index, count))

    for index, count in sorted(results, reverse=True):
        lst[index:index] = [np.nan] * count
    return lst


# --- Analysis & Export Functions ---

def analyseROC(matrix, modindexs, modlist, outdir):
    """Perform ROC analysis and export ROC curve points (FPR, TPR, threshold)."""

    data, ans = [], []

    # Collect all valid points
    for row in range(len(matrix)):
        for col in range(len(matrix[0])):
            val = matrix[row][col]
            if np.isnan(val):
                continue
            data.append(val)
            ans.append(1 if (row, col) in modindexs else 0)

    data = np.array(data)
    ans = np.array(ans)

    # --- ROC calculation ---
    fpr, tpr, thresholds = roc_curve(ans, data)
    auc_val = auc(fpr, tpr)

    # --- Export ROC curve points ---
    roc_curve_df = pd.DataFrame({
        "FPR": fpr,
        "TPR": tpr,
        "Threshold": thresholds
    })
    roc_curve_df.to_csv(
        os.path.join(outdir, "data_roc_curve_points.csv"),
        index=False
    )

    # --- Plot ---
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, label=f'AUC = {auc_val:.3f}')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    plt.savefig(os.path.join(outdir, "roc.png"), dpi=300)
    plt.close()


def analyseMod(matrix, modindexs, modlist, outdir):
    """Analyze modification statistics and export to CSV with reference names."""

    # Names corresponding to row indices (from the image)
    # Ensure the order matches the matrix row order
    row_names = [
        "Ala1B", "Ala2", "Arg2", "Arg3", "Arg4", "Arg5", "Asn", "Asp", "Cys",
        "fMet1", "fMet2", "Gln1", "Gln2", "Glu", "Gly1", "Gly2", "Gly3", "His",
        "Ile1", "Ile2", "Ile2v", "Leu1", "Leu1_P", "Leu2", "Leu3", "Leu4", "Leu5",
        "Lys", "Met", "Phe", "Pro1", "Pro2", "Pro3", "Sec", "Ser1", "Ser2", "Ser3",
        "Ser5", "Thr1", "Thr2", "Thr3", "Thr4", "Trp", "Tyr1", "Tyr2", "Val1", "Val2A", "Val2B"
    ]

    d = {}
    modcols = set()
    all_stat_rows = []

    # 1. Process Modified Sites
    for n, (row, col) in enumerate(modindexs):
        mod_type = modlist[n]
        val = matrix[row][col]
        res = analyzeResolution(matrix, modindexs, row, col)
        modcols.add(col)

        if mod_type not in d: d[mod_type] = []
        d[mod_type].append((val, res))

        # Get name from list (fallback to index if out of range)
        name = row_names[row] if row < len(row_names) else f"row_{row}"

        all_stat_rows.append({
            'Type': mod_type,
            'Name': name,  # Converted from row index
            'Col': col,
            'Val': val,
            'Res': res,
            'Is_Mod': 1
        })

    # 2. Process Unmodified Sites
    unmodlist = []
    for n in range(min(76, matrix.shape[1])):
        if n not in modcols:
            for m in range(len(matrix)):
                val = matrix[m][n]
                if not np.isnan(val) and val >= 0:
                    res = analyzeResolution(matrix, modindexs, m, n)
                    unmodlist.append((val, res))

                    name = row_names[m] if m < len(row_names) else f"row_{m}"

                    all_stat_rows.append({
                        'Type': 'unmod',
                        'Name': name,  # Converted from row index (m)
                        'Col': n,
                        'Val': val,
                        'Res': res,
                        'Is_Mod': 0
                    })

    # Export to CSV
    stats_df = pd.DataFrame(all_stat_rows)
    stats_df.to_csv(os.path.join(outdir, "data_mod_stats.csv"), index=False)

    # Plot results
    if unmodlist:
        avenonmod = np.mean([x[0] for x in unmodlist])
        nonmodsd = np.std([x[0] for x in unmodlist], ddof=1)
        plotScater(d, outdir, avenonmod, nonmodsd)


def analyzeResolution(matrix, modindexs, row, col):
    """Determine the base resolution (spread of signal)."""
    val = matrix[row][col]
    if val < 0.2: return 1
    thres = val / 2
    left, right = 0, 0
    # Left spread
    for m in range(col - 1, 0, -1):
        v = matrix[row][m]
        if v < thres or np.isnan(v): break
        if (row, m) not in modindexs: left += 1
    # Right spread
    for m in range(col + 1, min(76, matrix.shape[1])):
        v = matrix[row][m]
        if v < thres or np.isnan(v): break
        if (row, m) not in modindexs: right += 1
    return left + right + 1


def plotScater(d, outdir, avenonmod, nonmodsd):
    """Save the comparison plot of modifications."""
    # (Plotting logic simplified for execution)
    plt.figure()
    for mod, values in d.items():
        y, x = zip(*values)
        plt.scatter(x, y, label=mod)
    plt.axhline(y=avenonmod, color='r', linestyle='--')
    plt.legend()
    plt.savefig(os.path.join(outdir, "all_mods_scatter.png"))


# --- Main Entry Point ---

def compare(bam1, bam2, ref, matrixref, outdir, ez):
    if not os.path.exists(outdir): os.makedirs(outdir)

    samfile1 = pysam.AlignmentFile(bam1, "rb")
    samfile2 = pysam.AlignmentFile(bam2, "rb")

    records = list(SeqIO.parse(ref, "fasta"))
    records = sorted(records, key=lambda x: x.id.lower())
    ids = [r.id for r in records]

    rangemap = getRangeMap(matrixref)
    modindexs, modlist = getModIndex(rangemap, ids)

    matrix = []
    for record in records:
        
        id, seq = record.id, record.seq.replace("U", "T")
        a = getratio(samfile1, id, seq)
        b = getratio(samfile2, id, seq)
        log_ratio = getDiffLog(a, b)

        # Clipping adapters
        fpadapterlen = 3 if id == "His" else 4
        log_ratio = log_ratio[fpadapterlen: len(log_ratio) - 50]
        log_ratio = insertNaN(log_ratio, rangemap[id])
        matrix.append(log_ratio)

    # Align matrix for CSV and Heatmap
    max_len = max(len(row) for row in matrix)
    padded_matrix = [row + [np.nan] * (max_len - len(row)) for row in matrix]
    np_matrix = np.array(padded_matrix)

    # 1. Export Heatmap Matrix to CSV
    df_matrix = pd.DataFrame(np_matrix, index=ids)
    df_matrix.to_csv(os.path.join(outdir, ez+"_data_heatmap_matrix.csv"))

    # 2. Generate Heatmap Plot
    plt.figure(figsize=(20, 10))
    sns.heatmap(np_matrix, mask=np.isnan(np_matrix), cmap='bwr', center=0)
    plt.yticks(ticks=np.arange(len(ids)) + 0.5, labels=ids, rotation=0)
    plt.savefig(os.path.join(outdir, "heatmap.png"), dpi=300)

    # 3. Analyze and export data points
    analyseMod(np_matrix, modindexs, modlist, outdir)
    analyseROC(np_matrix, modindexs, modlist, outdir)

    samfile1.close()
    samfile2.close()

bam2 = "/mnt/share/ueda/trna_data/bam/rcc_sorted.bam"
bam1 = "/mnt/share/ueda/trna_data/bam/ivt_sorted.bam"
ref =  "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"
# matrixref = "/mnt/share/ueda/trna_data/tRNAEcoli.csv"
# matrixref = "/mnt/share/ueda/trna_data/ecoliTrna_paper.csv"
matrixref = "/mnt/share/ueda/trna_data/ecoliTrna_paper_25.csv"
outdir = "/home/ueda/ivt_2025"

compare(bam2,bam1,ref,matrixref,outdir,"rcc_ivt")


ref =  "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"


enzimes = ["thii_wt","trmb_wt","trub_wt"]

for ez in enzimes:
    # bam1 = "/mnt/share/ueda/TyCooNNPub/"+ez+"/sample_sorted.bam"
    bam2 = "/mnt/share/ueda/TyCooNNPub/wt/sample_sorted.bam"
    bam2 = "/mnt/share/ueda/trna_data/bam/rcc_sorted.bam"
    outdir = "/home/ueda/ivt_2025/"+ez
    compare(bam2,bam1,ref,matrixref,outdir,ez)
