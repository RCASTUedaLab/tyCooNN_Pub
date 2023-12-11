from Bio import SeqIO
import pysam

# bam1 = "/mnt/share/ueda/minimapmapping/gly3_truB/ko/gly3_ko.bam"
# bam2 = "/mnt/share/ueda/minimapmapping/gly3_truB/wt/gly3_wt.bam"
# ref =  "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"
# outdir = "/mnt/share/ueda/minimapmapping/"

bam2 = "/mnt/share/ueda/trna_data/bam/rcc_sorted.bam"
bam1 = "/mnt/share/ueda/trna_data/bam/ivt_sorted.bam"
ref =  "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"
matrixref = "/mnt/share/ueda/trna_data/tRNAEcoli.csv"
outdir = "/mnt/share/ueda/minimapmapping/"

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def getratio(samfile, id, seq):

    ret = np.zeros((len(seq),2))

    end = len(seq)
    for pileupcolumn in samfile.pileup(id, 1, end, min_mapping_quality=0,min_base_quality=0,stepper='nofilter'):

        ref = seq[pileupcolumn.pos]
        depth = pileupcolumn.n
        matchcnt = 0
        miscount = 0
        a,t,c,g = 0,0,0,0
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:

                query_position = pileupread.query_position
                base = pileupread.alignment.query_sequence[query_position]
                if (ref == base):
                    matchcnt+=1
                if base =="A":
                    a+=1
                if base =="T":
                    t+=1
                if base =="C":
                    c+=1
                if base =="G":
                    g+=1

            else:
                depth-=1


        miscount = depth - matchcnt
        # print(pileupcolumn.pos, ref, matchcnt, miscount, a, t, c, g)
        ret[pileupcolumn.pos][0] = matchcnt
        ret[pileupcolumn.pos][1] = miscount

    return ret

import scipy.stats as stats
def getDiffLog(a,b):

    ret = []
    for n in range(len(a)):

        if n >= len(b):
            continue

        aa = a[n]
        bb = b[n]

        table = [[aa[0],aa[1]],  # group1f
                 [bb[0],bb[1]]]  # group2f
        t1 = aa[0]+aa[1]
        t2 = bb[0] + bb[1]
        diffr = 0
        lg = 0
        if (t1 > 10 and t2 > 10):
            r1 = aa[1] / t1
            r2 = bb[1] /t2

            if r2 > 0 and max(r1,r2)> 0.2:
                diffr = r1/r2
                if diffr > 0:
                    lg = np.log2(diffr)


        ret.append(lg)
    return ret

# def getRange(id):
#
#     lst = []
#     # 1st
#     if id !="His":
#         lst.append((0,1))
#
#     #  2nd
#     if id in ["Arg5","Cys","Gln1","Gln2","Glu","Gly1","Gly2","Ile2","Ile2v","Sec","Ser1","Ser2","Ser3","Ser5"]:
#         lst.append((16, 2))
#     elif id in ["fMet1","fMet2","Pro1","Pro2","Pro3"]
#         pass
#     else:
#         lst.append((17,1))
#
#
#
#         lst.append((20, 2))
#         lst.append((45, 19))
#     return lst

def getRange(row):

    results = []
    count = 0
    non_empty_count = 0

    for i, item in enumerate(row):
        if item == "":
            count += 1
            if count == 1:  # Anvf
                start_index = non_empty_count
        else:
            non_empty_count += 1
            if count > 0:  # AI
                results.append((start_index, count))
                count = 0

    # vf
    if count > 0:
        results.append((start_index, count))

    return results


def insertNaN(lst,row):

    index_counts = getRange(row)
    for index, count in sorted(index_counts, reverse=True):
        lst[index:index] = [np.nan] * count
    return lst

import csv
def getRangeMap(matrixref):

    result_dict = {}
    with open(matrixref, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            key = row[0]
            result_dict[key] = row[1:]
    return result_dict

def getModIndex(rangemap,ids):

    result_list = []
    modlist = []
    rowidx = 0
    for id in ids:
        rows = rangemap[id]
        colidx = 0
        for row in rows:

            if str(row) in ["A","U","G","C"]:
                pass
            elif str(row) not in ["nan","NAN"]:
                if len(row) > 0:
                    result_list.append((rowidx,colidx))
                    modlist.append(str(row))
            colidx+=1
        rowidx+=1

    return result_list,modlist


from matplotlib.colors import TwoSlopeNorm
import pandas as pd
import matplotlib.transforms as mtransforms
import matplotlib.patches as patches
def compare(bam1,bam2,ref,matrixref,outdir):

    samfile1 = pysam.AlignmentFile(bam1, "rb")
    samfile2 = pysam.AlignmentFile(bam2, "rb")

    matrix = []
    ids = []
    lst = []
    for record in SeqIO.parse(ref, "fasta"):

        id = record.id
        seq = record.seq.replace("U", "T")
        lst.append((id,seq))


    lst = sorted(lst)
    lst = sorted(lst, key=lambda x: x[0].lower())

    ids = list(map(lambda x: x[0],lst))
    rangemap = getRangeMap(matrixref)
    modindexs,modlist = getModIndex(rangemap,ids)

    for id,seq in lst:

        # id = record.id
        # ids.append(id)
        # seq = record.seq.replace("U","T")
        #
        print(id,len(seq))
        a = getratio(samfile1, id, seq)
        b = getratio(samfile2, id, seq)
        log_ratio = getDiffLog(a,b)

        fpadapterlen = 4
        if id == "His":
            fpadapterlen = 3
        print(id, len(log_ratio))
        log_ratio = log_ratio[fpadapterlen:]
        tpadapterlen = 50
        log_ratio = log_ratio[0:len(log_ratio)-tpadapterlen]
        log_ratio = insertNaN(log_ratio,rangemap[id])
        matrix.append(log_ratio)
        print(id,len(log_ratio))


    out = outdir+"/rcc_ivt.png"

    plt.figure(figsize=(30, 10),dpi=300)
    fig, ax = plt.subplots(figsize=(30, 10))
    max_length = max(len(sublist) for sublist in matrix)
    padded_data = [sublist + [0] * (max_length - len(sublist)) for sublist in matrix]
    matrix = np.array(padded_data)
    mask = np.isnan(matrix)


    non_nan_indices = [i for i in range(matrix.shape[1]) if not np.isnan(matrix[0,i])]
    zero_based_indices = list(range(len(non_nan_indices)))

    sns.heatmap(mask, cmap=['gray', 'gray'], cbar=False, alpha=0.4,linewidths=0,linecolor='none')
    ax = sns.heatmap(matrix, mask=mask, cmap='bwr', center=0, cbar=True,linewidths=0)
    tick = np.array(non_nan_indices[::5])-0.5
    ax.set_xticks(tick)
    ax.set_xticklabels(zero_based_indices[::5])

    # print(ids)
    # print(len(ids))
    # put y label
    plt.yticks(ticks=np.arange(48),labels=ids,fontsize=12,ha='left')
    ax.tick_params(axis='y', which='major', pad=40)
    # lower y axis
    offset = mtransforms.ScaledTranslation(0, -0.05, fig.dpi_scale_trans)
    for label in ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

    highlight_cells = modindexs
    for cell in highlight_cells:
        y, x = cell
        plt.gca().add_patch(patches.Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=1,linestyle='--'))

    plt.savefig(out)

    samfile1.close()
    samfile2.close()

    #analyse modification
    analyseMod(matrix,modindexs,modlist)

    #ROC analysis
    analyseROC(matrix, modindexs, modlist,outdir)

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import auc
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
def analyseROC(matrix, modindexs, modlist,outdir):

    n = 0
    data = []
    ans = []

    ydata = []
    yans = []
    m7Gdata = []
    m7Gans = []
    s4Udata = []
    s4Uans = []

    for row in range(len(matrix)):
        for col in range(len(matrix[0])):

            rc = (row,col)
            val = matrix[row][col]
            if np.isnan(val):
                continue

            data.append(val)
            isMod = rc in modindexs
            if isMod:
                ans.append(1)
                idx = modindexs.index(rc)
                mod_kind = modlist[idx]

                if mod_kind == "Y":
                    ydata.append(val)
                    yans.append(1)
                if mod_kind == "m7G":
                    m7Gdata.append(val)
                    m7Gans.append(1)
                if mod_kind == "s4U":
                    s4Udata.append(val)
                    s4Uans.append(1)

            else:
                ans.append(0)
                yans.append(0)
                m7Gans.append(0)
                s4Uans.append(0)

                ydata.append(val)
                m7Gdata.append(val)
                s4Udata.append(val)


    fpr_all, tpr_all, _ = roc_curve(ans, data)
    fpr_y, tpr_y, _ = roc_curve(yans, ydata)
    fpr_m7G, tpr_m7G, _ = roc_curve(m7Gans, m7Gdata)
    fpr_s4U, tpr_s4U, _ = roc_curve(s4Uans, s4Udata)


    auc_all = auc(fpr_all, tpr_all)
    auc_y = auc(fpr_y, tpr_y)
    auc_m7G = auc(fpr_m7G, tpr_m7G)
    auc_s4U = auc(fpr_s4U, tpr_s4U)

    # ROCvbg
    plt.figure(figsize=(8, 6))

    plt.plot(fpr_all, tpr_all, color='blue', lw=2, label=f'All modification (AUC = {auc_all:.2f})')
    psi = "\u03A8"
    plt.plot(fpr_y, tpr_y, color='green', lw=2, label=psi+f' (AUC = {auc_y:.2f})')
    plt.plot(fpr_m7G, tpr_m7G, color='red', lw=2, label=f'm7G (AUC = {auc_m7G:.2f})')
    plt.plot(fpr_s4U, tpr_s4U, color='red', lw=2, label=f's4U (AUC = {auc_s4U:.2f})')

    plt.plot([0, 1], [0, 1], color='darkgray', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) - Multiple Models')
    plt.legend(loc="lower right")


    out = outdir + "/roc.png"
    plt.savefig(out)


def analyzeResolution(matrix,modindexs,row,col):

    thres = 0.2
    val = matrix[row][col]
    if val < thres:
        return 1

    m = col-1
    left = 0
    while m > 0:
        val = matrix[row][m]

        if val < thres:
            break
        else:
            if not np.isnan(val):
                if (row,m) not in modindexs:
                    left+=1
        m -= 1

    m = col+1
    right = 0
    while m < 75:
        val = matrix[row][m]

        if val < thres:
            break
        else:
            if not np.isnan(val):
                if (row, m) not in modindexs:
                    right += 1
        m += 1

    return right+left+1

def analyseMod(matrix,modindexs,modlist):

    d = {}
    n = 0
    for mod in modlist:

      row,col = modindexs[n]
      n+=1
      # print(mod,matrix[row][col],row,col)
      val = matrix[row][col]
      resolution = analyzeResolution(matrix,modindexs,row,col)
      if mod in d:
          lst = d[mod]
      else:
          lst = []
          d[mod] = lst
      lst.append((val, resolution))

    for k in d:
        pk = k
        if k == "Y":
            pk = chr(0x03A8)
        print(pk,d[k])

    plotScater(d)

def replace_negatives_with_zero(lst):
    return [0 if x < 0 else x for x in lst]

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import auc
from matplotlib.gridspec import GridSpec
def plotScater(d):

    #psude
    psude = d["Y"]
    m7G = d["m7G"]
    s4U = d["s4U"]

    y1, x1 = zip(*psude)
    y2, x2 = zip(*m7G)
    y3, x3 = zip(*s4U)

    y1 = replace_negatives_with_zero(y1)
    y2 = replace_negatives_with_zero(y2)
    y3 = replace_negatives_with_zero(y3)


    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4)
    ax_scatter = fig.add_subplot(gs[1:4, 0:3])
    ax_histx = fig.add_subplot(gs[0, 0:3], sharex=ax_scatter)
    ax_histy = fig.add_subplot(gs[1:4, 3], sharey=ax_scatter)
    ax_scatter.scatter(x1, y1, color="blue", label=chr(0x03A8)+"")
    ax_scatter.scatter(x2, y2, color="red", label="m7G")
    ax_scatter.scatter(x3, y3, color="green", label="s4U")
    sns.kdeplot(x1, ax=ax_histx, color="blue", fill=True)
    sns.kdeplot(x2, ax=ax_histx, color="red", fill=True)
    sns.kdeplot(x3, ax=ax_histx, color="green", fill=True)
    sns.kdeplot(y1, ax=ax_histy, color="blue", fill=True, vertical=True)
    sns.kdeplot(y2, ax=ax_histy, color="red", fill=True, vertical=True)
    sns.kdeplot(y3, ax=ax_histy, color="green", fill=True, vertical=True)

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_scatter.set_xlabel("Base resolution")
    ax_scatter.set_ylabel("Misbasecall ratio (log2)")
    ax_scatter.legend()
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    plt.savefig("/mnt/share/ueda/minimapmapping/threemoddist.png")

    plt.clf()
    x = []
    y = []
    labels=[]
    for key in d:

        lst = d[key]
        if key == "Y":
            key = chr(0x03A8)
        labels.append(key)
        ye, xe = zip(*lst)
        x.append(np.mean(xe))
        y.append(np.mean(ye))

    plt.scatter(x, y)  # Uz}
    zerotxt = ""
    for i in range(len(x)):

        if y[i] > 0:
            if labels[i] == "acp3U":
                plt.text(x[i] + 0.05, y[i]+0.05, labels[i])
            elif  labels[i] == "Gm":
                plt.text(x[i] + 0.05, y[i]+0.1, labels[i])
            else:
                plt.text(x[i]+0.05, y[i], labels[i])
        else:
            zerotxt = zerotxt +"," + labels[i]

    plt.text(1.05, 0, zerotxt[1:])
    plt.xlabel("Base resolution")
    plt.ylabel("Misbasecall ratio (log2)")

    plt.savefig("/mnt/share/ueda/minimapmapping/all.png")
    plt.clf()

    # plot bar
    categories = []
    trials = []
    successes = []

    for key in d:

        lst = d[key]
        if key == "Y":
            key = chr(0x03A8)
        categories.append(key)
        trials.append(len(lst))
        ye, xe = zip(*lst)
        threshold = 0.5
        count = len([x for x in ye if x >= threshold])
        successes.append(count)

    sorted_indices = sorted(range(len(trials)), key=lambda i: trials[i], reverse=True)
    categories = [categories[i] for i in sorted_indices]
    trials = [trials[i] for i in sorted_indices]
    successes = [successes[i] for i in sorted_indices]

    success_rates = [(s / t * 100) for s, t in zip(successes, trials)]
    plt.bar(categories, trials, label='detected', color='blue', edgecolor='white', linewidth=1.5)
    plt.bar(categories, successes, label='total', color='green', edgecolor='black', linewidth=1.5)

    for i, rate in enumerate(success_rates):
        if int(rate) < 100 and int(rate) > 0:

            plt.text(i + 0.2, trials[i], f'{int(rate)}%', ha='center', va='bottom')
            # plt.text(categories[i], trials[i], f'{int(rate)}%', ha='center', va='bottom')

    plt.xticks(rotation=90)
    plt.ylabel("Detected counts/Number of modified position")
    plt.savefig("/mnt/share/ueda/minimapmapping/bar.png")

    #ROC analysis
    all = []




    # plt.legend()


compare(bam2,bam1,ref,matrixref,outdir)