import os
from Bio import SeqIO
import pysam
import numpy as np

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
        if (t1 > 100 and t2 > 100):
            r1 = aa[1] / t1
            r2 = bb[1] /t2

            if r2 > 0 and max(r1,r2)> 0.2:
                diffr = r1/r2
                if diffr > 0:
                    lg = np.log2(diffr)


        ret.append(lg)
    return ret

import matplotlib.pyplot as plt

def compare(bam1,bam2,ref):


    samfile1 = pysam.AlignmentFile(bam1, "rb")
    samfile2 = pysam.AlignmentFile(bam2, "rb")

    lst = []
    for record in SeqIO.parse(ref, "fasta"):

        id = record.id
        if id == "Gly3":
            seq = record.seq.replace("U", "T")

    id = "Gly3"
    a = getratio(samfile1, id, seq)
    b = getratio(samfile2, id, seq)
    log_ratio = getDiffLog(a, b)

    fpadapterlen = 4
    log_ratio = log_ratio[fpadapterlen:]
    tpadapterlen = 50
    log_ratio = log_ratio[0:len(log_ratio) - tpadapterlen]
    #
    samfile1.close()
    samfile2.close()

    return log_ratio

ref = "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"

bam0 = "/mnt/share/ueda/minimapmapping/gly3_truB/ko/gly3_ko.bam"

bam1 = "/mnt/share/ueda/minimapmapping/gly3_truB/wt/gly3_wt.bam"
bam2 = "/mnt/share/ueda/minimapmapping/gly3_truB/wt3ko1/gly3_wt3ko1.bam"
bam3 = "/mnt/share/ueda/minimapmapping/gly3_truB/wt1ko1/gly3_wt1ko1.bam"
bam4 = "/mnt/share/ueda/minimapmapping/gly3_truB/wt1ko3/gly3_wt1ko3.bam"
bams = [bam1,bam2,bam3,bam4]
labels = ["WT","25%","50%","75%"]

savef = "/mnt/share/ueda/minimapmapping/gly3_truB.png"

data = []
n = 0
for bam in bams:

    log_ratio =  compare(bam, bam0, ref)
    x = np.arange(start=0, stop= len(log_ratio), step=1)
    lb = labels[n]
    plt.plot(x, log_ratio, label=lb)
    n +=1

plt.legend()

plt.savefig(savef,dpi=450)

print(data)

def compare(bam1,bam2,ref):


    samfile1 = pysam.AlignmentFile(bam1, "rb")
    samfile2 = pysam.AlignmentFile(bam2, "rb")

    lst = []
    for record in SeqIO.parse(ref, "fasta"):

        id = record.id
        if id == "Gly3":
            seq = record.seq.replace("U", "T")

    id = "Gly3"
    a = getratio(samfile1, id, seq)
    b = getratio(samfile2, id, seq)
    log_ratio = getDiffLog(a, b)

    fpadapterlen = 4
    log_ratio = log_ratio[fpadapterlen:]
    tpadapterlen = 50
    log_ratio = log_ratio[0:len(log_ratio) - tpadapterlen]
    #
    samfile1.close()
    samfile2.close()

    return log_ratio


def countmismatchratio(bam,position,ref):

    samfile1 = pysam.AlignmentFile(bam, "rb")
    for record in SeqIO.parse(ref, "fasta"):

        id = record.id
        if id == "Gly3":
            seq = record.seq.replace("U", "T")

    id = "Gly3"
    a = getratio(samfile1, id, seq)
    fpadapterlen = 4
    a = a[fpadapterlen:]

    r = a[position]
    return r

bamall = [bam0,bam1,bam2,bam3,bam4]
position = 54
for bam in bamall:

    mr = countmismatchratio(bam,position,ref)
    print(mr)
    r = mr[1]/ sum(mr)
    print(r)