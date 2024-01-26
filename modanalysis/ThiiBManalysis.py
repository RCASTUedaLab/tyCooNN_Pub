


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

def countmismatchratio(bam,position,ref,idd):

    samfile1 = pysam.AlignmentFile(bam, "rb")
    for record in SeqIO.parse(ref, "fasta"):

        id = record.id
        if id == idd:
            seq = record.seq.replace("U", "T")

    id = idd
    a = getratio(samfile1, id, seq)
    fpadapterlen = 4
    a = a[fpadapterlen:]

    r = a[position]
    return r

import random
def getR(mr0,mr1,r,nmax):

    cntref = 0
    cntdiff = 0
    for n in range(nmax):

        random_number = random.random()
        takefromko = random_number < r
        if takefromko:

            nm = mr0[0]+mr0[1]
            random_number = random.randint(1, nm)
            if random_number < mr0[0]:
                cntref+=1
            else:
                cntdiff += 1
        else:
            nm = mr1[0] + mr1[1]
            random_number = random.randint(1, nm)
            if random_number < mr1[0]:
                cntref += 1
            else:
                cntdiff += 1

    return (cntdiff/nmax) * 100

def printR(mr0,mr1):

    rr = [0,0.25,0.5,0.75,1]
    for r in rr:

        ratio = getR(mr0,mr1,r,300)
        print(ratio)

ref = "/home/ueda/project/tRex/referencetest/ecolitRNA_unmod_full.fa"
# bam0 = "/mnt/share/ueda/TyCooNNPub/wt/sample_sorted.bam"
bam0 = "/mnt/share/ueda/trna_data/bam/rcc_sorted.bam"
bam1 = "/mnt/share/ueda/TyCooNNPub/thii/sample_sorted.bam"

position = 7
mr0 = countmismatchratio(bam0, position, ref,"Asn")
print(mr0)

position = 7
mr1 = countmismatchratio(bam1, position, ref,"Asn")
print(mr1)

printR(mr0,mr1)

position = 7
mr0 = countmismatchratio(bam0, position, ref,"Pro2")
print(mr0)

position = 7
mr1 = countmismatchratio(bam1, position, ref,"Pro2")
print(mr1)

printR(mr0,mr1)

position = 7
mr0 = countmismatchratio(bam0, position, ref,"Pro3")
print(mr0)

position = 7
mr1 = countmismatchratio(bam1, position, ref,"Pro3")
print(mr1)

printR(mr0,mr1)

position = 7
mr0 = countmismatchratio(bam0, position, ref,"Tyr1")
print(mr0)

position = 7
mr1 = countmismatchratio(bam1, position, ref,"Tyr2")
print(mr1)

printR(mr0,mr1)