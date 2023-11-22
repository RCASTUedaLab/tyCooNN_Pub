from Bio import SeqIO
from Bio import pairwise2
import numpy as np
from multiprocessing import Pool
from functools import partial
from tqdm import *
import sys

trna = sys.argv[1]

default_confusion = 0.1

try:
    c = float(sys.argv[2])
    default_confusion = c
except:
    pass

print("confusion in score: ",default_confusion)

#label_file = trna + "/" + trna + ".label"
#filtered_reads = {}
#total = 0
#with open(label_file) as handle:
#    for line in handle:
#        line = line.strip()
#        values = line.split()
#        map_species = values[4]
#        if map_species != 'None':
#            read_id = values[0]
#            filtered_reads[read_id] = 0
#        total += 1
#
#n = len(list(filtered_reads.keys()))
#print("No. of filtered reads: ",n," out of ",total)

fastq_file = trna + "/fastq/" + trna + ".fastq"
records = list(SeqIO.parse(fastq_file,'fastq'))
seq_data = {}
for record in records:
    rid = record.id
    #if rid in filtered_reads:
    #    seq = str(record.seq)
    #    seq_data[rid] = seq
    seq = str(record.seq)
    seq_data[rid] = seq
n = len(list(seq_data.keys()))
print("No. of fastq reads: ",n)

ref_file = "ref/ref_" + trna + ".fa"
ref_rec = list(SeqIO.parse(ref_file,'fasta'))
reference = {}
for ref in ref_rec:
    reference[ref.id] = str(ref.seq)

def break_confusion_fmet(name_tup,confused):
    name_tup_result = []
    if trna == 'fMet':
        name1 = name_tup[0][0]
        name2 = name_tup[1][0]
        score1 = name_tup[0][1]
        score2 = name_tup[1][1]
        seq1 = (name_tup[0][2],name_tup[0][3])
        seq2 = (name_tup[1][2],name_tup[1][3])
        boundary1 = (name_tup[0][4],name_tup[0][5])
        boundary2 = (name_tup[1][4],name_tup[1][5])
        if name1 == 'fMet1':
            name_tup_result = name_tup
        elif name1 == 'fMet2':
            name_tup_result.append((name2,score2,seq2[0],seq2[1],boundary2[0],boundary2[1]))
            name_tup_result.append((name1,score1,seq1[0],seq1[1],boundary1[0],boundary1[1]))
        return name_tup_result,False
    else:
        return name_tup,confused

def apply_pairwise2(R,Q,match=2,mismatch=-1,gapeopen=-2.7,gapex=-0.4):
    return pairwise2.align.localms(Q,R,match,mismatch,gapeopen,gapex)

def get_best_scoring_alignment(aln,confusion=default_confusion):
    name_tup = []
    for name in aln.keys():
        name_tup.append((name,aln[name].score,aln[name].seqA,aln[name].seqB,aln[name].start,aln[name].end))
    name_tup = sorted(name_tup,key=lambda x:x[1],reverse=True)
    s1 = name_tup[0][1]
    s2 = name_tup[1][1]
    delta_score = s1 - s2
    confused = False
    if delta_score <= confusion:
        confused = True
    if confused:
        name_tup,confused = break_confusion_fmet(name_tup,confused)
    #print(s1,s2,delta_score,confused)
    best = {}
    best['id'] = name_tup[0][0]
    best['score'] = name_tup[0][1]
    best['seq']  = (name_tup[0][2],name_tup[0][3])
    best['boundary'] = (name_tup[0][4],name_tup[0][5])
    best['is_confused'] = confused
    return best

def get_aln_mp(inp,rf):
    rid = inp[0]
    q   = inp[1]
    Aln = {}
    for name in rf.keys():
        r = rf[name]
        alignment = apply_pairwise2(r,q)
        Aln[name] = alignment[0]
    best = get_best_scoring_alignment(Aln)
    if not best['is_confused']:
        return (rid,best['id'])
    else:
        return (rid,None)

pack_input = []
for rid in seq_data.keys():
    Q = seq_data[rid]
    pack_input.append((rid,Q))
with Pool(16) as pool:
    results = []
    with tqdm(total=n) as pbar:
        for r in pool.imap_unordered(partial(get_aln_mp,rf=reference),pack_input):
            results.append(r)
            pbar.update()
    #results = pool.map(partial(inp,rf=reference),pack_input)
#results = []
#for p in pack_input:
#    r = get_aln_mp(p,reference)
#    results.append(r)

read_split = {}
for result in results:
    if result[1] is not None:
        read_split[result[0]] = result[1]

# --------------------------------------------------------
#read_split = {}
#for k,rid in enumerate(seq_data.keys()):
#    Q = seq_data[rid]
#    Aln = {}
#    for name in reference.keys():
#        R = reference[name]
#        alignment = apply_pairwise2(Q,R)
#        Aln[name] = alignment[0]
#    best = get_best_scoring_alignment(Aln)
#    if not best['is_confused']:
#        read_split[rid] = best['id']
#    if k % 100 == 0:
#        print("\rDone %d out of %d" % (k,n),end='')
#print("\rDone %d out of %d" % (k,n))
# --------------------------------------------------------

out_lab = trna + "/"  + trna + "_map.lab"
n1 = len(list(read_split.keys()))
print("n = ",n1)
with open(out_lab,"w") as ohandle:
    for rid in read_split.keys():
        name = read_split[rid]
        ohandle.write("%s %s\n" % (rid,name))

