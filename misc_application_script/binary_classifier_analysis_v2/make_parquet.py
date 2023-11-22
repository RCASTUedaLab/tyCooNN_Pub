import sys,os
import glob
import pathlib

sys.path.append("/home/bhaskar/work/tyCooNN_2206/")
import preprocess.TrimAndNormalize as tn
import utils.tyUtils as ut
import pandas as pd
from multiprocessing import Pool
import itertools

trna = sys.argv[1]

paramPath = '/home/bhaskar/work/tyCooNN_2206/test/makePQ.yaml'
param = ut.get_parameter(paramPath)

fast5s = glob.glob(trna + "/fast5/*.fast5")
print("Reading from %d files" % len(fast5s))
with Pool(param.ncore) as p:
    reads = p.map(ut.get_fast5_reads_from_file, fast5s)
    reads = list(itertools.chain.from_iterable(reads))
print("Total reads: %d" % len(reads))

split_list = trna + "/" + trna + "_concordant.lab"
split_list2= trna + "/" + trna + "_discordant.lab"
sp_read = {}
sp_read2= {}
split_species = {}
split_species2_0 = {}
with open(split_list) as fs:
    for line in fs:
        line = line.strip()
        rid, sp = line.split(" ")
        sp_read[rid] = sp
        split_species[sp] = []
        split_species2_0[sp] = []

with open(split_list2) as fs:
    for line in fs:
        line = line.strip()
        rid, sp1, sp2 = line.split(" ")
        sp_read2[rid] = [sp1,sp2]

for read in reads:
    rid = read.read_id
    if rid in sp_read:
        name = sp_read[rid]
        split_species[name].append(read)
    elif rid in sp_read2:
        name1,name2 = sp_read2[rid]
        split_species2_0[name1].append(read)

def make_df(reads,name):
    print("Total reads in %s: %d" % (name,len(reads)))
    trimmed_filterFlgged_read = tn.trimAdaptor(reads,param)
    filtered_reads = [read for read in trimmed_filterFlgged_read if read.filterFlg == 0]
    print("Total reads after filtering: %d" % len(filtered_reads))
    format_reads = tn.formatSignal(filtered_reads,param)

    datalist = []
    for read in format_reads:
        tp = (read.read_id,trna,read.formatSignal)
        datalist.append(tp)

    df = pd.DataFrame(datalist, columns=['read_id', 'trna', 'trimsignal'])
    return df

for name in split_species.keys():
    s_reads = split_species[name]

    outpq = trna + "/" + name + "_con12.pq"
    df = make_df(s_reads,name)
    df.to_parquet(outpq)

    s_reads = split_species2_0[name]
    outpq = trna + "/" + name + "_dis01.pq"
    df = make_df(s_reads,name)
    df.to_parquet(outpq)


