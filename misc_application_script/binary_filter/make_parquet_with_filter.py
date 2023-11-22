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

paramPath = 'makePQ.yaml'
param = ut.get_parameter(paramPath)

fast5s = glob.glob(trna + "/fast5/*.fast5")
print("Reading from %d files" % len(fast5s))
with Pool(param.ncore) as p:
    reads = p.map(ut.get_fast5_reads_from_file, fast5s)
    reads = list(itertools.chain.from_iterable(reads))
print("Total reads: %d" % len(reads))

filter_list = trna + "/" + trna + "_filter.label"
sp_read = {}
split_species = {}
with open(filter_list) as fs:
    for line in fs:
        line = line.strip()
        rid, sp, tt = line.split(" ")
        sp_read[rid] = sp
        split_species[sp] = []

for read in reads:
    rid = read.read_id
    if rid in sp_read:
        name = sp_read[rid]
        split_species[name].append(read)

for name in split_species.keys():
    s_reads = split_species[name]
    print("Total reads in %s: %d" % (name,len(s_reads)))
    trimmed_filterFlgged_read = tn.trimAdaptor(s_reads,param)
    filtered_reads = [read for read in trimmed_filterFlgged_read if read.filterFlg == 0]
    print("Total reads after filtering: %d" % len(filtered_reads))
    format_reads = tn.formatSignal(filtered_reads,param)

    datalist = []
    for read in format_reads:
        tp = (read.read_id,trna,read.formatSignal)
        datalist.append(tp)
    
    outpq = trna + "/" + name + ".pq"

    df = pd.DataFrame(datalist, columns=['read_id', 'trna', 'trimsignal'])
    df.to_parquet(outpq)

