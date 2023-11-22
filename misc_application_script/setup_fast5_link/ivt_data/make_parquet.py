import sys,os
import glob

sys.path.append("/home/bhaskar/work/tyCooNN_2206/")
import preprocess.TrimAndNormalize as tn
import utils.tyUtils as ut
import pandas as pd
from multiprocessing import Pool
import itertools

trna = sys.argv[1]
outpq = trna + "/" + trna + ".pq" 

paramPath = '/home/bhaskar/work/tyCooNN_2206/test/makePQ.yaml'
param = ut.get_parameter(paramPath)

fast5s = glob.glob(trna + "/*.fast5")
print("Doing %s" % trna)
print("Reading from %d files" % len(fast5s))
with Pool(param.ncore) as p:
    reads = p.map(ut.get_fast5_reads_from_file, fast5s)
    reads = list(itertools.chain.from_iterable(reads))
print("Total reads: %d" % len(reads))
trimmed_filterFlgged_read = tn.trimAdaptor(reads,param)
filtered_reads = [read for read in trimmed_filterFlgged_read if read.filterFlg == 0]
print("Total reads after filtering: %d" % len(filtered_reads))
format_reads = tn.formatSignal(filtered_reads,param)

datalist = []
for read in format_reads:
    tp = (read.read_id,trna,read.formatSignal)
    datalist.append(tp)

df = pd.DataFrame(datalist, columns=['read_id', 'trna', 'trimsignal'])
df.to_parquet(outpq)

