import sys,os
import glob
from multiprocessing import Pool
import itertools
import util as ut

indir = sys.argv[1] + "/fast5/" 
f5list = list(sorted(glob.glob(indir + '*.fast5',recursive=True)))
print("No. of files:",len(f5list))

with Pool(8) as p:
    reads = p.map(ut.get_fast5_reads_from_file, f5list)
    reads = list(itertools.chain.from_iterable(reads))

N = len(reads)
print("No of reads:",N)

outdir = sys.argv[1]
outfile = sys.argv[1] + ".label"
outfile = outdir + "/" + outfile
print("Output at:",outfile)

with open(outfile,"w") as fw:

    for read in reads:
        read_id = read.read_id
        attrs = read.map_attrs
        filterflg = attrs['filterflg']
        filterpass = attrs['filterpass']
        trimSuccess = attrs['trimSuccess']
        tRNA_infer = attrs['tRNA_infer']
        tRNAIndex = attrs['tRNAIndex']
        softmax_prob = attrs['softmax_prob']
        fw.write("%s %d %d %d %6s %d %.4f\n" % (read_id,filterflg,filterpass,trimSuccess,
                                                tRNA_infer,tRNAIndex,softmax_prob))

