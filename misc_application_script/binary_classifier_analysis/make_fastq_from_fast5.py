import sys,os
import glob
from multiprocessing import Pool
import itertools
import util as ut

trna = sys.argv[1]

indir = trna + "/fast5/" 
f5list = list(sorted(glob.glob(indir + '*.fast5',recursive=True)))
print("No .of files:",len(f5list))

rid_filter = {}
with open(trna + "/" + trna + ".label") as fl:
    for line in fl:
        line = line.strip()
        array = line.split(" ")
        rid = array[0]
        map_info = int(array[1])
        if map_info != -1:
            rid_filter[rid] = 0

with Pool(8) as p:
    reads = p.map(ut.get_fast5_reads_from_file, f5list)
    reads = list(itertools.chain.from_iterable(reads))

N = len(reads)
print("No of reads:",N)

outdir = sys.argv[1] + "/fastq"
isExist = os.path.exists(outdir)
if not isExist:
    os.makedirs(outdir)
outfile = sys.argv[1] + ".fastq"
outfile = outdir + "/" + outfile
print("Output at:",outfile)

with open(outfile,"w") as fw:

    for read in reads:
        if read.read_id in rid_filter:
            fastq = read.fastq
            fw.write(fastq)

