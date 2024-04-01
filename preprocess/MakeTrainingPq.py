import preprocess.TrimAndNormalize as tn
import utils.tyUtils as ut
import pandas as pd
import csv
from collections import Counter

def generatePqForTrainingAll(paramPath,listOfIOPath,takeCount=12000):

    f = open(listOfIOPath, 'r')
    tsv = csv.reader(f, delimiter='\t')

    # tRNAlabel indir   outpq
    for row in tsv:

        tRNALabel, indirs, outpq, outstat = row[0],row[1],row[2],row[3]
        outpq = outpq.strip(" ")
        outpq = outpq.strip("\t")
        print("doing..",tRNALabel,outpq)
        genaratePqForTraining(paramPath, tRNALabel, indirs, outpq, outstat,takeCount)

def genaratePqForTraining(paramPath,tRNALabel,indirs,outpq,outstat,takeCount=12000):

    param = ut.get_parameter(paramPath)  # setting of file path and max core

    indirs = indirs.split(",")
    reads = ut.get_fast5_reads_dirs(indirs, param.ncore)
    
    trimmed_filterFlgged_read = tn.trimAdaptor(reads,param)

    print_trim_stat(trimmed_filterFlgged_read,tRNALabel,outstat)

    #
    filtered_reads = [read for read in trimmed_filterFlgged_read \
                      if read.filterFlg == 0]
    if takeCount is not None:
        filtered_reads = filtered_reads[0:takeCount]
    format_reads = tn.formatSignal(filtered_reads,param)
    #
    datalist = []
    for read in format_reads:
        tp = (read.read_id,tRNALabel,read.formatSignal)
        datalist.append(tp)

    df = pd.DataFrame(datalist, columns=['read_id', 'trna', 'trimsignal'])
    df.to_parquet(outpq)

    return reads

def print_trim_stat(reads,tname,outstat):

    fs = open(outstat,'a')
    flags = [read.filterFlg for read in reads]
    flag_dic = dict(Counter(flags))
    flag_name = {0:'pass',1:'MeanQ',2:'Maxsignallen',3:'MaxdurationRate',
                 4:'DeltaMin',5:'DeltaMax',6:'ReadlenMin',7:'ReadlenMax',
                 8:'TraimFail'}
    fresult = [0] * 9
    for fkey in sorted(flag_name.keys()):
        if fkey in flag_dic:
            fcount = flag_dic[fkey]
            fresult[fkey] = fcount
    tot = sum(fresult)
    fs.write("STAT ] %s: N: %d " % (tname,tot))
    for fkey in sorted(flag_name.keys()):
        fs.write("%s: %d " % (flag_name[fkey],fresult[fkey]))
    fs.write("\n")
    fs.close()
