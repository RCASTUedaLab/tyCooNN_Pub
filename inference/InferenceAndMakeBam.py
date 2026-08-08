import pandas as pd
from tensorflow import keras
import nnmodels.CNNWavenet as cnnwavenet
import multiprocessing
import csv
import numpy as np
import utils.tyUtils as ut
import os
from ont_fast5_api.fast5_interface import get_fast5_file
from inference.ExCounter import Counter
from inference.ExCounter import MiniCounter
import preprocess.TrimAndNormalize as tn
import tensorflow as tf
import pathlib


def getTRNAlist(trnapath):
    trnas = []
    with open(trnapath) as f:
        l = f.readlines()
        for trna in l:
            if len(trna) > 0:
                trna = trna.replace('\n', '')
                trna = trna.replace('\"', '')
                trnas.append(trna)
    return trnas

def getSQ(fasta_file):

    lengths_list = []
    current_name = None
    current_sequence = ""

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    lengths_list.append({'LN': len(current_sequence), 'SN': current_name})
                current_name = line[1:].split()[0]
                current_sequence = ""
            else:
                current_sequence += line

        if current_name:
            lengths_list.append({'LN': len(current_sequence), 'SN': current_name})

    print(lengths_list)
    return lengths_list

def getHeader(fasta_file):

    header = {'HD': {'VN': '1.0'}}
    header['SQ'] = getSQ(fasta_file)
    return header


def toUniqueName(tRNA):

    tRNA = tRNA.replace("_rcc","").replace("_ivt","")
    return tRNA

def splitBytRNA(alldata,threshold):

    ret_dict = {}
    for minic in alldata:

        tRNA = minic.tRNA
        if tRNA is None:
            continue
        if minic.maxval < threshold:
            continue

        #tRNAIdx = minic.tRNAIdx
        trna = toUniqueName(tRNA)
        if "spike_in" in trna:
            continue

        if trna in ret_dict:
            ls = ret_dict[trna]
        else:
            ls = []
            ret_dict[trna] = ls
        ls.append(minic.fastq)


    return ret_dict


import pysam
import mappy as mp
maxcnt = 1000
def maptoref(bam,ref,refdir,alldata,threshold):

    bam_file = pysam.AlignmentFile(bam, 'wb', header=getHeader(ref))
    readdict = splitBytRNA(alldata,threshold)

    for trna in readdict.keys():

        if "_rcc" in trna:
            trna = trna.replace("_rcc","")

        ref = refdir + "/" + trna + ".fasta"

        aligner = mp.Aligner(ref, n_threads=10, min_dp_score=15, w=4, bw=1, k=1,
                             best_n=1, min_cnt=1, min_chain_score=1)  # load or build index

        cnt = 0
        mappedc = 0
        fastqlist = readdict[trna]
        for fastq in fastqlist:
            fqs = fastq.split("\n")
            sequence = fqs[1]
            quality = fqs[3]

            for hit in aligner.map(sequence, MD=True):

                if hit.strand == -1:
                    continue

                a = pysam.AlignedSegment()
                a.query_name = fqs[0]
                a.query_sequence = sequence[hit.q_st:hit.q_en]
                a.flag = 0
                a.reference_id = bam_file.get_tid(hit.ctg)
                a.reference_start = hit.r_st
                a.mapping_quality = hit.mapq
                a.cigarstring = hit.cigar_str


                a.template_length = hit.r_en - hit.r_st
                a.query_qualities = pysam.qualitystring_to_array(quality[hit.q_st:hit.q_en])
                a.tags = (("NM", hit.NM), ("MD", hit.MD))
                bam_file.write(a)
                mappedc += 1
                break
            cnt += 1
            if cnt % 100 == 0:
                print(cnt, mappedc)
            if mappedc == maxcnt:
                break


    bam_file.close()
    sortbam = bam[0:len(bam)-4]+"_sorted.bam"
    pysam.sort('-o',sortbam, bam )
    pysam.index(sortbam)

def extend_dict(signaldict,allsignaldict):

    keys = signaldict.keys()
    for key in keys:
        if key in allsignaldict:
            allsignaldict[key].extend(signaldict[key])
        else:
            allsignaldict[key] = signaldict[key]


def ensure_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

import os.path
import random
def evaluate(paramPath,indirs,outdir,outpath,ref, refdir,threshold=0.75):


    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

    outweight = outdir + "learent_arg_weight.h5"
    if not os.path.isfile(outweight):
        outweight = outdir + "/learent_arg_weight.h5"

    param = ut.get_parameter(paramPath)
    # indirs = indirs.split(",")
    f5list = []
    for dir in indirs:
        f5list.extend(ut.get_fast5_files_in_dir(dir, param.ncore))

    random.shuffle(f5list)

    # maxf5 = 30
    # f5list = f5list[0:maxf5]

    trnapath = outdir + '/tRNAindex.csv'
    trnas = getTRNAlist(trnapath)
    print(len(f5list))
    print("trna", trnas)
    model = cnnwavenet.build_network(shape=(None, param.trimlen, 1), num_classes=len(trnas))
    model.load_weights(outweight)
    print("finish load wight")

    totalcounter = Counter(trnas, threshold=threshold)
    cnt = 0
    fqpath = outpath + "/trna.fastq"
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    fq = open(fqpath, mode='w')
    alldata = []

    allsignaldict = {}
    for f5file in f5list:

        counter,datadict,signaldict =  evaluateEach(param, f5file, outpath, model, trnas, cnt, threshold)
        extend_dict(signaldict,allsignaldict)

        datalist = list(datadict.values())
        alldata.extend(datalist)

        totalcounter.sumup(counter)
        cnt += 1
        print(f5file)
        print("done..{}/{}".format(cnt, len(f5list)))


    fq.close()

    # output result
    csvout = outpath + "/count.csv"
    data = []
    data.append(totalcounter.passfilterCnt)
    data.append(totalcounter.allCnt)

    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(csvout)
    #
    filtercsv = outpath + "/filer.csv"
    data = []
    data.append(totalcounter.filterFlgCnt)
    filterlabel = ["pass", "meanqlow", "siglen", "adap2fail", "deltalow", "delthigh", "readlenlow", "readlenhigh",
                   "trimfail"]
    df = pd.DataFrame(data, columns=filterlabel)
    df.to_csv(filtercsv)

    #write parquet
    keys = allsignaldict.keys()
    for k in keys:
        dlist = allsignaldict[k]
        pqpathdir = outpath + "/pq/"
        ensure_directory(pqpathdir)
        pqpath = pqpathdir+k+".pq"
        df = pd.DataFrame(dlist, columns=['readid', 'signal'])
        df.to_parquet(pqpath)

    bamfile = outpath + "/sample.bam"
    maptoref(bamfile,ref,refdir,alldata,threshold)



from Bio import SeqIO


def fastaToDict(fasta):
    seqdict = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        seqdict[record.id] = record.seq.replace('U', 'T')

    return seqdict


# do it file by file
def evaluateEach(param, f5file, outpath, model, trnas, cnt_file, threshold):

    print(f5file)
    reads = ut.get_fast5_reads_from_file(f5file)
    trimmed_filterFlgged_read = tn.trimAdaptor(reads, param)

    flagCount = {}
    for read in trimmed_filterFlgged_read:
        flag = read.filterFlg
        if flag in flagCount:
            flagCount[flag] += 1
        else:
            flagCount[flag] = 1
    flagcsv = outpath + "/flag.csv"
    filepath = pathlib.Path(flagcsv)
    if filepath.exists():
        if cnt_file == 0:
            fileflag = open(filepath, 'w')
        else:
            fileflag = open(filepath, 'a')
    else:
        fileflag = open(filepath, 'w')
    filterlabel = ["pass", "meanqlow", "siglen", "adap2fail", "deltalow", "delthigh", "readlenlow", "readlenhigh",
                   "trimfail"]
    for flg in sorted(flagCount.keys()):
        flg_lab = filterlabel[flg]
        flg_cnt = flagCount[flg]
        fileflag.write("%4d %2d %12s %12d\n" % (cnt_file, flg, flg_lab, flg_cnt))
    fileflag.close()
    # print(flagCount)
    format_reads = tn.formatSignal(trimmed_filterFlgged_read, param)
    datalabel = []
    data = []
    datadict = {}

    for read in format_reads:

        # print(read.read_id)
        if (read.filterFlg == 0):
            datadict[read.read_id] = MiniCounter(read.filterFlg, read.trimSuccess,read.fastq)
            datalabel.append(read.read_id)
            data.append(read.formatSignal)

    print("Number of trimmed reads: ", len(datalabel))

    data_r = np.reshape(data, (-1, param.trimlen, 1))
    prediction = model.predict(data_r, batch_size=None, verbose=0, steps=None)
    print(data_r.shape, prediction.shape)

    signaldict = {}
    cnt = -1
    for row in prediction:

        # incriment
        cnt += 1
        rdata = np.array(row)
        maxidxs = np.where(rdata == rdata.max())
        # unique hit with more than zero Intensity
        if len(maxidxs) == 1 and rdata.max() >= 0:

            maxidx = int(maxidxs[0])
            maxv = rdata.max()
            maxtrna = trnas[maxidx]
            readid = datalabel[cnt]
            minicnt = datadict[readid]
            minicnt.addInference(maxtrna, maxidx, maxv)

            if rdata.max() > threshold:

                if maxtrna in signaldict:
                    dlist  = signaldict[maxtrna]
                    dlist.append((readid,data[cnt]))
                else:
                    dlist = []
                    dlist.append((readid, data[cnt]))
                    signaldict[maxtrna] = dlist

                #
    counter = Counter(trnas, threshold=threshold)
    for key in datadict:
        minicnt = datadict[key]
        counter.inc(minicnt)


    return counter,datadict,signaldict


def getDummyQual(seqlen):
    return ''.join(['A' for i in range(seqlen)])


def getFastq(read_id, seqdict, tRNA, seqlen):
    if tRNA not in seqdict:
        # print(tRNA)
        return None

    seq = seqdict[tRNA]
    if len(seq) <= seqlen:
        seqlen = len(seq)

    hang = 5
    start = (len(seq) - seqlen) - hang
    if start < 0:
        start = 0
    # seq = seq[start:len(seq)]
    qual = getDummyQual(len(seq))
    fq = str(read_id) + " \n" + str(seq) + "\n" + "+" + "\n" + str(qual)
    # print(fq)
    return fq


import logging
import os
import shutil
from ont_fast5_api.fast5_file import Fast5File, Fast5FileTypeError
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.compression_settings import GZIP
import ont_fast5_api.conversion_tools.multi_to_single_fast5 as multi_to_single_fast5
import h5py
import sys

if sys.version_info[0] > 2:
    unicode = str

import time


