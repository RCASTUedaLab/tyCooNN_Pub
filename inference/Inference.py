import os
import os.path
import h5py
import sys
import pathlib
import shutil
import glob
from Bio import SeqIO
import numpy as np
import pandas as pd
from pyarrow import parquet as pq
from ont_fast5_api.multi_fast5 import MultiFast5File
import ont_fast5_api.conversion_tools.multi_to_single_fast5 as multi_to_single_fast5
import utils.tyUtils as ut
import nnmodels.CNNWavenet as cnnwavenet
from inference.ExCounter import Counter
from inference.ExCounter import MiniCounter
import preprocess.TrimAndNormalize as tn
if sys.version_info[0] > 2:
    unicode = str

def getTRNAlist(trnapath):

    trnas = []
    with open(trnapath) as f:
        l = f.readlines()
        for trna in l:
            if len(trna) > 0:
                trna = trna.replace('\n','')
                trna = trna.replace('\"', '')
                trnas.append(trna)
    return trnas

def infer(input, modeldir, outpath, ref, fast5fmt, threshold, parampath):

    print("Softmax post-filter threshold: ",threshold)
    
    modelweight = modeldir + "learent_arg_weight.h5"
    if not os.path.isfile(modelweight):
        print("No model with augmented training: fall back to learent_weight.h5")
        modelweight = modeldir + "/learent_weight.h5"

    param = ut.get_parameter(parampath)
    input = input.split(",")
    f5list = []
    for dir in input:
        f5list.extend(ut.get_fast5_files_in_dir(dir,param.ncore))
    print("Number of fast5 files: %d" % len(f5list))

    trnapath = modeldir + '/tRNAindex.csv'
    trnas = getTRNAlist(trnapath)
    print("tRNAs:\n",np.array(trnas))

    model = cnnwavenet.build_network(shape=(None, param.trimlen, 1), num_classes=len(trnas))
    model.load_weights(modelweight)

    totalcounter = Counter(trnas,threshold=threshold)
    cnt = 0
    fqpath = outpath + "/trna.fastq"
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    fq = open(fqpath, mode='w')
    for f5file in f5list:
        counter = evaluateEach(param,f5file,outpath,model,trnas,ref,fast5fmt,cnt,fq,threshold)
        totalcounter.sumup(counter)
        cnt +=1
        print("done..{}/{}".format(cnt,len(f5list)))

    fq.close()

    #output result
    csvout = outpath + "/count.csv"
    data = []
    data.append(totalcounter.passfilterCnt)
    data.append(totalcounter.allCnt)

    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(csvout)
    
    filtercsv = outpath + "/filer.csv"
    data = []
    data.append(totalcounter.filterFlgCnt)
    filterlabel = ["pass","meanqlow","siglen","adap2fail","deltalow","delthigh","readlenlow","readlenhigh","trimfail"]
    df = pd.DataFrame(data, columns=filterlabel)
    df.to_csv(filtercsv)

def fastaToDict(fasta):

    seqdict = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        seqdict[record.id]  = record.seq.replace('U','T')

    return seqdict

# do it file by file
def evaluateEach(param,f5file,outpath,model,trnas,ref,fast5fmt,cnt_file,fq,threshold):

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
            fileflag = open(filepath,'w')
        else:
            fileflag = open(filepath,'a')
    else:
        fileflag = open(filepath,'w')
    filterlabel = ["pass","meanqlow","siglen","adap2fail","deltalow","delthigh","readlenlow","readlenhigh","trimfail"]
    for flg in sorted(flagCount.keys()):
        flg_lab = filterlabel[flg]
        flg_cnt = flagCount[flg]
        fileflag.write("%4d %2d %12s %12d\n" % (cnt_file,flg,flg_lab,flg_cnt))
    fileflag.close()
    #print(flagCount)

    format_reads = tn.formatSignal(trimmed_filterFlgged_read, param)
    datalabel = []
    data = []
    datadict = {}

    seqdict = fastaToDict(ref)

    fast5dir = outpath +"/fast5"
    if not os.path.exists(fast5dir):
        os.makedirs(fast5dir)
    fast5out = fast5dir+"/"+  os.path.basename(f5file)

    for read in format_reads:

        if (read.filterFlg == 0):
            datadict[read.read_id] = MiniCounter(read.filterFlg,read.trimSuccess)
            datalabel.append(read.read_id)
            data.append(read.formatSignal)

    print("Number of trimmed reads: ",len(datalabel))
    data = np.reshape(data, (-1, param.trimlen, 1))
    prediction = model.predict(data, batch_size=None, verbose=0, steps=None)
    print(data.shape,prediction.shape)

    cnt = -1
    for row in prediction:

        cnt += 1
        rdata = np.array(row)
        maxidxs = np.where(rdata == rdata.max())
        #unique hit with more than zero Intensity
        if len(maxidxs) == 1 and rdata.max() >= 0:
            maxidx = int(maxidxs[0])
            maxv = rdata.max()
            maxtrna = trnas[maxidx]
            readid = datalabel[cnt]
            minicnt =  datadict[readid]
            minicnt.addInference(maxtrna,maxidx,maxv)

    counter = Counter(trnas,threshold=threshold)
    for key in datadict:
        minicnt = datadict[key]
        counter.inc(minicnt)

    singlefast5dir = outpath + "/single_fast5"
    #output fast5
    copyWithAdddata(f5file,fast5out,datadict,seqdict,fast5fmt,singlefast5dir,cnt_file,fq)

    return counter

def getDummyQual(seqlen):

    return ''.join(['A' for i in range(seqlen)])

def getFastq(read_id,seqdict,tRNA,seqlen):

    if tRNA not in seqdict:
        #print(tRNA)
        return None

    seq = seqdict[tRNA]
    if len(seq) <= seqlen:
        seqlen = len(seq)

    hang = 5
    start = (len(seq)-seqlen)-hang
    if start < 0:
        start = 0
    qual = getDummyQual(len(seq))
    fq = str(read_id)+ " \n"  + str(seq) +"\n" +"+" + "\n" + str(qual)
    return fq

def copyWithAdddata(f5file,fast5out,datadict,seqdict,fast5fmt,singlefast5dir,cnt,fq):

    #copy first
    shutil.copyfile(f5file, fast5out)

    with MultiFast5File(fast5out, 'a') as multi_f5:
        rcnt = -1
        for read in multi_f5.get_reads():

            rcnt += 1
            component = "basecall_1d"
            group_name = "Basecall_1D_099"
            dataset_name = "BaseCalled_template"

            basecall_run = read.get_latest_analysis("Basecall_1D")
            fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
            seqlen = len(fastq.split("\n")[1])

            if read.read_id in datadict:

                minicnt = datadict[read.read_id]
                fstline = fastq.split("\n")[0]
                fastqadd = getFastq(fstline,seqdict, minicnt.tRNA, seqlen)

                if fastqadd is not None:

                    fq.write(fastqadd)
                    fq.write("\n")

                    attrs = {
                        "tRNA": minicnt.tRNA,
                        "tRNAIndex": minicnt.tRNAIdx,
                        "value": minicnt.maxval,
                        "filterpass": (minicnt.filterFlg == 0),
                        "filterflg": minicnt.filterFlg,
                        "trimSuccess": minicnt.trimSuccess
                    }
                    read.add_analysis(component, group_name, attrs)
                    path = 'Analyses/{}/'.format(group_name)
                    read.handle[path].create_group(dataset_name)
                    path = 'Analyses/{}/{}'.format(group_name, dataset_name)

                    read.handle[path].create_dataset(
                        'Fastq', data=str(fastqadd),
                        dtype=h5py.special_dtype(vlen=unicode))


    multi_f5.close()

    if fast5fmt == "S":
        print('print single5 output to',singlefast5dir,str(cnt+1))
        multi_to_single_fast5.convert_multi_to_single(fast5out, singlefast5dir,str(cnt+1))

