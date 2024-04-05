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

def evaluate(input, modeldir, outcsv, outcsv2, threshold):

    print(input)
    fs = glob.glob(input + "/*.pq*")
    print(fs)
    trnas = []

    X_test = []
    Y_test = []
    wlen = 0
    for f in fs:

        print(f)
        pqt = pq.read_table(f,
                            columns=['read_id', 'trna', 'trimsignal'])

        dfp = pqt.to_pandas()
        cnt = 0
        wlen = 0
        for idx, row in dfp.iterrows():
            trna = row[1]
            signal = row[2]
            if wlen == 0:
                wlen = len(signal)

            if not cnt % 12 >= 2:
                X_test.append(signal)
                Y_test.append(trna)

            cnt += 1

        trna = dfp["trna"].unique()
        trnas.append(trna)

    trnas = sorted(trnas)
    # name to index
    Y_test = list(map(lambda trna: trnas.index(trna), Y_test))
    num_classes = np.unique(Y_test).size

    test_x = np.reshape(X_test, (-1, wlen, 1))

    modelweight = modeldir + "learent_arg_weight.h5"
    if not os.path.isfile(modelweight):
        print("No model with augmented training: fall back to learent_weight.h5")
        modelweight = modeldir + "/learent_weight.h5"

    model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=len(trnas))
    model.load_weights(modelweight)

    prediction = model.predict(test_x, batch_size=None, verbose=0, steps=None)

    retdict = {}
    cnt = -1
    prob = []
    for row in prediction:

        cnt += 1
        rdata = np.array(row)
        maxidxs = np.where(rdata == rdata.max())
        prob.append(rdata.max())
        ans = Y_test[cnt]
        if len(maxidxs) > 1:
            continue  # multiple hit
        maxidx = maxidxs[0]

        if ans in retdict:
            ridxs = retdict[ans]
            ridxs[maxidx] = ridxs[maxidx] + 1
        else:
            ridxs = np.zeros(num_classes)
            retdict[ans] = ridxs
            ridxs[maxidx] = ridxs[maxidx] + 1

    print("average prob.: ", np.mean(prob))
    print("std: ", np.std(prob))

    data = [list(retdict[i]) for i in range(num_classes)]

    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(outcsv, index=False)

    probthres = threshold
    retdict = {}
    cnt = -1
    prob = []
    for row in prediction:

        cnt += 1
        rdata = np.array(row)
        if rdata.max() < probthres:
            continue
        maxidxs = np.where(rdata == rdata.max())
        prob.append(rdata.max())
        ans = Y_test[cnt]
        if len(maxidxs) > 1:
            continue  # multiple hit
        maxidx = maxidxs[0]

        if ans in retdict:
            ridxs = retdict[ans]
            ridxs[maxidx] = ridxs[maxidx] + 1
        else:
            ridxs = np.zeros(num_classes)
            retdict[ans] = ridxs
            ridxs[maxidx] = ridxs[maxidx] + 1

    print("trimmed average prob.: ", np.mean(prob))
    print("trimmed std: ", np.std(prob))

    data = [list(retdict[i]) for i in range(num_classes)]
    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(outcsv2, index=False)

def evaluatepq(input, modeldir, outcsv, threshold):

    modelweight = modeldir + "learent_arg_weight.h5"
    if not os.path.isfile(modelweight):
        modelweight = modeldir + "/learent_weight.h5"

    data = pd.read_parquet(input)
    signal_data = data['trimsignal'].to_list()

    wlen = len(signal_data[0])

    trnapath = modeldir + '/tRNAindex.csv'
    trnas = getTRNAlist(trnapath)
    print("trna", np.array(trnas))

    model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=len(trnas))
    model.load_weights(modelweight)

    data = np.reshape(signal_data, (-1, wlen, 1))
    prediction = model.predict(data, batch_size=None, verbose=0, steps=None)
    print(data.shape, prediction.shape)

    label = {}
    label_threshold = {}
    for t in trnas:
        label[t] = 0
        label_threshold[t] = 0

    for row in prediction:

        rdata = np.array(row)
        maxidxs = np.where(rdata == rdata.max())
        # unique hit with more than zero Intensity
        if len(maxidxs) == 1 and rdata.max() >= 0:
            maxidx = int(maxidxs[0])
            maxtrna = trnas[maxidx]
            maxv = rdata.max()
            label[maxtrna] += 1
            if maxv >= threshold:
                label_threshold[maxtrna] += 1

    with open(outcsv, "w") as fw:
        header_string = "," + ",".join(sorted(trnas))
        fw.write("%s\n" % header_string)
        values = [];
        values_cut = []
        for k in sorted(label.keys()):
            values.append(str(label[k]))
            values_cut.append(str(label_threshold[k]))
        v1 = "," + ",".join(values)
        v2 = "," + ",".join(values_cut)
        fw.write("%s\n%s\n" % (v2, v1))

