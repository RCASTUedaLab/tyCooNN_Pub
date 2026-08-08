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

    """
    - Inference on any parquet file (just for a handy feature)

    note this will not generate any fast5 file however it can be applied
    to a pre-made parquet file


    python tyCooNN.py handyevaluate

       input:     input directory for fast5 files
       modeldir:  location of trained model
       outcsv:    file-name for abundances of species
       threshold: post-filter threshold

    """

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

