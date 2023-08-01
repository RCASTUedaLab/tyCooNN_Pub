from pyarrow import parquet as pq
from itertools import chain
import numpy as np
from multiprocessing import Pool
import random
import glob
import pandas as pd
from tensorflow import keras
import nnmodels.CNNWavenet as cnnwavenet
from numba import jit
import multiprocessing
import csv
import numpy as np
import os.path as path
import os
import tensorflow as tf
import pickle

def prepare_data(inp_prepare_data):

    fid,filename = inp_prepare_data

    X_test = []
    Y_test = []
    trnas  = []
    pqt = pq.read_table(filename,
                        columns=['read_id', 'trna','trimsignal'])
    print("Load %3d: %s" % (fid+1,filename))
    dfp = pqt.to_pandas()
    cnt = 0
    wlen = 0
    for idx, row in dfp.iterrows():
        trna = row[1]
        signal = row[2]
        if wlen ==0:
            wlen = len(signal)

        if not cnt % 12 >= 2:
            X_test.append(signal)
            Y_test.append(trna)

        cnt+=1

    trna = list(dfp["trna"].unique())
    trnas.append(trna)

    return {'trna': trnas, 'wlen': wlen, 'x': X_test, 'y': Y_test }


def evaluate(opts):

    dirpath = opts['inp_loc']
    outdir = opts['model_loc']
    csvout = opts['csvout']
    csvout2 = opts['csvout2']
    csvdir = os.path.split(csvout)[0]

    if 'threshold' in opts:
        threshold=opts['threshold']
    else:
        threshold=0.75
    if 'max_core' in opts:
        max_core = opts['max_core']
    else:
        max_core = 4

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

    use_mult_gpu = False
    use_gpu = False
    if 'gpu' in opts:
        gpu_select = str(opts['gpu'])
        print("Using %s" % gpu_select)
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_select
        if gpu_select == '':
            use_gpu = False
        else:
            use_gpu = True
            num_gpu = len(gpu_select.split(','))
            if num_gpu > 1:
                use_mult_gpu = True
    if 'gpu_memory_limit' in opts:
        gpu_memory_limit = 1024 * opts['gpu_memory_limit']
        gpu_logical_set = True

    else:
        gpu_logical_set = False

    # Setting computing device
    if use_gpu:
        if use_mult_gpu:
            gpus = tf.config.list_physical_devices('GPU')
            if gpu_logical_set:
                for gpu in gpus:
                    tf.config.set_logical_device_configuration(gpu,[tf.config.LogicalDeviceConfiguration(memory_limit=gpu_memory_limit)])
                gpus = tf.config.list_logical_devices('GPU')
            else:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
        else:
            if gpu_logical_set:
                gpus = tf.config.list_physical_devices('GPU')
                tf.config.set_logical_device_configuration(gpus[0],[tf.config.LogicalDeviceConfiguration(memory_limit=gpu_memory_limit)])



    fs = sorted(glob.glob(dirpath + "/*.pq*"))
    trnas = []

    X_test = []
    Y_test = []

    wlen = 0
    with Pool(processes=max_core) as pool:
        inp_prepare = [(fid,f) for fid,f in enumerate(fs)]
        result = pool.map(prepare_data,inp_prepare)
        for r in result:
            trnas.append(r['trna'][0])
            X_test.extend(r['x'])
            Y_test.extend(r['y'])
            wlen = r['wlen']

    trnas = list(chain.from_iterable(trnas))
    trnas = sorted(trnas)
    # name to index
    Y_test = list(map(lambda trna: trnas.index(trna), Y_test))
    num_classes = len(np.unique(Y_test))

    test_x = np.reshape(X_test, (-1, wlen, 1))
    testy_noncategorical = Y_test

    outweight = outdir + "learent_arg_weight.h5"
    if not path.isfile(outweight):
        outweight = outdir + "/learent_arg_weight.h5"
    if use_mult_gpu and gpu_logical_set:
        strategy = tf.distribute.MirroredStrategy(gpus)
        with strategy.scope():
            model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=len(trnas))
            model.load_weights(outweight)
    else:
        model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=len(trnas))
        model.load_weights(outweight)

    prediction = model.predict(test_x, batch_size=None, verbose=0, steps=None)

    with open(csvdir + '/predict_data.pkl', 'wb') as f:
        pickle.dump({'test':test_x,'prediction':prediction}, f,protocol=pickle.HIGHEST_PROTOCOL)

    retdict = {}
    cnt = -1
    prob = []
    for row in prediction:

        cnt += 1
        rdata = np.array(row)
        maxidxs = np.where(rdata == rdata.max())
        prob.append(rdata.max())
        ans = testy_noncategorical[cnt]
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
    print("average prob.: ",np.mean(prob))
    print("std: ",np.std(prob))

    data = [list(retdict[i]) for i in range(num_classes)]

    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(csvout,index=False)

    #probthres = np.mean(prob) #- np.std(prob)
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
        ans = testy_noncategorical[cnt]
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

    print("trimmed average prob.: ",np.mean(prob))
    print("trimmed std: ",np.std(prob))

    data = [list(retdict[i]) for i in range(num_classes)]

    df = pd.DataFrame(data, columns=trnas)
    df.to_csv(csvout2,index=False)



