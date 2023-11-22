from pyarrow import parquet as pq
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
from tensorflow.keras.callbacks import ModelCheckpoint
import utils.tyUtils as ut
from training.SignalGenerator import ArgumentlGenerator
from training.SignalGenerator import BatchIterator
from training.SignalGenerator import TrainDataGenerator,TestDataGenerator,CustomDataGenerator
import tensorflow as tf
import tensorflow_addons as tfa
from pathlib import Path
import os
from tensorflow.keras.utils import Sequence

def load_data(dirpath):

    print(dirpath)
    fs = glob.glob(dirpath + "/*.pq*")
    #fs = fs[0:3] #for debug
    print(fs)
    trnas = []

    X_train = []
    Y_train = []
    X_test = []
    Y_test = []
    wlen = 0
    for f in fs:

        print(f)
        pqt = pq.read_table(f,
                            columns=['read_id', 'trna','trimsignal'])

        dfp = pqt.to_pandas()
        cnt = 0
        wlen = 0
        for idx, row in dfp.iterrows():
            trna = row[1]
            signal = row[2]
            if wlen == 0:
                wlen = len(signal)

            if cnt % 12 >= 2:
                X_train.append(signal)
                Y_train.append(trna)
            else:
                X_test.append(signal)
                Y_test.append(trna)

            cnt+=1

        trna = dfp["trna"].unique()
        trnas.append(trna)

    print("size of train: ",len(X_train))

    trnas = sorted(trnas)
    # name to index
    Y_train = list(map(lambda trna: trnas.index(trna), Y_train))
    Y_test = list(map(lambda trna: trnas.index(trna), Y_test))
    num_classes = np.unique(Y_train).size

    return X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes

def load_structured_data(input):

    print("start")
    trnas = []
    X_train = []
    Y_train = []
    X_test  = []
    Y_test  = []
    wlen = 0
    #limit_train = 2000
    #limit_test  = 400
    for label in input['train'].keys():
        print("Load %s from train" % label)
        input_pq = input['train'][label]
        pqt = pq.read_table(input_pq,columns=['read_id', 'trna','trimsignal'])
        dfp = pqt.to_pandas()
        #dfp = dfp.iloc[:limit_train].copy()
        for idx, row in dfp.iterrows():
            signal = row[2]
            if wlen == 0:
                wlen = len(signal)
            X_train.append(signal)
            Y_train.append(label)
        trnas.append(label)
    for label in input['test'].keys():
        print("Load %s from test" % label)
        input_pq = input['test'][label]
        pqt = pq.read_table(input_pq,columns=['read_id', 'trna','trimsignal'])
        dfp = pqt.to_pandas()
        #dfp = dfp.iloc[:limit_test].copy()
        for idx, row in dfp.iterrows():
            signal = row[2]
            if wlen == 0:
                wlen = len(signal)
            X_test.append(signal)
            Y_test.append(label)
        trnas.append(label)

    trnas = np.unique(trnas)
    print("size of train: ",len(X_train))
    trnas = sorted(trnas)

    # name to index
    Y_train = list(map(lambda trna: trnas.index(trna), Y_train))
    Y_test = list(map(lambda trna: trnas.index(trna), Y_test))
    num_classes = np.unique(Y_train).size

    return X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes

def train_structuredInput(input,outdir,epoch = 50,gpu = None,data_argument = 0,
                          loss_fn='categorical_crossentropy',jumpstart=False,jump_version=1):

    print("Init")
    X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes = load_structured_data(input)
    train_with_input(X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes,outdir,
                     epoch=epoch,gpu=gpu,data_argument=data_argument,loss_fn=loss_fn,
                     jumpstart=jumpstart,jump_version=jump_version)

def train(dirpath,outdir,epoch = 50,gpu= None,data_argument = 0,loss_fn='categorical_crossentropy',
          jumpstart=False,jump_version=1):

    X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes = load_data(dirpath)
    train_with_input(X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes,outdir,epoch=epoch,gpu=gpu,
                     data_argument=data_argument,loss_fn=loss_fn,jumpstart=jumpstart,jump_version=jump_version)

def train_with_input(X_train,Y_train,X_test,Y_test,wlen,trnas,num_classes,outdir,epoch = 50,
                     gpu=None,data_argument = 0,loss_fn = 'categorical_crossentropy',jumpstart=False,jump_version=1):


    if data_argument == 0:

        #lr = 0.0008
        lr = 0.0001
        model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=num_classes)
        optim = keras.optimizers.Adam(learning_rate=lr, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
        optim = tfa.optimizers.SWA(optim)
        batch_size = 128

    if data_argument > 0:
        lr = 0.00001
        model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=num_classes,do_r =0.1)
        optim = keras.optimizers.Adam(learning_rate=lr, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0,
                                      amsgrad=False)
        optim = tfa.optimizers.SWA(optim)
        batch_size = 128


    model.summary()
    model.compile(loss=loss_fn, optimizer=optim, metrics=['accuracy'])

    Path(outdir).mkdir(parents=True, exist_ok=True)
    outweight = outdir + "/learent_weight.h5"
    historypath = outdir + '/history.csv'
    graphpath = outdir + "/learent_graph.png"
    if jumpstart:
        if data_argument > 0:
            outweight = outdir + "/learent_arg_weight.h5"
            model.load_weights(outweight)
            outweight = outdir + "/learent_arg_weight_%d.h5" % jump_version
            historypath = outdir + '/history_arg_%d.csv' % jump_version
            graphpath = outdir + "/learent_arg_graph_%d.png" % jump_version
        else:
            model.load_weights(outweight)
            outweight = outdir + "/learent_weight_%d.h5" % jump_version
            historypath = outdir + '/history_%d.csv' % jump_version
            graphpath = outdir + "/learent_graph_%d.png" % jump_version
    else:
        if data_argument > 0:
            model.load_weights(outweight)
            outweight = outdir + "/learent_arg_weight.h5"
            historypath = outdir + '/history_arg.csv'
            graphpath = outdir + "/learent_arg_graph.png"

    modelCheckpoint = ModelCheckpoint(filepath=outweight,
                                      monitor='val_accuracy',
                                      verbose=1,
                                      save_best_only=True,
                                      save_weights_only=True,
                                      mode='max',
                                      period=1)

    if data_argument == 0:
        
        train_gen = CustomDataGenerator(X_train, Y_train, batch_size, wlen, num_classes,epoch=epoch,shuffle=True)
        test_gen  = CustomDataGenerator(X_test,Y_test,batch_size,wlen,num_classes,epoch=epoch,shuffle=False)

        history = model.fit(train_gen,validation_data=test_gen,epochs=epoch,callbacks=[modelCheckpoint])

    else:

        train_gen = CustomDataGenerator(X_train, Y_train, batch_size, wlen, num_classes,
                                        epoch=epoch,shuffle=False,aug=data_argument)
        test_gen  = CustomDataGenerator(X_test,Y_test,batch_size,wlen,num_classes,epoch=epoch,shuffle=False)

        history = model.fit(train_gen,validation_data=test_gen,epochs=epoch,callbacks=[modelCheckpoint])

    FIG_SIZE_WIDTH = 12
    FIG_SIZE_HEIGHT = 10
    FIG_FONT_SIZE = 25

    #tRNAs
    trnapath = outdir + '/tRNAindex.csv'
    with open(trnapath,'w') as file:
        for t in trnas:
            file.write("%s\n" % t)

    hist_df = pd.DataFrame(history.history)
    hist_df.to_csv(historypath)

    ut.plot_history(history,
                 save_graph_img_path=graphpath,
                 fig_size_width=FIG_SIZE_WIDTH,
                 fig_size_height=FIG_SIZE_HEIGHT,
                 lim_font_size=FIG_FONT_SIZE)


