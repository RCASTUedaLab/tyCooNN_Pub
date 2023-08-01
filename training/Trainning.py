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
import tensorflow as tf
import tensorflow_addons as tfa

# @click.option('--inp')
# @click.option('--out')
# @click.option('--extn', default='.pq')
# @click.option('--test', type=float, default=0.1)
# @click.option('--seed', type=int, default=100)
# @click.option('--limit', type=int, default=10000)
# @click.option('--ngpu', type=int, default=2)
# @click.option('--epoch', type=int, default=500)
# @click.option('--batch', type=int, default=64)
# def train(inp, out, extn, test, seed, limit, ngpu, epoch, batch):
#     CallTrainer(inp, out, extn=extn, testFraction=test, seed=seed, limit=limit,
#                 wlen=4096, pseudocount=0.1, gpu_count=ngpu, epoch=epoch, batch_size=batch)


def prepare_data(df_inp,trnas):

    X = []
    Xtrace = []
    Y = []

    # pqt = pq.read_table(f,
    #   columns=['read_id', 'trna','meanq','normdelta','countCp','trimsignal', 'trimtrace','fastq'])

    for idx, row in df_inp.iterrows():

        trna = row[1]
        signal = row[5]
        trace = row[6]

        index = trnas.index(trna)

        X.append(signal)
        Xtrace.append(trace)
        Y.append(index)

    X = np.array(X)
    Y = np.array(Y)
    return X,Xtrace,Y

def formatX(X,wlen):
   return np.reshape(X, (-1, wlen, 1))

def formatY(Y,num_classes):
   Y = np.reshape(Y, (-1, 1,))
   return keras.utils.to_categorical(Y, num_classes)


def train(dirpath,outdir,epoch = 50,data_argument = 0):

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

    trnas = sorted(trnas)
    # name to index
    Y_train = list(map(lambda trna: trnas.index(trna), Y_train))
    Y_test = list(map(lambda trna: trnas.index(trna), Y_test))
    num_classes = np.unique(Y_train).size
    #

    if data_argument == 0:

        lr = 0.0008
        model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=num_classes)
        optim = keras.optimizers.Adam(learning_rate=lr, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
        batch_size = 256

    if data_argument > 0:
        lr = 0.00001
        model = cnnwavenet.build_network(shape=(None, wlen, 1), num_classes=num_classes,do_r =0.1)
        optim = keras.optimizers.Adam(learning_rate=lr, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0,
                                      amsgrad=False)
        optim = tfa.optimizers.SWA(optim)
        batch_size = 128

    model.summary()
    model.compile(loss='categorical_crossentropy', optimizer=optim, metrics=['accuracy'])


    outweight = outdir + "/learent_weight.h5"
    historypath = outdir + '/history.csv'
    graphpath = outdir + "/learent_graph.png"
    if data_argument > 0:
        # do training without data argumentation 50 epoch
        # then data argumentation with smaller learning rate
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


    test_x = formatX(X_test, wlen)
    test_y = formatY(Y_test, num_classes)
    if data_argument == 0:

        train_x = formatX(X_train,wlen)
        train_y = formatY(Y_train,num_classes)
        print('test_x shape:', test_x.shape)
        print('test_y shape:', test_y.shape)

        history = model.fit(train_x, train_y, epochs=epoch, batch_size=batch_size, verbose=1,
                  shuffle=True, validation_data=(test_x, test_y),callbacks=[modelCheckpoint])
    else:

        signalgen = ArgumentlGenerator(X_train, Y_train, batch_size,wlen,num_classes,data_argument,epoch)
        batchgen = BatchIterator(X_test, Y_test, batch_size,wlen,num_classes,epoch)
        # history = model.fit(signalgen.flow(),steps_per_epoch=signalgen.numbatch(), verbose=1,epochs=epoch,
        #           shuffle=False, validation_data=(test_x, test_y),callbacks=[modelCheckpoint])
        history = model.fit_generator(signalgen.flow(),steps_per_epoch=signalgen.numbatch(),
                                      validation_data=batchgen.flow(),validation_steps=batchgen.numbatch(),epochs=epoch,callbacks=[modelCheckpoint])

    FIG_SIZE_WIDTH = 12
    FIG_SIZE_HEIGHT = 10
    FIG_FONT_SIZE = 25

    #tRNAs
    trnapath = outdir + '/tRNAindex.csv'
    with open(trnapath, 'w', newline='') as file:
        writer = csv.writer(file, quoting=csv.QUOTE_ALL, delimiter=',')
        writer.writerows(trnas)

    hist_df = pd.DataFrame(history.history)
    hist_df.to_csv(historypath)

    ut.plot_history(history,
                 save_graph_img_path=graphpath,
                 fig_size_width=FIG_SIZE_WIDTH,
                 fig_size_height=FIG_SIZE_HEIGHT,
                 lim_font_size=FIG_FONT_SIZE)



