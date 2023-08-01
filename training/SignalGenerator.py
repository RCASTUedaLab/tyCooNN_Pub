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
import training.DataArgumentation as da


def shuffle_samples(X, y):

    # zipped = list(zip(X, y))
    # np.random.shuffle(zipped)
    # X_result, y_result = zip(*zipped)
    # return np.asarray(X_result), np.asarray(y_result)

    data_size = len(X)
    shuffle_indices = np.random.permutation(np.arange(data_size))
    shuffled_data = X[shuffle_indices]
    shuffled_labels = y[shuffle_indices]
    return shuffled_data,shuffled_labels

def formatX(X,wlen):
   return np.reshape(X, (-1, wlen, 1))

def formatY(Y,num_classes):
   Y = np.reshape(Y, (-1, 1,))
   return keras.utils.to_categorical(Y, num_classes)

import matplotlib.pyplot as plot
class ArgumentlGenerator(object):

    def __init__(self,x, y, batch_size,signal_size, class_count, augmentation_factor,epoch):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.augmentation_factor = augmentation_factor
        self.epoch = epoch

    def numbatch(self):

        return  int((len(self.x)*self.augmentation_factor - 1) / self.batch_size) + 1

    def flow(self):

        for n in range(self.epoch+1):
            augmented_signals, augmented_labels \
                = da.augment_data(self.x, self.y,self.signal_size, self.augmentation_factor)
            #augmented_signals, augmented_labels = shuffle_samples(augmented_signals, augmented_labels)


            num_batches_per_epoch = self.numbatch()
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * self.batch_size
                end_index = min((batch_num + 1) * self.batch_size, len(augmented_signals))
                batch_X = augmented_signals[start_index: end_index]
                batch_Y = augmented_labels[start_index: end_index]

                # plot.plot(batch_X[0])
                # plot.title(str(batch_Y[0]))
                # plot.savefig("/share/trna/tyCooNNTest/testfig/epoch"+str(n)+"_"+str(batch_num)+".png")
                # plot.clf()

                batch_X = formatX(batch_X,self.signal_size)
                batch_Y = formatY(batch_Y,self.class_count)
                yield (batch_X,batch_Y)



class BatchIterator(object):

    def __init__(self,x, y, batch_size,signal_size,class_count,epoch):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.epoch = epoch


    def numbatch(self):
        return int((len(self.x) - 1) / self.batch_size) + 1

    def flow(self):

        for n in range(self.epoch+1):
            num_batches_per_epoch = self.numbatch()
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * self.batch_size
                end_index = min((batch_num + 1) * self.batch_size, len(self.x))
                batch_X = self.x[start_index: end_index]
                batch_Y = self.y[start_index: end_index]
                batch_X = formatX(batch_X,self.signal_size)
                batch_Y = formatY(batch_Y,self.class_count)
                yield (batch_X,batch_Y)