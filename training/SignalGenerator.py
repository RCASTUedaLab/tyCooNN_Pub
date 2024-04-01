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
import training.DataAugmentation as da
from tensorflow.keras.utils import Sequence
import matplotlib.pyplot as plot
import gc

def shuffle_samples(X, y):

    # zipped = list(zip(X, y))
    # np.random.shuffle(zipped)
    # X_result, y_result = zip(*zipped)
    # return np.asarray(X_result), np.asarray(y_result)

    data_size = len(X)
    shuffle_indices = np.random.permutation(np.arange(data_size))
    X = np.array(X)
    y = np.array(y)
    shuffled_data = X[shuffle_indices]
    shuffled_labels = y[shuffle_indices]
    return shuffled_data,shuffled_labels

def formatX(X,wlen):
   return np.reshape(X, (-1, wlen, 1))

def formatY(Y,num_classes):
   Y = np.reshape(Y, (-1, 1,))
   return keras.utils.to_categorical(Y, num_classes)

class AugmentationGenerator(object):

    def __init__(self,x, y, batch_size,signal_size, class_count, augmentation_factor,epoch,ncore=8):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.augmentation_factor = augmentation_factor
        self.epoch = epoch
        self.ncore = ncore

    def numbatch(self):

        return  int((len(self.x)*self.augmentation_factor - 1) / self.batch_size) + 1

    def flow(self):

        for n in range(self.epoch):

            gc.collect()

            augmented_signals, augmented_labels \
                = da.augment_data(self.x, self.y,self.signal_size, self.augmentation_factor,self.ncore)
                #= da.augment_data_shm(self.x, self.y,self.signal_size, self.augmentation_factor,self.ncore)
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

class TrainDataGenerator(object):

    def __init__(self,x, y, batch_size,signal_size, class_count,epoch):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.epoch = epoch
        print("Number of training data ",len(self.x))
        print("Number of training batches ",self.numbatch())

    def numbatch(self):

        #return  int((len(self.x) - 1) / self.batch_size) + 1
        return  int(np.ceil(len(self.x)/self.batch_size))

    def flow(self):

        for n in range(self.epoch+1):
            X,Y = shuffle_samples(self.x,self.y)
            num_batches_per_epoch = self.numbatch()
            print("Training batches per epoch:",num_batches_per_epoch)
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * self.batch_size
                end_index = min((batch_num + 1) * self.batch_size, len(self.x))
                batch_X = X[start_index: end_index]
                batch_Y = Y[start_index: end_index]

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

        for n in range(self.epoch):
            num_batches_per_epoch = self.numbatch()
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * self.batch_size
                end_index = min((batch_num + 1) * self.batch_size, len(self.x))
                batch_X = self.x[start_index: end_index]
                batch_Y = self.y[start_index: end_index]
                batch_X = formatX(batch_X,self.signal_size)
                batch_Y = formatY(batch_Y,self.class_count)
                yield (batch_X,batch_Y)

class TestDataGenerator(object):

    def __init__(self,x, y, batch_size,signal_size,class_count,epoch):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.epoch = epoch
        print("Number of validaion data ",len(self.x))
        print("Number of training batches ",self.numbatch())

    def numbatch(self):
        #return int((len(self.x) - 1) / self.batch_size) + 1
        return  int(np.ceil(len(self.x)/self.batch_size))

    def flow(self):

        for n in range(self.epoch+1):
            num_batches_per_epoch = self.numbatch()
            print("Test batches per epoch:",num_batches_per_epoch)
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * self.batch_size
                end_index = min((batch_num + 1) * self.batch_size, len(self.x))
                batch_X = self.x[start_index: end_index]
                batch_Y = self.y[start_index: end_index]
                batch_X = formatX(batch_X,self.signal_size)
                batch_Y = formatY(batch_Y,self.class_count)
                yield (batch_X,batch_Y)


class CustomDataGenerator(Sequence):

    def __init__(self,x, y, batch_size,signal_size,class_count,epoch,shuffle,aug=0,ncore=8):

        self.x = np.array(x,np.float32)
        self.y = y
        self.batch_size = batch_size
        self.signal_size = signal_size
        self.class_count = class_count
        self.epoch = epoch
        self.shuffle = shuffle
        self.aug = aug
        self.ncore = ncore
        if self.aug > 0:
            self.n = len(self.x) * self.aug
            self.augmented_signals, self.augmented_labels \
                = da.augment_data(self.x, self.y,self.signal_size, self.aug,self.ncore)
        else:
            self.n = len(self.x)

        print("Number of data-point: ",self.n)

        self.this_epoch = 1

    def on_epoch_end(self):

        self.this_epoch += 1
        if self.this_epoch < self.epoch:
            if self.shuffle:
                self.x,self.y = shuffle_samples(self.x,self.y)
            if self.aug > 0:
                self.augmented_signals, self.augmented_labels \
                    = da.augment_data(self.x, self.y,self.signal_size, self.aug,self.ncore)

    def __getitem__(self, index):

        start_index = index * self.batch_size
        end_index   = (index+1) * self.batch_size

        if self.aug > 0:
            batch_X = self.augmented_signals[start_index: end_index]
            batch_Y = self.augmented_labels[start_index: end_index]
        else:
            batch_X = self.x[start_index: end_index]
            batch_Y = self.y[start_index: end_index]

        batch_X = formatX(batch_X,self.signal_size)
        batch_Y = formatY(batch_Y,self.class_count)
        return batch_X,batch_Y

    def __len__(self):
        #return self.n // self.batch_size
        return int(np.ceil(self.n/self.batch_size))
