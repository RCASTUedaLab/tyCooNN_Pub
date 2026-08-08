import pandas as pd
from tensorflow import keras
import nnmodels.CNNWavenet as cnnwavenet
from numba import jit
import multiprocessing
import csv
import numpy as np
from tensorflow.keras.callbacks import ModelCheckpoint
import utils.tyUtils as ut
from multiprocessing import Pool
import psutil,os

def augment_data(signals, labels, signal_size,augmentation_factor,ncore=8):
    if augmentation_factor <= 1:
        print('Not performing data augmentation')
        signals,labels = suffle(signals,labels)
        return signals, labels
    print('start data argumentaion')
    n_cores = ncore
    #n_cores = multiprocessing.cpu_count() // 5
    print('n_cores', n_cores)
    pool = Pool(n_cores)
    o_results = []
    results = []
    j = 0
    for signal, label in zip(signals, labels):
        o_results.append((signal, label))
        for _ in range(augmentation_factor - 1):
            results.append(pool.apply_async(modify_signal_l, args=(signal, label)))
            j += 1
    pool.close()
    pool.join()
    results = [p.get() for p in results]
    results.extend(o_results)
    print('end data argumentaion')
    augmented_signals, augmented_labels = margeAndSuffle(results,signal_size)
    print('end data merge')
    return augmented_signals, augmented_labels

#@jit(nopython=True)
def margeAndSuffle(results,signal_size):

    datasize = len(results)
    shuffle_indices = np.random.permutation(np.arange(datasize))
    # print("si",shuffle_indices)
    x = np.empty([datasize, signal_size], dtype=np.float32)
    y = np.empty(datasize, dtype=int)
    augmented_signals,augmented_labels = _merge(x, y, results, shuffle_indices)
    # print("al",augmented_labels)
    return augmented_signals,augmented_labels

def suffle(signals,labels):

    datasize = len(signals)
    labels = np.array(list(labels))
    print(labels)
    shuffle_indices = np.random.permutation(np.arange(datasize))
    sufflex = signals[shuffle_indices]
    suffley = labels[shuffle_indices]
    return sufflex,suffley

#@jit(nopython=True)
def _merge(x,y,results,shuffle_indices):

    for idx in range(len(shuffle_indices)):
        i = shuffle_indices[idx]
        x[idx] = results[i][0]
        y[idx] = results[i][1]
    return x,y


@jit
def modify_signal_l(signal, label):
    sig = modify_signal(signal)
    return (sig, label)

@jit(nopython=True)
def modify_signal(signal):

    modification_count = len(signal) * 0.3
    modification_count = int(round(modification_count / 2)) * 2 # to int
    #modification_positions = random.sample(range(len(signal)), k=modification_count)
    modification_positions = np.random.choice(len(signal),modification_count)

    half = int(modification_count / 2)
    duplication_positions = set(modification_positions[:half])
    deletion_positions = set(modification_positions[half:])

    new_signal = np.empty(len(signal), np.float32)
    j = 0
    for i, val in enumerate(signal):
        if i in duplication_positions and j < len(signal):

            new_signal[j] = val
            j += 1
            new_signal[j] = val
            j += 1
        elif i in deletion_positions:
            pass
        elif j < len(signal):
            new_signal[j] = val
            j += 1

    return new_signal

