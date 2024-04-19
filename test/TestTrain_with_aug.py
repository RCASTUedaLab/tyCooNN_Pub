# This script is for testing. Usually
# one should train by python TyCooNN.py train ...
# functionality

import glob
import os
import sys
sys.path.append("../")
import training.Trainning as training

def testTrain():

    '''
    Note: This training is after usual training 
          without data augmentation. It starts from
          learned model without data augmentation
          and further train.
    '''

    input = "/path/to/rcc/parquet/files"
    outdir = "/path/to/trained/model/"
    epoch = 10
    training.train(input, outdir, epoch,data_augment =3)

testTrain()

