# This script is for testing. Usually 
# one should train by python TyCooNN.py train ...
# functionality

import glob
import os
import sys
sys.path.append("../")
import training.Trainning as training

def testTrain():

    input = "/path/to/rcc/parquet/files"
    outdir = "/path/to/trained/model/"
    epoch = 100
    training.train(input, outdir, epoch)

testTrain()

