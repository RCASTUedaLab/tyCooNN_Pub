import glob
import os
import training.GenaratePqForTrainning as pqg
import sys
import training.Trainning as traning

def testTrain():

    input = "/share/trna/tyCooNNTest/trim12000/"
    outdir = "/share/trna/tyCooNNTest/testout"
    epoch = 100
    traning.train(input, outdir, epoch)
    epoch = 10
    #traning.train(input, outdir, epoch,data_argument =3)

testTrain()

