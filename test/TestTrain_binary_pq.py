import glob
import os
import sys
sys.path.append("../")
import training.Trainning_binary as training

def testTrain(trna):

    #input = {'train': {'wt': '/mnt/share/bhaskar/pq_binary/WT_TruB/train/WT/' + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TruB/train/KO/' + trna + '.pq'},
    #         'test':  {'wt': '/mnt/share/bhaskar/pq_binary/WT_TruB/test/WT/'  + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TruB/test/KO/'  + trna + '.pq'}}
    #outdir = "/mnt/share/bhaskar/pq_binary/WT_TruB/train/model/" + trna
    #input = {'train': {'wt': '/mnt/share/bhaskar/pq_binary/WT_TrmB/train/WT/' + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TrmB/train/KO/' + trna + '.pq'},
    #         'test':  {'wt': '/mnt/share/bhaskar/pq_binary/WT_TrmB/test/WT/'  + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TrmB/test/KO/'  + trna + '.pq'}}
    #outdir = "/mnt/share/bhaskar/pq_binary/WT_TrmB/train/model/" + trna
    #input = {'train': {'wt': '/mnt/share/bhaskar/pq_binary/WT_TapT/train/WT/' + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TapT/train/KO/' + trna + '.pq'},
    #         'test':  {'wt': '/mnt/share/bhaskar/pq_binary/WT_TapT/test/WT/'  + trna + '.pq',
    #                   'ko': '/mnt/share/bhaskar/pq_binary/WT_TapT/test/KO/'  + trna + '.pq'}}
    #outdir = "/mnt/share/bhaskar/pq_binary/WT_TapT/train/model/" + trna
    input = {'train': {'wt': '/mnt/share/bhaskar/pq_binary/WT_TapT/train/WT_1/' + trna + '.pq',
                       'ko': '/mnt/share/bhaskar/pq_binary/WT_TapT/train/KO_1/' + trna + '.pq'},
             'test':  {'wt': '/mnt/share/bhaskar/pq_binary/WT_TapT/test/WT_1/'  + trna + '.pq',
                       'ko': '/mnt/share/bhaskar/pq_binary/WT_TapT/test/KO_1/'  + trna + '.pq'}}
    outdir = "/mnt/share/bhaskar/pq_binary/WT_TapT/train/model_1/" + trna
    epoch = 200
    training.train_structuredInput(input, outdir, epoch)
    epoch = 10
    training.train_structuredInput(input, outdir, epoch,data_argument =3)
    #epoch = 25
    #training.train_structuredInput(input, outdir, epoch,data_argument =5)

trna = sys.argv[1]
testTrain(trna)

