import glob
import os
import sys
sys.path.append("../")
import training.Trainning_gen_simple as training

def testTrain(trna1,trna2,trna):

    input = {'train': {}, 'test': {}}

    train_dir = "/mnt/share/bhaskar/rcc_data_extend/rcc_split_data/train/"
    test_dir  = "/mnt/share/bhaskar/rcc_data_extend/rcc_split_data/test/"
    names = [trna1 + '_merge',trna2 + '_merge']
    train_dic = {}
    test_dic  = {}
    for name in names:
        train_dic[name] = train_dir + name + ".pq"
        test_dic[name]  = test_dir  + name + ".pq"
    input['train'] = train_dic
    input['test']  = test_dic

    outdir = "../resource/model_ivt_only_pub/" + trna + "_v2/"
    epoch = 20
    training.train_structuredInput(input, outdir, epoch, data_argument = 3,loss_fn='binary_crossentropy')

if __name__ == '__main__':
    trna1 = sys.argv[1]
    trna2 = sys.argv[2]
    trna  = sys.argv[3]
    testTrain(trna1,trna2,trna)

