import glob
import os
import sys
sys.path.append("../")
import training.Trainning_gen_simple as training

def testTrain(name_db,model_loc):

    input = {'train': {}, 'test': {}}

    train_dir = "/mnt/share/bhaskar/rcc_data_extend/6k/train/"
    test_dir  = "/mnt/share/bhaskar/rcc_data_extend/6k/test/"
    names = []
    with open(name_db) as fi:
        for line in fi:
            name = line.strip()
            names.append(name)
    train_dic = {}
    test_dic  = {}
    for name in names:
        train_dic[name] = train_dir + name + ".pq"
        test_dic[name]  = test_dir  + name + ".pq"
    input['train'] = train_dic
    input['test']  = test_dic

    outdir = model_loc
    epoch = 100
    training.train_structuredInput(input, outdir, epoch,loss_fn='categorical_crossentropy')

if __name__ == '__main__':
    name_db = sys.argv[1]
    model_loc = sys.argv[2]
    testTrain(name_db,model_loc)

