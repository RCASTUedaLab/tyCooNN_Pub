import pandas as pd
import numpy as np
import glob
import os
from multiprocessing import Pool
from functools import partial

def decide_split_train_test(N,valid=False):
    cut = (7300, 6000)
    take = None
    if N >= cut[0]:
        N_train = cut[1]
        N_test  = cut[0] - N_train 
        take = cut[0]
    else:
        N_train = np.floor(N * 0.82)
        N_test  = N - N_train
        take = N

    train_index = np.arange(N_train)
    test_index  = np.arange(N_train,take)
    if valid:
        valid_index = np.arange(N_train,N)
    else:
        valid_index = None

    return train_index,test_index,valid_index


rcc_path = "/mnt/share/bhaskar/rcc_data_extend/rcc/"
ivt_path = "/mnt/share/bhaskar/ivt_data/pq/"

files = {'rcc':sorted(glob.glob(rcc_path + "*.pq")),
         'ivt':sorted(glob.glob(ivt_path + "*.pq"))}

def split_data(tup_f,ftype,num):
    i = tup_f[0]
    f = tup_f[1]
    name,extn = os.path.splitext(os.path.basename(f))
    given_name = name + "_" + ftype
    dat = pd.read_parquet(f,engine='fastparquet')
    num_read = len(dat)
    train_index,test_index,valid_index = decide_split_train_test(num_read,valid=True)
    train = dat.iloc[train_index].copy()
    test  = dat.iloc[test_index].copy()
    valid = dat.iloc[valid_index].copy()
    l1 = len(train); l2 = len(test); l3 = len(valid);
    print("Type: %s Reading %d of %d: N=%20s Size: train: %5d test: %5d valid: %5d" % (ftype,i+1,num,given_name,l1,l2,l3))
    train.to_parquet("6k/train/" +  given_name + ".pq")
    test.to_parquet("6k/test/" + given_name + ".pq")
    valid.to_parquet("6k/valid/" + given_name + ".pq")
    return num_read

for ftype in sorted(files.keys()):
    flist = files[ftype]
    num = len(flist)
    tup_f = [(i,f) for i,f in enumerate(flist)]
    with Pool(8) as pool:
        size = pool.map(partial(split_data,ftype=ftype,num=num),tup_f)
