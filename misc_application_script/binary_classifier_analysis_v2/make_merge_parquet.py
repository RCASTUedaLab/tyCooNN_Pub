import pandas as pd
import numpy as np
import glob
import os,sys

trna1 = sys.argv[1]
trna2 = sys.argv[2]
trna  = sys.argv[3]
sizes = [int(s) for s in sys.argv[4].split(",")]

class makeData:

    def __init__(self,trna,trna1,trna2):
        self.trna = trna
        self.trna1 = trna1
        self.trna2 = trna2

    def get_data(self,trnas,isodecoders,path,extn_suffix=".pq"):
        filepaths = [path + "/" + t for t in trnas]
        trnas = [trna1,trna2]
        files = [filepaths[i] + "/" + t + extn_suffix for i,t in enumerate(isodecoders)]
        dat = [pd.read_parquet(file) for file in files]
        return dat

    def set_con(self,dat):
        self.dat_con = dat
        for i in range(len(self.dat_con)):
            print("Number of concordant reads [%d]: %d " % (i,len(self.dat_con[i])))

    def set_ivt(self,dat):
        self.dat_ivt = dat
        for i in range(len(self.dat_ivt)):
            print("Number of ivt readst reads [%d]: %d " % (i,len(self.dat_ivt[i])))

    def make_merge(self,i,r1,r2):
        d1 = self.dat_ivt[i]
        d2 = self.dat_con[i]
        t1 = d1.iloc[r1].copy()
        t2 = d2.iloc[r2].copy()
        return pd.concat([t1,t2])

    def save_to_parquet(self,dat,file):
        dat.to_parquet(file)
        return

d = makeData(trna,trna1,trna2)
isodecoder = (trna1,trna2)
trna_list  = (trna,trna)
d.set_con(d.get_data(trna_list,isodecoder,
                     path="/mnt/share/bhaskar/rcc_data_extend/rcc_split_data",
                     extn_suffix="_con12.pq"))
trna_list = (trna1,trna2)
d.set_ivt(d.get_data(trna_list,isodecoder,
                path="/mnt/share/bhaskar/ivt_data"))

N_train = sizes[0]
N_map1_train = sizes[1]; N_ivt1_train = N_train - N_map1_train;
N_map2_train = sizes[2]; N_ivt2_train = N_train - N_map2_train;

print(trna, trna1, trna2)
print(" TRAIN: ")
print(trna1," RCC(by IVT NN):",N_ivt1_train," RCC(by map):",N_map1_train)
print(trna2," RCC(by IVT NN):",N_ivt2_train," RCC(by map):",N_map2_train)
train_merge1 = d.make_merge(0,np.arange(N_ivt1_train),np.arange(N_map1_train))
train_merge2 = d.make_merge(1,np.arange(N_ivt2_train),np.arange(N_map2_train))

N_test = sizes[3]
N_map1_test = sizes[4]; N_ivt1_test = N_test - N_map1_test;
N_map2_test = sizes[5]; N_ivt2_test = N_test - N_map2_test;

print(" TEST: ")
print(trna1," RCC(by IVT NN):",N_ivt1_test," RCC(by map):",N_map1_test)
print(trna2," RCC(by IVT NN):",N_ivt2_test," RCC(by map):",N_map2_test)
test_merge1 = d.make_merge(0,np.arange(N_ivt1_train,N_ivt1_train+N_ivt1_test),np.arange(N_map1_train,N_map1_train+N_map1_test))
test_merge2 = d.make_merge(1,np.arange(N_ivt2_train,N_ivt2_train+N_ivt2_test),np.arange(N_map2_train,N_map2_train+N_map2_test))

d.save_to_parquet(train_merge1,"train/" + isodecoder[0] + "_merge.pq")
d.save_to_parquet(train_merge2,"train/" + isodecoder[1] + "_merge.pq")

d.save_to_parquet(test_merge1,"test/" + isodecoder[0] + "_merge.pq")
d.save_to_parquet(test_merge2,"test/" + isodecoder[1] + "_merge.pq")
