# this script is used to split to train,test reads for training.
# Prevously we extracted reads to parquet files. Now using them
# reads are separated.
# Note, during training itself train and test data can be prepared,
# as done in training/Trainning.py, splitting data into 10000,2000
# ratio for train and test. Howeever, in combined training or binary
# training, the ratio varies, so training/Trainning.py needs to be 
# modified case-by-case. Instead, in this approach we split initial
# parquet files into train,test according to out need, and then train
# by using training/Trainning_v2.py, which uses train_structuredInput

import pandas as pd
import pathlib

class split_data:
    
    def __init__(self):
        self.names_db = None

    def make_train_test(self,pathIn,pathOut,extn,
                        train_size=6000,test_size=1200):
    '''
    pathIn: This is same where you first make parquet files from fast5
            for rcc by using make_pq_rcc_ec.py or make_pq_rcc_sc.py
    pathOut: under this train and test directory will be made and output
             parquet files are placed
    train_size,test_size: number of reads in train and test dataset
    '''
        pathlib.Path(pathOut + "train/").mkdir(parents=True, exist_ok=True)
        pathlib.Path(pathOut + "test/").mkdir(parents=True, exist_ok=True)
        names = self.names_db
        for rname in sorted(names.keys()):
            gname = names[rname]
            filename = path_rcc + "/" + rname + ".pq"
            dat = pd.read_parquet(filename)
            dat_train = dat.iloc[:train_size].copy()
            dat_test  = dat.iloc[train_size:(train_size+test_size)].copy()
            train_path = pathOut + "/train/" + gname + "_" + extn + ".pq"
            test_path  = pathOut + "/test/"  + gname + "_" + extn + ".pq"
            dat_train.to_parquet(train_path)
            dat_test.to_parquet(test_path)
            print(gname)

class split_data_rcc(split_data):

    def __init__(self,rcc_name_map_file):
        self.names_db = self.get_rcc_name_map(self.rcc_name_map_file)
        
    def get_rcc_name_map(self,rcc_name_map_file):

    '''
    rcc_name_map_file: includes space-separated names of rcc isodecoders
                       at cololumn 1, the name used in sequencing data, like
                       ala1, ala2 etc. At column 2 a given name is included.
                       The given names are same as names of isodecoders used 
                       the reference, see ../resource/trna_ref.fa
                       Although the given names are same as reference sequence
                       names, for 4 special isodecoders this is not correct.
                       These are fMet, Leu1, Ile2 and Tyr. Actually, fMet rcc 
                       fraction is composed of fMet1 and fMet2, Leu1 fraction is
                       composed of Leu1 and Leu1(P), Ile2 is composed of Ile2 and 
                       Ile2v, and Tyr is composed of Tyr1 and Tyr2. Here for
                       simplicity fmet sequence data is given a name fMet, ile2 
                       is named as Ile2 (actually containd Ile2 and Ile2v), leu1
                       is named as Leu1 (actually contained Leu1 and Leu1{P}), and
                       tyr is named as Tyr.
                       For more correct name mapping see 
                       ../resource/species_trna_map
    '''
        rcc_names = {}
        with open(rcc_name_map_file,'r') as fr:
            for line in fr:
                if line[0] == '#':
                    continue
                line = line.strip()
                rname, gname = line.split(' ')
                rcc_names[rname] = gname

        return rcc_names

class split_data_ivt(split_data):
    def __init__(self,ivt_names):
        self.names_db = self.get_ivt_name_map(ivt_names)

    def get_ivt_name_map(self,ivt_names):
    '''
    ivt_names: list of ivt names in a single-column file,
               because ivt-names are directly related to reference names,
               no mapping is required, so it is single column.
               Note, here fMet1, fMet2, etc are separated.
    '''
        ivt_names = {}
        with open(ivt_names,'r') as fi:
            for line in fi:
                if line[0] == '#':
                continue
                iname = line.strip()
                ivt_names[iname] = iname
        return ivt_names

s = split_data_rcc("../resource/rcc_names") 
s.make_train_test(pathIn="/mnt/share/bhaskar/pq_db/rcc_12k",
                  pathOut="/mnt/share/bhaskar/combined_train_rcc_ivt",extn="rcc",
                  train_size=6000,test_size=1200)
s = split_data_ivt("../resource/ivt_names")
s.make_train_test(pathIn="/mnt/share/bhaskar/pq_db/ivt_8k",
                  pathOut="/mnt/share/bhaskar/combined_train_rcc_ivt",extn="ivt",
                  train_size=6000,test_size=1200)
