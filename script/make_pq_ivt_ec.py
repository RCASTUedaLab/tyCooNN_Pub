import glob
import os
import sys
import fnmatch
sys.path.append("../")
import training.GenaratePqForTrainning as pqg

class genPQ:
    LIST_FILE_LOC = "/mnt/share/bhaskar/pq_db/ivt_8k/"
    OUTPQ_DIR = "/mnt/share/bhaskar/pq_db/ivt_8k/"
    OUTSTAT_LOC   = "/mnt/share/bhaskar/pq_db/ivt_8k/"
    INP_DIR   = "/data/suzukilab/seqdata/basecall/split/"
    SETTINGS  = "../resource/makePQ.yaml"
    
    def __init__(self,settings_file = SETTINGS,takeCount = 12000,batch_num = None, 
                 outpq_dir = OUTPQ_DIR, inp_dir = INP_DIR, 
                 list_file_loc = LIST_FILE_LOC, outstat_loc = OUTSTAT_LOC):
        self.settings = settings_file
        self.outpq_dir = outpq_dir
        self.inp_dir = inp_dir
        self.takeCount = takeCount
        self.list_file_loc = list_file_loc
        self.outstat_loc = outstat_loc
        if batch_num is not None:
            self.set_batch(batch_num)
    
    def set_batch(self,batch_num):
        # fast5 file for each ivt batch are at: inp_dir + prefix +batchnum +rest_of_path,
        # Yo shold cnage accordingly.
        batch_loc = self.inp_dir + "/220820_ecoli_tRNA_IVT_batch" + str(batch_num) + "/fast5/species/"
        
        self.species = []
        with open(batch_loc + "sp.list","r") as fl: # sp.list is the list of tRNA in batch_loc
            for line in fl:
                line = line.strip()
                self.species.append(line)
        self.out_file = self.list_file_loc + "inputs" + str(batch_num) + ".txt"
        self.stat_file = self.outstat_loc + "stat" + str(batch_num) + ".txt"
        self.batch_loc = batch_loc
    
    def writeIOFile(self):

        with open(self.out_file,'w') as fw:
            for trna in self.species:
            
                lst = self.batch_loc + trna
                pqpath = self.outpq_dir + trna + ".pq"
                fw.write(trna + "\t" + lst + "\t" + pqpath + "\t" + self.stat_file + "\n")
        
    def generatePq(self):
        print("Writting from: ",self.out_file)
        print("Use settings from: ",self.settings)
        return pqg.generatePqForTrainingAll(paramPath=self.settings, 
                                            listOfIOPath=self.out_file, 
                                            takeCount=self.takeCount)

for nb in range(1,4):
    g = genPQ(takeCount=8400,batch_num=nb)
    g.writeIOFile()
    stat = g.generatePq()
