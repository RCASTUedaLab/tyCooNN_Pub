import glob
import os
import sys
import fnmatch
sys.path.append("../")
import training.GenaratePqForTrainning as pqg

class genPQ:

    LIST_FILE = "/mnt/share/bhaskar/pq_db/rcc_12k/yeast/inputs.txt"
    OUTPQ_DIR = "/mnt/share/bhaskar/pq_db/rcc_12k/yeast/"
    OUTSTAT   = "/mnt/share/bhaskar/pq_db/rcc_12k/yeast/stat.txt"
    INP_DIR   = "/data/share/spikein_fast5/basecalled"
    SETTINGS  = "../resource/makePQ.yaml"

    def __init__(self,settings_file = SETTINGS):
        self.settings = settings_file
        self.pqNames = {'phe_mature':'spikeMat','phe_i1':'spikeInt'}

    def writeIOFile(self,out_file=LIST_FILE,outpq_dir=OUTPQ_DIR,stat_file=OUTSTAT,inp_dir=INP_DIR):

        fw = open(out_file,'w')

        self.out_file = out_file
        self.outpq_dir = outpq_dir
        self.outstat = stat_file
        self.inp_dir = inp_dir

        files = glob.glob(inp_dir+"/*/")

        for f in files:
            dirlist = {}
            for root, dirnames,filenames in os.walk(f):
                for filename in fnmatch.filter(filenames,"*.fast5"):
                    dirlist[root] = 1
            dirlist = list(dirlist.keys())

            # Change accordingly the next line
            basename = os.path.basename(os.path.dirname(f)).replace("200710_yeast_tx_","").replace("_V7","")

            basename = self.pqNames[basename]

            pqpath = outpq_dir + basename + ".pq"
            statpath = stat_file
            lst = ",".join(dirlist)

            fw.write(basename + "\t" + lst + "\t" + pqpath + "\t" + statpath + "\n")

        fw.close()

    def generatePq(self,takeCount=12000):
        print("Writting from: ",self.out_file)
        print("Use settings from: ",self.settings)
        return pqg.generatePqForTrainingAll(paramPath=self.settings, listOfIOPath=self.out_file, takeCount=takeCount)

g = genPQ()
g.writeIOFile()
stat = g.generatePq(takeCount=12000)
