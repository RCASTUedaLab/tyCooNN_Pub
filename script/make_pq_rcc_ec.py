import glob
import os
import sys
import fnmatch
sys.path.append("../")
import preprocess.MakeTrainingPq as mkpq

class genPQ:

    '''
    LIST_FILE = "path_to_fast5_list_file/inputs.txt"  # Output file
    OUTPQ_DIR = "path_to_parquet_directory"           # Output location 
    OUTSTAT   = "path_to_fast5_list_file/stat.txt"    # Preprocessing statistics, output
    INP_DIR   = "path_to_fast5_basecalled/"           # input
    SETTINGS  = "../resource/settings.yaml"             # Preprocessing options, input
    '''
    
    LIST_FILE = "/mnt/share/bhaskar/pq_db/rcc_12k/ec/inputs.txt"
    OUTPQ_DIR = "/mnt/share/bhaskar/pq_db/rcc_12k/ec/"
    OUTSTAT   = "/mnt/share/bhaskar/pq_db/rcc_12k/ec/stat.txt"
    INP_DIR   = "/mnt/share/suzukilab/rcc/"
    SETTINGS  = "../resource/settings.yaml"
    
    def __init__(self,settings_file = SETTINGS):
        self.settings = settings_file

    def writeIOFile(self,out_file=LIST_FILE,outpq_dir=OUTPQ_DIR,stat_file=OUTSTAT,inp_dir=INP_DIR):

        fw = open(out_file,'w')

        self.out_file = out_file
        self.outpq_dir = outpq_dir
        self.outstat = stat_file
        self.inp_dir = inp_dir

        files = sorted(glob.glob(inp_dir+"/*/"))

        for f in files:
            dirlist = {}
            for root, dirnames,filenames in os.walk(f):
                for filename in fnmatch.filter(filenames,"*.fast5"):
                    dirlist[root] = 1
            dirlist = list(dirlist.keys())
            print(dirlist)

            # Change accordingly the next two lines
            basename = os.path.basename(os.path.dirname(f))

            pqpath = outpq_dir + basename + ".pq"
            statpath = stat_file
            lst = ",".join(dirlist)

            fw.write(basename + "\t" + lst + "\t" + pqpath + "\t" + statpath + "\n")

        fw.close()

    def generatePq(self,takeCount=12000):
        print("Writting from: ",self.out_file)
        print("Use settings from: ",self.settings)
        return mkpq.generatePqForTrainingAll(paramPath=self.settings, listOfIOPath=self.out_file, takeCount=takeCount)

g = genPQ()
g.writeIOFile()
# 12000 reads taken per tRNA
#stat = g.generatePq(takeCount=12000)
