import glob
import os,errno
from pathlib import Path

def symlink_force(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)
        else:
            raise e

rcc_species = {}
with open("combined_train_rcc_ivt/species_trna_map") as fs:
    for line in fs:
        line = line.strip()
        names = line.split(" ")
        rcc_name = names[1]
        given_name = names[2]
        rcc_species[rcc_name] = given_name

inp_dir = "/data/share/spikein_fast5/basecalled"
files = glob.glob(inp_dir + "/*/")

for i,f in enumerate(files):
    dirlist = []
    fs = glob.glob(f+"/*/workspace")
    if len(fs)>0:
        dirlist = fs
    else:
        fs = glob.glob(f + "/*/*/*/workspace")
        dirlist = fs
    basename = os.path.basename(os.path.dirname(f)).replace("ecoli_rcc_","")
    basename = basename.split("_")[3] + "_" + basename.split("_")[4]
    if basename == 'phe_i1': 
        basename = 'pheNoMat'
    elif basename == 'phe_mature':
        basename = 'pheMat'
    else:
        continue
    given_name = rcc_species[basename]
    Path("rcc_data/" + given_name).mkdir(parents=True, exist_ok=True)
    print(i,basename, given_name)
    fcount = 0
    for d in dirlist:
        fast5 = glob.glob(d + "/*.fast5")
        print("\t",d)
        for f5 in fast5:
            fname = "%04d" % fcount
            linked_file = "rcc_data/" + given_name + "/" + fname + ".fast5"
            source_file = f5
            symlink_force(source_file,linked_file)
            print("%s" % linked_file)
            fcount += 1

