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

ivt_species = {}
with open("ivt_data/species_trna_map") as fs:
    for line in fs:
        line = line.strip()
        names = line.split(" ")
        ivt_name = names[0]
        given_name = names[0]
        ivt_species[ivt_name] = given_name

inp_dir = "/data/suzukilab/seqdata/basecall/split/"
batch = ["220820_ecoli_tRNA_IVT_batch1","220820_ecoli_tRNA_IVT_batch2","220820_ecoli_tRNA_IVT_batch3"]
for b in batch:
    d = inp_dir + b + "/fast5/species/"
    with open(d + "sp.list") as fl:
        for line in fl:
            line = line.strip()
            ivt_species[line] = d + line

for i,given_name in enumerate(sorted(ivt_species.keys())):
    path = ivt_species[given_name]
    print(given_name,path)
    fs = glob.glob(path + "/*.fast5")
    Path("ivt_data/" + given_name).mkdir(parents=True, exist_ok=True)
    print(i,given_name)
    fcount = 0
    for f in fs:
        fname = "%04d" % fcount
        linked_file = "ivt_data/" + given_name + "/" + fname + ".fast5"
        source_file = f
        symlink_force(source_file,linked_file)
        print("\t%s" % linked_file)
        fcount += 1

