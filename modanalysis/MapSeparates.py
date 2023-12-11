import mappy as mp
import glob

dir1 = "/mnt/share/ueda/trna_data/fastq/rcc/*.fastq"
dir2 = "/mnt/share/ueda/trna_data/fastq/ivt/*.fastq"

bam1 = "/mnt/share/ueda/trna_data/bam/rcc.bam"
bam2 = "/mnt/share/ueda/trna_data/bam/ivt.bam"
metrix1 = "/mnt/share/ueda/trna_data/metixrcc.txt"
metrix2  = "/mnt/share/ueda/trna_data/metixivt.txt"

ref = "/mnt/share/ueda/trna_data/ref/ecolitRNA_unmod_full.fa"
refdir ="/mnt/share/ueda/trna_data/ref/"

import os
import pysam



def getSQ(fasta_file):

    lengths_list = []
    current_name = None
    current_sequence = ""

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    lengths_list.append({'LN': len(current_sequence), 'SN': current_name})
                current_name = line[1:].split()[0]
                current_sequence = ""
            else:
                current_sequence += line

        if current_name:
            lengths_list.append({'LN': len(current_sequence), 'SN': current_name})

    print(lengths_list)
    return lengths_list

def getHeader(fasta_file):

    header = {'HD': {'VN': '1.0'}}
    header['SQ'] = getSQ(fasta_file)
    return header

def getlen(cigarlist):

    ll = 0
    for length,operation in cigarlist:


        if operation in [0, 1, 4, 7, 8]:
            ll += length


    return ll

def maptobam(dir,bam,ref,refdir,matrix):

    files = glob.glob(dir)
    mtx = open(matrix, 'w')

    bam_file = pysam.AlignmentFile(bam, 'wb', header=getHeader(ref))

    for file in files:
        print(file)
        filename_with_extension = os.path.basename(file)
        trna = os.path.splitext(filename_with_extension)[0]
        print(trna)
        if trna =="Ala1":
            trna = "Ala1B"

        ref = refdir+"/"+trna+".fasta"


        aligner = mp.Aligner(ref, n_threads=10,min_dp_score=15,w=4, bw=1, k=1,
                       best_n=1, min_cnt=1, min_chain_score=1)  # load or build index

        with pysam.FastxFile(file) as fh:
            cnt = 0
            mappedc = 0
            for entry in fh:
                for hit in aligner.map(entry.sequence, MD=True):

                    if hit.strand == -1:
                        continue

                    a = pysam.AlignedSegment()
                    a.query_name = entry.name
                    a.query_sequence = entry.sequence[hit.q_st:hit.q_en]
                    a.flag =  0
                    a.reference_id = bam_file.get_tid(hit.ctg)
                    a.reference_start = hit.r_st
                    a.mapping_quality = hit.mapq
                    a.cigarstring = hit.cigar_str

                    # cigar_length = getlen(hit.cigar)
                    # query_length = len(a.query_sequence)
                    # if cigar_length != query_length:
                    # print(f"Read ID: {entry.name}")
                    # print(f"CIGAR Length: {cigar_length}, Query Length: {query_length}")
                    # print(f"CIGAR String: {hit.cigar_str}")
                    # print(f"Query Sequence: {a.query_sequence}")
                        #something wrong
                    a.template_length = hit.r_en - hit.r_st
                    a.query_qualities = pysam.qualitystring_to_array(entry.quality[hit.q_st:hit.q_en])
                    # print(a.query_sequence)
                    # print(entry.quality[hit.q_st:hit.q_en])
                    a.tags = (("NM", hit.NM), ("MD", hit.MD))
                    bam_file.write(a)
                    mappedc+=1
                    break
                cnt+=1
                if cnt%100==0:
                    print(cnt,mappedc)
                if  mappedc==1000:
                    print(cnt, mappedc)
                    mtx.write(trna+"\t"+str(cnt)+"\t"+str(mappedc)+"\n")
                    break


    bam_file.close()
    mtx.close()

maptobam(dir1,bam1,ref,refdir,metrix1)
maptobam(dir2,bam2,ref,refdir,metrix2)


def split_fasta(fasta_file):
    """
    Split a FASTA file into separate files for each entry.

    :param fasta_file: Path to the input FASTA file.
    """
    with open(fasta_file, 'r') as file:
        current_file = None

        for line in file:
            if line.startswith('>'):
                if current_file:
                    current_file.close()

                # Extract the sequence ID from the header line
                seq_id = line.split()[0][1:]  # Remove '>' and split by whitespace
                filename = f"{seq_id}.fasta"  # Create a file name based on the sequence ID

                # Open a new file for writing
                current_file = open("/mnt/share/ueda/trna_data/ref/"+filename, 'w')

            if current_file:
                current_file.write(line)

        if current_file:
            current_file.close()

# gp
split_fasta(ref)