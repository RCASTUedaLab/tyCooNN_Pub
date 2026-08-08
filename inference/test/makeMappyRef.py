from Bio import SeqIO
import os

def split_fasta(fasta_file):
    directory = os.path.dirname(fasta_file)
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):

            record.seq = record.seq.upper().replace('U', 'T')
            file_name = f"{record.id}.fa"
            output_path = os.path.join(directory, file_name)
            with open(output_path, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")


fasta_file = "/mnt/share/ueda/trna_data/ref/ecolitRNA_unmod_full.fa"
split_fasta(fasta_file)