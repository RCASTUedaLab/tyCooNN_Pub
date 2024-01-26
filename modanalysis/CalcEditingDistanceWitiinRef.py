from Bio import SeqIO
import Levenshtein

def read_fasta(file_path):
    """FASTAt@Cz"""
    with open(file_path, "r") as file:
        return [record.seq for record in SeqIO.parse(file, "fasta")]

def find_closest_sequence(sequences, index):
    """wCfbNXzz"""
    min_distance = float('inf')
    closest_seq = None
    target_seq = sequences[index]

    for i, seq in enumerate(sequences):
        if i != index:
            distance = Levenshtein.distance(str(target_seq), str(seq))
            if distance < min_distance:
                min_distance = distance
                closest_seq = seq

    return closest_seq, min_distance

# FASTAt@C
# sequences = read_fasta("/share/reference/trna/ecolitRNA_unmod.fa")
sequences = read_fasta("/share/reference/trna/ecolitRNA_mod.fa")

# ezz
for index in range(len(sequences)):
    closest_seq, distance = find_closest_sequence(sequences, index)
    print(f"Distance: {distance}")
