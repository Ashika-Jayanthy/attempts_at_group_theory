from core import *
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
from scipy.linalg import expm

# Data: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

data_dir = "./Data"
fig_dir = "./Figures"

def Y(y,t):
    y_t = Sequence(sequence_array[t]).run()
    return expm(y_t)

# Load data
virus_genomic_sequences = []
virus_names = []
for record in SeqIO.parse(f"{data_dir}/viral.3.1.genomic.fna", "fasta"):
    virus_names.append(' '.join([i for i in record.description.split(" ")[1:]]))
    virus_genomic_sequences.append(record.seq.upper())

print(virus_names)

n_sequences = len(virus_genomic_sequences)
y_array = []

for n in range(n_sequences):
    seq = virus_genomic_sequences[n]
    y = expm(Sequence(seq).run())
    next_y = rkmk_step(Y,y,n)
    y_array.append(next_y)
