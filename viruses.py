from core import *
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
from scipy.linalg import expm


data_dir = "./Data"
fig_dir = "./Figures"

def Y(y,t):
    y_t = Sequence(sequence_array[t]).run()
    return expm(y_t)

# Load data
virus_genomic_sequences = []
virus_names = []
for record in SeqIO.parse(f"{data_dir}/viral.3.1.genomic.fna", "fasta"):
    virus_names.append(''.join([i[1:] for i in record.description.split(" ")]))
    virus_genomic_sequences.append(record.seq.upper())

n_sequences = len(virus_genomic_sequences)
y_array = []

for n in n_sequences:
    seq = virus_genomic_sequences[n]
    y = expm(Sequence(seq).run())
    next_y = rkmk_step(Y,y)
    y_array.append(next_y)



    y_array.append(y)
