from Bio import pairwise2
from Bio import SeqIO
from core import *
import glob
import sys
import numpy as np
import random

indir = "./Data/hmm_sequences/"
outdir = "./Data/yarrays"

def pairwise_alignment(s1,s2):
    score = pairwise2.align.globalxx(s1, s2, score_only=True, one_alignment_only=True)
    return score

files = glob.glob(f"{indir}/*")

ordered_sequences = []

for ff,file in enumerate(files):
    sequences = []
    for record in SeqIO.parse(file,"fasta"):
        sequences.append(record.seq)
    n = len(sequences)
    distances = np.zeros((n,n))

    for i in range(n):
        for j in range(i+1,n):
            align_score = pairwise_alignment(sequences[i],sequences[j])
            distances[i,j] = align_score
            distances[j,i] = align_score

    idx = random.choice(np.arange(n))
    distances[:,idx] = 0
    ordered_sequences.append(sequences[idx])
    while len(ordered_sequences)<len(sequences)-1:
        new_idx = np.argmax(distances[idx])
        distances[:,new_idx] = 0
        ordered_sequences.append(sequences[new_idx])
        idx = new_idx



    def Y(t):
        seq = ordered_sequences[t]
        y_t = expm(Sequence(seq).run())
        return y_t


    n_sequences = len(ordered_sequences)
    y_array = []

    for n in range(n_sequences):
        y = Y(n)
        next_y = rkmk_step(Y,y,n)
        y_array.append(next_y)
        print(f"{n} of {n_sequences}")

    for i in y_array:
        np.savetxt(f"{outdir}/f{ff}_{i}_yarray.txt", y_array[i], fmt='%.18e', delimiter=' ', newline='\n')
