from Bio import pairwise2
from Bio import SeqIO
from core import *
import glob
import sys
import numpy as np
import random

indir = "./Data/hmm_sequences/"

def pairwise_alignment(s1,s2):
    score = pairwise2.align.globalxx(s1, s2, score_only=True, one_alignment_only=True)
    return score

files = glob.glob(f"{indir}/*")

ordered_sequences = []

for file in files[0:10]:
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
    print(idx)
    indexes = []
    indexes.append(idx)
    ordered_sequences.append(sequences[idx])
    while len(ordered_sequences)<len(sequences):
        idx = np.argmax(distances[np.mask(distances[idx],indexes)])
        indexes.append(idx)
        print(idx)
        ordered_sequences.append(sequences[idx])
