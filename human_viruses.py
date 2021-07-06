from Bio import pairwise2
from Bio import SeqIO
from core import *
import glob
import sys

indir = "./Data/hmm_sequences/"

def pairwise_alignment(s1,s2):
    score = pairwise2.align.globalxx(s1, s2, score_only=True, one_alignment_only=True)
    return score

files = glob.glob(f"{indir}/*")


for file in files:
    sequences = []
    for record in SeqIO.parse(file,"fasta"):
        sequences.append(record.seq)
    distances
    for i in range(len(sequences)):
        for j in range(i+1,len(sequences)):
