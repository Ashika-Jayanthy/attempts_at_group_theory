from Bio import SeqIO
from core import *
import glob
import numpy as np
import random

indir = "./Data/hmm_sequences/"
outdir = "./Data/background_yarrays"



def DNA(l):
    return ''.join([random.choice('CGTA')])

files = glob.glob(f"{indir}/*")
lengths = []

for record in SeqIO.parse(files[0],"fasta"):
    lengths.append(len(record.seq))
n_sequences = len(lengths)

for ff in range(1000):
    print(f"{ff} of 1000")

    def Y(t):
        seq = DNA(lengths[t])
        y_t = expm(Sequence(seq).run())
        return y_t

    y_array = []

    for n in range(n_sequences):
        y = Y(n)
        next_y = rkmk_step(Y,y,n)
        y_array.append(next_y)


    for i in range(len(y_array)):
        np.savetxt(f"{outdir}/f{ff}_{i}_yarray.txt", y_array[i], fmt='%.18e', delimiter=' ', newline='\n')
