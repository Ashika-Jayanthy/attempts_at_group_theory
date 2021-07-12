from Bio import SeqIO
from core import *
import numpy as np
import random
import pickle as pkl

indir = "./Data/hmm_sequences/"
outdir = "./Data/background_yarrays"

y0 = np.array([[complex(1,0),complex(0,0)],[complex(0,0),complex(0,0)]])

def DNA(l):
    return ''.join([random.choice('CGTA') for i in range(l)])


lengths = []
for record in SeqIO.parse(f"./Data/hmm_sequences/1.fasta","fasta"):
    lengths.append(len(record.seq))
n_sequences = len(lengths)

for ff in range(1000):
    print(f"{ff} of 1000")

    def Y(y,t):
        seq = DNA(lengths[t])
        alg = Sequence(seq).run()
        y_t = matrix_multiply(y,expm(alg))
        return y_t

    y_array = []
    y = y0
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array.append(y)

    print(f"{ff} Writing output..")
    with open(f"{outdir}/f{ff}_backgroundyarray.pkl",'wb') as outfile:
        pkl.dump(y_array,outfile)
