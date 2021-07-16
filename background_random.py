from Bio import SeqIO
from core import *
import numpy as np
import random
import pickle as pkl

indir = "./Data/hmm_sequences/"
outdir = "./Data/background_yarrays"

y0 = np.array([[complex(0.5,0.5),complex(1,1)],[complex(0.5,0.5),complex(1,1)]])

def DNA(l):
    return ''.join([random.choice('CGTA') for i in range(l)])


lengths = []
for record in SeqIO.parse(f"./Data/hmm_sequences/1.fasta","fasta"):
    lengths.append(len(record.seq))
n_sequences = len(lengths)

all_y = np.zeros((1000,20,2,2),dtype="complex128")

for i in range(1,1001):
    print(f"{i} of 1000")

    def Y(y,t):
        seq = DNA(lengths[t])
        alg = Sequence(seq).run()
        y_t = matrix_multiply(y,expm(alg))
        return y_t

    y_array = np.zeros((20,2,2),dtype="complex128")
    y = y0
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array[n] = y

    all_y[i] = y_array


print(f"{ff} Writing output..")
with open(f"{outdir}/background_all_y.pkl",'wb') as outfile:
    pkl.dump(all_y,outfile)
