from Bio import SeqIO
from core import *
import glob
import numpy as np
import random

indir = "./Data/hmm_sequences/"
outdir = "./Data/background_yarrays"

y0 = np.array([[complex(0,1),complex(0,0)],[complex(0,0),complex(0,0)]])

def DNA(l):
    return ''.join([random.choice('CGTA') for i in range(l)])

files = glob.glob(f"{indir}/*")
lengths = []

for record in SeqIO.parse(files[0],"fasta"):
    lengths.append(len(record.seq))
n_sequences = len(lengths)

for ff in range(1000):
    print(f"{ff} of 1000")

    def Y(y,t):
        seq = DNA(lengths[t])
        y_t = matrix_multiply(y,expm(Sequence(seq).run()))
        return y_t

    y_array = []
    y = Y(y0,0)
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array.append(y)

    print(f"{file_num} Writing output..")
    with open(f"{outdir}/f{file_num}_backgroundyarray.pkl",'wb') as outfile:
        pkl.dump(y_array,outfile)
