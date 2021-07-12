from Bio import pairwise2
from Bio import SeqIO
from core import *
import numpy as np
import random
import threading
import pickle as pkl


indir = "./Data/ordered_sequences"
y0 = np.array([[complex(1,0),complex(0,0)],[complex(0,0),complex(0,0)]])


for i in range(1,1001):
    with open(f"{indir}/f{i}_ordered_sequence.pkl",'rb') as michaelscott:
        ordered_sequences = pkl.load(michaelscott)

    n_sequences = len(ordered_sequences)

    def Y(y,t):
        seq = ordered_sequences[t]
        alg = Sequence(seq).run()
        y_t = matrix_multiply(y,expm(alg))
        return y_t

    y_array = []
    y = y0
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array.append(y)

    print(y_array)

    #print(f"{file_num} Writing output..")
    #with open(f"{outdir}/f{file_num}_yarray.pkl",'wb') as outfile:
        #pkl.dump(y_array,outfile)
