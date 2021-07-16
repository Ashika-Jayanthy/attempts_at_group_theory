from core import *
import numpy as np
import pickle as pkl


indir = "./Data/ordered_sequences"
outdir = "./Data/yarrays"
y0 = np.array([[complex(0.5,0.5),complex(1,1)],[complex(0.5,0.5),complex(1,1)]])

all_y = np.zeros((1000,20,2,2),dtype="complex128")

for i in range(1,1001):
    print(f"{i} of 1000")

    with open(f"{indir}/f{i}_ordered_sequence.pkl",'rb') as michaelscott:
        ordered_sequences = pkl.load(michaelscott)

    n_sequences = len(ordered_sequences)

    def Y(y,t):
        seq = ordered_sequences[t]
        alg = Sequence(seq).run()
        y_t = matrix_multiply(y,expm(alg))
        return y_t

    y_array = np.zeros((20,2,2),dtype="complex128")
    y = y0
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array[n] = y

    all_y[i-1] = y_array

print(f"{i} Writing output..")
with open(f"{outdir}/virus_all_y.pkl",'wb') as outfile:
    pkl.dump(all_y,outfile)
