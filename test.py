from Bio import pairwise2
import numpy as np
import random
import threading
from core_debugging import *

y0 = np.array([[complex(1,0),complex(0,0)],[complex(0,0),complex(-1,0)]])

def DNA(l):
    return ''.join([random.choice('CGTA') for i in range(l)])

def pairwise_alignment(s1,s2):
    score = pairwise2.align.localxx(s1, s2, score_only=True, one_alignment_only=True)
    return score

def per_thread(start,stop):
    n_sequences = 5
    sequences = [DNA(30) for i in range(n_sequences)]
    distances = np.zeros((n_sequences,n_sequences))


    n_sequences = len(sequences)
    distances = np.zeros((n_sequences,n_sequences))

    for i in range(n_sequences):
        for j in range(i+1,n_sequences):
            align_score = pairwise_alignment(sequences[i],sequences[j])
            distances[i,j] = align_score
            distances[j,i] = align_score


    ordered_sequences = []
    idx = random.choice(np.arange(n_sequences))
    distances[:,idx] = 0
    ordered_sequences.append(sequences[idx])
    while len(ordered_sequences)<n_sequences:
        idx = np.argmax(distances[idx])
        distances[:,idx] = 0
        ordered_sequences.append(sequences[idx])

    def Y(y,t):
        seq = ordered_sequences[t]
        alg = Sequence(seq).run()
        y_t = matrix_multiply(y,expm(alg))
        return y_t

    print(ordered_sequences)
    y_array = []
    y = y0
    for n in range(n_sequences):
        y = rkmk_step(Y,y,n)
        y_array.append(y)
    return y_array

out = per_thread(1,2)
for i in out:
    print(i)
