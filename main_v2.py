import numpy as np
from scipy.linalg import expm
import math

####
s = 4
n_sequences = 10
####

class Sequence:
    def __init__(self,sequence):
        self.dict = {
        'A': np.array([[0,complex(0,1)],[complex(0,1),0]]),
        'T': np.array([[0,-1],[1,0]]),
        'G': np.array([[complex(0,1),0],[0,-(complex(0,1))]]),
        'C': np.array([[complex(0,1),0],[0,complex(0,1)]])
        }
        self.sequence = sequence

    def run(self):
        self.seqgroup = np.array([self.dict[s] for s in self.sequence])
        self.algebra = np.sum([i*j for i,j in zip(self.seqgroup, np.arange(0,len(self.seqgroup)))],axis=0)
        return self.algebra

def random_sequence_evolve(sequence,l_replacement=1):
    replacement = ''.join(random.choice('CGTA') for _ in range(l_replacement))
    replacement_index = random.choice(np.arange(0,len(sequence)))
    new_sequence = sequence[:replacement_index] + replacement + sequence[replacement_index + l_replacement:]
    return new_sequence

def Y(t,y):
    g1 = Sequence(sequence_array[int(t)]).run()
    gg1 = expm(g1)
    g_t = np.array([gg1[0,0], gg1[1,0]])
    return g_t


sequence_array = []
init_sequence = ''.join(random.choice('CGTA') for _ in range(7))
ss = random_sequence_evolve(init_sequence)
sequence_array.append(ss)
for t in range(T):
    ss = random_sequence_evolve(ss)
    sequence_array.append(ss)


y0 = Sequence(init_sequence).run()

for n in range(n_sequences):
    k1 = Y(y0)
    for i in range(2,s):
        
