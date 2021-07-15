from scipy.linalg import expm
import numpy as np
import random
import collections
from Bio import pairwise2

####
s = 4
A = np.array(
    [
        [0,   0,   0,   0],
        [0.5, 0,   0,   0],
        [0,   0.5, 0,   0],
        [0,   0,   1.0, 0]
    ]
)
b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])
c = np.array([0, 0.5, 0.5, 1.0])

m1, m2, m3 = 2, 2, -1


generators = np.array([np.array([[complex(0.,0.),complex(0.,1.)],[complex(0.,1.),complex(0.,0.)]]),
np.array([[complex(0.,0.),complex(1.,0.)],[complex(-1.,0.),complex(0.,0.)]]),
np.array([[complex(0.,1.),complex(0.,0.)],[complex(0.,0.),complex(0.,-1.)]]),
np.array([[complex(0.,1.),complex(0.,0.)],[complex(0.,0.),complex(0.,1.)]])])
####


def pairwise_alignment(s1,s2):
    score = pairwise2.align.localxx(s1, s2, score_only=True, one_alignment_only=True)
    return score

def DNA(l):
    return ''.join(random.choice('CGTA') for _ in range(l))

def matrix_multiply(a,b):
    n1,m1 = a.shape
    n2,m2 = b.shape
    ans = np.zeros((n1,m2),dtype="complex128")
    for row in range(n1):
        row_value = a[row]
        for column in range(m2):
            column_value = b.T[column]
            ans[row,column] = np.vdot(row_value,column_value)
    return ans

def commutator(a,b):
    return matrix_multiply(a,b) - matrix_multiply(b,a)

def condition_check(val, type="matrix"):
    if type == "matrix":
        a,b,c,d = np.real(val[0,0]), np.imag(val[0,0]), np.real(val[1,0]), np.imag(val[1,0])
    elif type == "vector":
        a,b,c,d = np.real(val[0]), np.imag(val[0]), np.real(val[1]), np.imag(val[1])
    return a*a + b*b + c*c + d*d

def rkmk_step(Y,y,n,h=1e-7):
    k = np.zeros((s,2,2), dtype="complex128")
    I1 = Y(y,n)
    k[0] = Y(y,n)
    for i in range(1,s):
        u = h * np.sum([A[i,j] * k[j] for j in range(i)], axis=0)
        u_tilda = u + (((c[i] * h) / 6) * commutator(I1, u))
        k[i] = Y(matrix_multiply(y, expm(u_tilda)), n)

    I2 = ((m1 * (k[1] - I1)) + (m2 * (k[2] - I1)) + (m3 * (k[3] - I1))) / h
    v = h * np.sum([b[j] * k[j] for j in range(s)], axis=0)
    v_tilda = v + ((h / 4) * commutator(I1,v)) + ((h**2 / 24) * commutator(I2,v))
    y = matrix_multiply(y, expm(v_tilda))
    print(condition_check(y))
    if not np.isclose(condition_check(y),1.):
        return 'NaN'
    else:
        return y

def random_sequence_evolve(sequence,l_replacement=1):
    replacement = ''.join(random.choice('CGTA') for _ in range(l_replacement))
    replacement_index = random.choice(np.arange(0,len(sequence)))
    new_sequence = sequence[:replacement_index] + replacement + sequence[replacement_index + l_replacement:]
    return new_sequence


class Sequence:
    def __init__(self,sequence):
        self.dict = {
        'A': np.array([[complex(0.,0.),complex(0.,1.)],[complex(0.,1.),complex(0.,0.)]]),
        'T': np.array([[complex(0.,0.),complex(1.,0.)],[complex(-1.,0.),complex(0.,0.)]]),
        'G': np.array([[complex(0.,1.),complex(0.,0.)],[complex(0.,0.),complex(0.,-1.)]]),
        'C': np.array([[complex(0.,1.),complex(0.,0.)],[complex(0.,0.),complex(0.,1.)]]),
        }
        self.sequence = sequence

    def run(self):
        self.seqgroup = np.array([self.dict[s] if s in self.dict.keys() else random.choice(generators) for s in self.sequence]) #check source of non-N read errors
        self.algebra = np.sum([i*j for i,j in zip(self.seqgroup, np.arange(0,len(self.seqgroup)))],axis=0)
        return self.algebra
