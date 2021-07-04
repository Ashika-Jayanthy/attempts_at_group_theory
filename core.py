from scipy.linalg import expm
import numpy as np
import random


####
s = 4
A = np.array(
    [
        [0,   0,   0,   0],
        [0.5, 0,   0,   0],
        [0,   0.5, 0,   0],
        [0,   0,   1.0, 0]
    ], dtype="complex128"
)
b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6], dtype="complex128")
c = np.array([0, 0.5, 0.5, 1.0])

m1, m2, m3 = 2, 2, -1

generators = np.array([np.array([[0,complex(0,1)],[complex(0,1),0]]),
np.array([[0,-1],[1,0]]),
np.array([[complex(0,1),0],[0,-(complex(0,1))]]),
np.array([[complex(0,1),0],[0,complex(0,1)]])])
####

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
    elif type == "c2":
        a,b,c,d = np.real(val[0]), np.imag(val[0]), np.real(val[1]), np.imag(val[1])
    return a*a + b*b + c*c + d*d

def rkmk_step(Y,y,n,h=1e-10):
    k = np.zeros((s,2,2), dtype="complex128")
    I1 = Y(n)
    k[0] = Y(n)

    for i in range(0,2):
        u = np.zeros((s-2,2,2), dtype="complex128")
        u_tilda = np.zeros((s-2,2,2), dtype="complex128")

        u[i] = h * np.sum([A[i+2,j] * k[j] for j in range(i)], axis=0)
        u_tilda[i] = u[i] + (((c[i+2] * h) / 6) * commutator(I1, u[i]))
        #k[i+1] = Y(matrix_multiply(y, expm(u_tilda[i])), n)
        k[i+1] = Y(n)

    I2 = ((m1 * (k[1] - I1)) + (m2 * (k[2] - I1)) + (m3 * (k[3] - h))) / h
    v = h * np.sum([b[j] * k[j] for j in range(s)], axis=0)
    v_tilda = v + ((h / 4) * commutator(I1,v)) + ((h**2 / 24) * commutator(I2,v))
    y = matrix_multiply(y, expm(v_tilda))
    if not np.isclose(condition_check(y),1.):
        raise ValueError
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
        'A': np.array([[0,complex(0,1)],[complex(0,1),0]]),
        'T': np.array([[0,-1],[1,0]]),
        'G': np.array([[complex(0,1),0],[0,-(complex(0,1))]]),
        'C': np.array([[complex(0,1),0],[0,complex(0,1)]]),
        'N': random.choice(generators)
        }
        self.sequence = sequence

    def run(self):
        self.seqgroup = np.array([self.dict[s] for s in self.sequence])
        self.algebra = np.sum([i*j for i,j in zip(self.seqgroup, np.arange(0,len(self.seqgroup)))],axis=0)
        return self.algebra
