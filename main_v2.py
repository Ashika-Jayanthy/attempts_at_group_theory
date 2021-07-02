import numpy as np
from scipy.linalg import expm

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
    return matrix_multiply(a,b)

####
s = 4
h = 1
n_sequences = 10

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

def Y(y,t):
    y_t = matrix_multiply(y,Sequence(sequence_array[t]).run())
    return y_t


sequence_array = []
init_sequence = ''.join(random.choice('CGTA') for _ in range(7))
ss = random_sequence_evolve(init_sequence)
sequence_array.append(ss)
for t in range(T):
    ss = random_sequence_evolve(ss)
    sequence_array.append(ss)


y0 = Sequence(init_sequence).run()
y_array = []

for n in range(n_sequences):
    k = np.zeros(s)
    I1, k[0] = Y(y,n)

    for i in range(2,s):
        u = np.zeros(s-2,2,2)
        u_tilda = np.zeros(s-2,2,2)
        u[i] = h * np.sum([A[i,j] * k[j] for j in range(i-1)])
        u_tilda[i] = u[i] + (((c[i] * h) / 6) * commutator(I1, u[i]))
        k[i-1] = Y(matrix_multiply(y, expm(u_tilda[i])), n)

    v = np.zeros(s-2,2,2)
    v_tilda = np.zeros(s-2,2,2)
    I2 = ((m1 * (k[1] - I1)) + (m2 * (k[2] - I1)) + (m3 * (k[3] - h))) / h
    v = h * np.sum([b[j] * k[j] for j in range(s)])
    v_tilda = v + ((h / 4) * commutator(I1,v)) + ((h**2 / 24) * commutator(I2,v)
    y = matrix_multiply(y, expm(v_tilda))
    y_array.append(y)

print(y_array)
