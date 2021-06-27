import numpy as np
from main import *
import math
import matplotlib.pyplot as plt
from scipy.linalg import expm

genomic_sequences = ['ATG','TGC','AAA']

def random_image_array(num_images = 10,num_pixels = 20):
    return np.array([np.random.choice((0,255),n_pixels) for a in range(n_images)])


def X():
    return

def G(t,g):
    alg = Sequences(genomic_sequences[t])
    grp = expm(alg)
    g_t = np.array([grp[0,0], grp[1,0]])
    return g_t


g0 = np.array([complex(0,0),complex(0,1)])
t_i = 0
t_f = 2
h = 1

solution = solve(G, g0, t_i, t_f, h)
print(solution[0])

for x in solution[0]:
    plt.polar([0,np.angle(x[1])],[0,np.abs(x[1])],marker='o')
plt.show()
