import numpy as np
from scipy.linalg import expm
import math

class LieGroup:
    def __init__(self):
        pass

    def action(self,g,u,type):
        if type == "left":
            return np.array([i*u for i in g])
        elif type == "right":
            return np.array([u*i for i in g])

class LieAlgebra:
    def __init__(self):
        pass

    def exp(self,y):
        return expm(y)

    def commutator(self,a,b,matrix=True):
        if matrix:
            return np.matmul(a,b) - np.matmul(b,a)
        else:
            return a*b - b*a

    def dexpinv(self,a,b,order=4):
        B = bernoulli(order)
        k = 0
        stack = b
        out = B[k] * b
        k=1
        stack = self.commutator(a,stack)
        out += B[k] * stack
        k+=1

        while k<order:
            stack = self.commutator(a,stack)
            out += B[k] / math.factorial(k) * stack
            k+=1
        return out

class Sequences(LieGroup,LieAlgebra):
    def __init__(self,sequences):
        self.dict = {
        'A': np.array([[0,complex(0,1)],[complex(0,1),0]]),
        'T': np.array([[0,1],[-1,0]]),
        'G': np.array([[complex(0,1),0],[0,-(complex(0,1))]]),
        'C': np.array([[complex(0,1),0],[0,complex(0,1)]])
        }
        self.sequence = sequence

    def sequence2group(self):
        self.group = np.array([self.dict[s] for s in self.sequence])
        return

    def group2algebra(self):
        self.algebra = np.sum([i*j for i,j in zip(self.group, np.arange(0,len(self.group)))])
        return

    def run(self):
        self.sequence2group()
        self.group2algebra()
        return self.algebra


class S3Sphere:
    def __init__(self, y=np.array([0, 0, 0, 1])):
        self.n = y.size
        self.y = y

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        if not np.isclose(np.inner(value, value), 1.0):
            raise ValueError(f"y does not lie on the 3-sphere. y^T . y should be one, was {np.inner(value, value)}")
        self._y = value


class RKMK4(LieGroup,LieAlgebra):
    def __init__(self):
        self.a = np.array(
            [
                [0,   0,   0,   0],
                [0.5, 0,   0,   0],
                [0,   0.5, 0,   0],
                [0,   0,   1.0, 0]
            ]
        )
        self.b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])
        self.c = np.array([0, 0.5, 0.5, 1.0])
        self.order = 4
        self.s = 4

    def step(self, f, t, y, h):

        n = y.size
        k = np.zeros((n, self.s))

        for i in range(self.s):
            u = np.zeros(n)

            for j in range(i):
                u += self.a[i, j] * k[:, j]

            u *= h
            k[:, i] = self.dexpinv(u, f(t + self.c[i] * h, self.action(self.exp(u), y), "left"), self.order)

        v = np.zeros(n)
        for i in range(self.s):
            v += self.b[i] * k[:, i]

        return self.action(self.exp(h * v), y,"left")

def solve(f,y0,t_init,t_final,h):

    manifold = S3Sphere(y0)
    timestepper = RKMK4()
    n_steps, last_step = divmod((t_final - t_init), h)
    n_steps = int(n_steps)

    t_array = [t_init + i * h for i in range(n_steps + 1)]

    number_of_cols = n_steps + 1 if np.isclose(last_step, 0) else n_steps + 2

    y_array = np.zeros((len(y0), number_of_cols))

    y_array[:, 0] = y0

    for i in range(1, n_steps + 1):
        manifold.y = timestepper.step(f, t_array[i - 1], manifold.y, h)
        y_array[:, i] = manifold.y

    if not np.isclose(last_step, 0):
        manifold.y = timestepper.step(f, t_array[-1], manifold.y, last_step)
        y_array[:, -1] = manifold.y
        t_array.append(t_final)

    return y_array, t_array
