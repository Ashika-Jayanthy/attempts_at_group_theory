import numpy as np
from scipy.linalg import expm
import math
from scipy.special import bernoulli


def c2_to_matrix(val):
    z1,z2 = val[0],val[1]
    return np.array([[z1,-np.conjugate(z2)],[z2,np.conjugate(z1)]])

def matrix_to_c2(m):
    return np.array([m[0,0],m[1,0]])

def condition_check(val, type="matrix"):
    if type == "matrix":
        a,b,c,d = np.real(val[0,0]), np.imag(val[0,0]), np.real(val[1,0]), np.imag(val[1,0])

    elif type == "c2":
        a,b,c,d = np.real(val[0]), np.imag(val[0]), np.real(val[1]), np.imag(val[1])

    return a*a + b*b + c*c + d*d

class LieGroup:
    def __init__(self):
        pass

    def action(self,g,u,type):
        if type == "left":
            return np.array([i.dot(u) for i in g])
        elif type == "right":
            return np.array([i.dot(u) for i in g.T])


class LieAlgebra:
    def __init__(self):
        pass

    def exp(self,y):
        return expm(y)

    def commutator(self,a,b):
        return np.matmul(a,b) - np.matmul(b,a)


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

class SO3(LieAlgebra):
    def __init__(self):
        pass

    def Amatrix(self, y):
        b, c, d = y[0], y[1], y[2]
        return np.array([[0, -d, c], [d, 0, -b], [-c, b, 0]])

    def Bmatrix(self,y):
        b, c, d = y[0], y[1], y[2]
        return np.array([[b*b, b*c, b*d], [b*c, c*c, c*d], [b*d, c*d, d*d]])

    def rodrigues_exp(self, y):
        theta = np.linalg.norm(y)
        A = self.Amatrix(y)
        B = self.Bmatrix(y)
        rexp = np.cos(theta)*np.eye(3) + (np.sin(theta) / theta) * A + ((1 - np.cos(theta)) / theta ** 2) * B
        return rexp

    def dexpinv(self, u, v):
        theta = np.linalg.norm(u)
        vv = np.array([v[2, 1], v[0, 2], v[1, 0]])
        uu = self.Amatrix(u)
        lhs = np.eye(3) - 0.5 * uu + (2 - theta / np.tan(0.5 * theta)) / (2 * theta ** 2) * uu
        return self.action(lhs, vv, "left")




class Sequence:
    def __init__(self,sequence):
        self.dict = {
        'A': np.array([[0,complex(0,1)],[complex(0,1),0]]),
        'T': np.array([[0,-1],[1,0]]),
        'G': np.array([[complex(0,1),0],[0,-(complex(0,1))]]),
        'C': np.array([[complex(0,1),0],[0,complex(0,1)]])
        }
        self.sequence = sequence

    def sequence2group(self):
        self.group = np.array([self.dict[s] for s in self.sequence])
        return

    def group2algebra(self):
        self.algebra = np.sum([i*j for i,j in zip(self.group, np.arange(0,len(self.group)))],axis=0)
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


class C2:
    def __init__(self, y = np.array([complex(0,0), complex(0,1)])):
        self.n = y.size
        self.y = y

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        if not np.isclose(np.square(np.abs(value[0])) + np.square(np.abs(value[1])), 1.0):
            raise ValueError
        self._y = value


class RKMK4(LieGroup,LieAlgebra):
    def __init__(self):
        self.a = np.array(
            [
                [0,   0,   0,   0],
                [0.5, 0,   0,   0],
                [0,   0.5, 0,   0],
                [0,   0,   1.0, 0]
            ], dtype="complex128"
        )
        self.b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6], dtype="complex128")
        self.c = np.array([0, 0.5, 0.5, 1.0])
        self.order = 4
        self.s = 4

    def step(self, func, t, y, h):

        n = y.size
        k = np.zeros((self.s,n,n),dtype="complex128")

        for i in range(self.s):
            u = np.zeros((n,n),dtype="complex128")
            for j in range(i):
                u += self.a[i, j] * k[j, :]
            u *= h
            k[i:] = self.dexpinv(u, func(t + self.c[i] * h, self.action(self.exp(u), y, "left")), self.order)
        v = np.zeros((n,n),dtype="complex128")
        for i in range(self.s):
            v += self.b[i] * k[i, :]

        return self.action(self.exp(h * v), y, "left")

def solve(func,y0,t_init,t_final,h):

    manifold = C2(y0)
    timestepper = RKMK4()
    n_steps, last_step = divmod((t_final - t_init), h)
    n_steps = int(n_steps)

    t_array = [t_init + i * h for i in range(n_steps + 1)]

    number_of_cols = n_steps + 1 if np.isclose(last_step, 0) else n_steps + 2

    y_array = np.zeros((number_of_cols, len(y0)),dtype="complex128")

    y_array[0,:] = y0

    for i in range(1, n_steps + 1):
        g = manifold.y
        manifold.y = timestepper.step(func, t_array[i - 1], manifold.y, h)
        y_array[i,:] = manifold.y

    if not np.isclose(last_step, 0):
        manifold.y = timestepper.step(func, t_array[-1], manifold.y, last_step)
        y_array[-1,:] = manifold.y
        t_array.append(t_final)

    return y_array, t_array
