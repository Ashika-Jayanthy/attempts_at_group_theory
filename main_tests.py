import numpy as np
from main import *
import math

num_images = 10
num_pixels = 20
d_lambda = 0.01
lambda_values = np.linspace(0,1,1/d_lambda)

def random_image_array(n_images,n_pixels):
    return np.array([np.random.choice((0,255),n_pixels) for a in range(n_images)])

def tangent_bundle(vector_space):
    return

def gamma(lmbda):
    return np.array([complex(),complex()])

def X_mugamma():
    X = tangent_bundle(random_image_array(num_images,num_pixels))
    return

def ordered_exponential():
    return

def G(t,g,g0):
    g_t = (ordered_exponential(-np.sum([d_lambda * g(gamma(lmbda) * X_mugamma(gamma(lmbda))) for lmbda in lmbda_values])))*g0
    return g_t


g0 = np.array([complex(0,0),complex(0,1)])
t_i = 0
t_f = 5
h = 0.01

solution = solve(G, g0, t_i, t_f, h)
print(solution)
