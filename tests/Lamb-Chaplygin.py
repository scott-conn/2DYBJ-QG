from ../model import *
import numpy as np

N = 512
L = 4e5*np.pi
Ls = [L,L]
ns = [N,N]
p = [
    5e8,
    1e8,
    1e-4,
    5e-3,
    2*np.pi/325,
    2*np.pi/84e3,
    0.001,
    30,
    200,
    0.01
]

init_file = 'LC_init_file.h5'
run_sim(Ls,ns,p,init_file)
