import integrator_na_Fraction as IN
from copy import deepcopy as DPC
import sympy as SP
import numpy as np
from fractions import Fraction
import time 

n = m = 3
m0 = 4
n0 = 3
n1 = n0+n+1
m1 = m0+m+1
nm = (n + 1) * (m + 1)
CCC = None
G = []
relu = lambda x: max(x, 0)


triangles=[
    [[0, 0], [0, 1/4], [1/4, 1/4]], #
    [[0, 0], [1/4, 1/4], [1/4, 0]], #
    [[1/4, 0], [1/4, 1/4], [1/2, 1/4]],#
    [[1/4, 0], [1/2, 1/4], [1/2, 0]],#
    [[1/2, 0], [1/2, 1/4], [3/4, 1/4]],#
    [[1/2, 0], [3/4, 1/4], [3/4, 0]],#
    [[3/4, 0], [3/4, 1/4], [1, 1/4]], #
    [[3/4, 0], [1, 1/4], [1, 0]], #
    [[0, 1/4], [0, 1/2], [1/4, 1/2]], #
    [[0, 1/4], [1/4, 1/2], [1/4, 1/4]], #
    [[1/4, 1/4], [1/4, 1/2], [1/2, 1/2]],#
    [[1/4, 1/4], [1/2, 1/2], [1/2, 1/4]],#
    [[1/2, 1/4], [1/2, 1/2], [3/4, 1/2]],#
    [[1/2, 1/4], [3/4, 1/2], [3/4, 1/4]],#
    [[3/4, 1/4], [3/4, 1/2], [1, 1/2]], #
    [[3/4, 1/4], [1, 1/2], [1, 1/4]], #
    [[0, 1/2], [0, 3/4], [1/4, 3/4]], #
    [[0, 1/2], [1/4, 3/4], [1/4, 1/2]], #
    [[1/4, 1/2], [1/4, 3/4], [1/2, 3/4]],#
    [[1/4, 1/2], [1/2, 3/4], [1/2, 1/2]],#
    [[1/2, 1/2], [1/2, 3/4], [3/4, 3/4]],#
    [[1/2, 1/2], [3/4, 3/4], [3/4, 1/2]],#
    [[3/4, 1/2], [3/4, 3/4], [1, 3/4]], #
    [[3/4, 1/2], [1, 3/4], [1, 1/2]], #
    [[0, 3/4], [0, 1], [1/4, 3/4]], #
    [[0, 3/4], [1/4, 1], [1/4, 1/2]], #
    [[1/4, 3/4], [1/4, 1], [1/2, 1]],#
    [[1/4, 3/4], [1/2, 1], [1/2, 3/4]],#
    [[1/2, 3/4], [1/2, 1], [3/4, 1]],#
    [[1/2, 3/4], [3/4, 1], [3/4, 3/4]],#
    [[3/4, 3/4], [3/4, 1], [1, 1]], #
    [[3/4, 3/4], [1, 1], [1, 3/4]] #  
]




for iterat, tr in enumerate(triangles):
    inti = IN.Integrator(tr)
    G_k = [[0 for j in range(nm)] for i in range(nm)]
    G.append(G_k)
    if CCC == None:
        CCC = inti.Cnk_
    else:
        inti.Cnk_ = CCC
    start_time = time.perf_counter()
    for i1 in range(n0, n1):
        for j1 in range(m0, m1):
            for i2 in range(n0, n1):
                for j2 in range(m0, m1):
                    S = [0, 0, 0]
                    if i2 >= 4:
                        S[0] = inti.integrate(i1+i2-4, j1+j2) * (i2 - 3) * (i2 - 2) * (i2 - 1) * (i2)
                    if i2 >= 2 and j2 >= 2:
                        S[1] = 2 * inti.integrate(i1+i2-2, j1+j2-2) * (i2) * (i2 - 1) * (j2 - 1) * (j2)
                    if j2 >= 4:
                        S[0] = inti.integrate(i1+i2, j1+j2-4) * (j2 - 3) * (j2 - 2) * (j2 - 1) * (j2)
                    S = S[0] + S[1] + S[2]
                    G_k[(i1 - n0) * m + (j1 - m0)][(i2 - n0) * m + (j2 - m0)] = S
    end_time = time.perf_counter()
    print(f'Iteration: {iterat}')
    print(f'Iteration time: {end_time - start_time:.6f}')
    print("="*20)
    print(np.array(G_k))
print("="*10, '\n', np.array(G))                    