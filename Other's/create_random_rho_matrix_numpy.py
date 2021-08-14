import numpy as np
from numpy import random, true_divide
from scipy import linalg
import time

start = time.time()
t = True
j = True
k = True
l = True
m = True
while j or k or l or m or t:
    x = random.random((4,4)) #! random matrix between 0 - 1
    rounded_matrix = np.round_(x, decimals = 2)
    y = np.trace(rounded_matrix)
    rounded_matrix[1,0] = rounded_matrix[0,1]
    rounded_matrix[2,0] = rounded_matrix[0,2]
    rounded_matrix[3,0] = rounded_matrix[0,3]
    rounded_matrix[2,1] = rounded_matrix[1,2]
    rounded_matrix[3,1] = rounded_matrix[1,3]
    rounded_matrix[3,2] = rounded_matrix[2,3]
    a,b = linalg.eig(rounded_matrix)
    if y == 1 :
        t = False
        # print(y,": ",t)
    else:
        continue
    if np.real(a[0]) >= 0:
        # print("a[0] = ",np.real(a[0]))
        j = False
    else:
        continue
    if np.real(a[1]) >= 0:
        # print("a[1] = ",np.real(a[1]))
        k = False
    else:
        continue
    if np.real(a[2]) >= 0:
        # print("a[2] = ",np.real(a[2]))
        l = False
    else:
        continue
    if np.real(a[3]) >= 0:
        # print("a[3] = ",np.real(a[3]))
        m = False
    else:
        continue
    print(a)
print("a = ", a)   
print(rounded_matrix)
rho = rounded_matrix
rho_2 = np.dot(rho,rho)
tr_rho_2 = np.trace(rho_2)
print("tr_rho = ",y)
print("tr_rho_2 = ", tr_rho_2)
print("Run Time: " + str( time.time() - start ))


