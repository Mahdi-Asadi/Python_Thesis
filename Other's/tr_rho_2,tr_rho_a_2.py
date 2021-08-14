from os import write
from scipy import linalg
import numpy as np
import re
import matplotlib.pyplot as plt
from math import sqrt
f_in = "E:\\1\\JC_Model_Schrodinger_py.txt" #! address file for input
f1 = open(f_in,"r+")   #! open data file
f_out_S = "E:\\1\\rho^2.txt" #! Von-Neuman entropy address file for output
f2 = open(f_out_S,"w+") #! open output file
m = 9999 #! m = Number of divisions between x1 and x2 in the fortran program
t = np.arange(0,30,0.0003)

s_list = [] # for plot s on based t
l_list = [] # for plot l on based t
c_list = [] # for plot c on based t
tr_rho_2_list = []
tr_rho_a_2_list = []
for i in range(m):
    a = t[i]
    #? ----------------------------------------------------------------------------
    print(i)
    x = f1.readline().split("     ") #! read data
    z  = []
    for i in x: #! delete space in number
        x = re.sub(" ", "", i)
        z.append(x)
    #? ----------------------------------------------------------------------------
    rho = list(z) #! change to list
    # print(rho)
    rho.pop(0) #! delete first column
    # print(rho)
    rho.pop(4)#! delete last column
    # print(rho)
    rho_f = []
    for i in rho:
        rho_f.append(float(i))
    # print(type(rho_f[0]))
    rho = rho_f
    # print(rho)
    #? ----------------------------------------------------------------------------
    c_1 = complex(rho[0],rho[1])
    c_2 = complex(rho[2],rho[3])
    # print(c_1,c_2)
    comp_c_1 = np.conj(c_1)
    comp_c_2 = np.conj(c_2)
    # print(comp_c_1,comp_c_2)
    rho = np.array([[(c_1.real**2 + c_1.imag**2), c_1 * comp_c_2],[c_2 * comp_c_1, (c_2.real**2 + c_2.imag**2)]])
    # print(rho)
    # ----------------------------------------------------------------------------
    rho_2 = np.dot(rho,rho)
    tr_rho_2 = np.trace(rho_2)
    # print(tr_rho_2)
    
    x = np.real(tr_rho_2)
    x = round(x,2)
    tr_rho_2_list = tr_rho_2_list + [x]
    # print(tr_rho_2_list)
    # -----------------------------------------------------------------------------
    rho_a = np.array([[(c_2.real**2 + c_2.imag**2),0],[0, (c_1.real**2 + c_1.imag**2)]])
    rho_a_2 = np.dot(rho_a,rho_a)
    tr_rho_a_2 = np.trace(rho_a_2)
    x1 = round(tr_rho_a_2,2)
    tr_rho_a_2_list = tr_rho_a_2_list + [x1]    
    # ----------------------------------------------------------------------------
    A = rho_a[0,0]*rho_a[0,0]
    B = rho_a[1,1]*rho_a[1,1] 
    S = -A*np.log2(A)-B*np.log2(B)
    L = 2*(1-(A**2)-(B**2))
    # ----------------------------------------------------------------------------
    
    f2.write(str(x))
    f2.write("     ")
    f2.write(str(a))
    f2.write("     ")
    f2.write(str(S))
    f2.write("     ")
    f2.write(str(L))
    f2.write("\n")
    #? ----------------------------------------------------------------------------  
    #? ----------------------------------------------------------------------------
  
# # # print(s_list)
plt.plot(t,tr_rho_2_list,label = "tr_rho_2")
plt.plot(t,tr_rho_a_2_list,label = "tr_rho_a_2")
# # plt.plot(t,c_list,label = "Concurrence")
plt.legend()
plt.xlabel("T")
plt.ylabel("Entanglement Measure")
plt.show()