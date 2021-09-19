from os import write
from scipy import linalg
import numpy as np
import re
import matplotlib.pyplot as plt
f_in = "E:\\1\\JC_Model_Rho_Mixed.txt" #! address file for input
f1 = open(f_in,"r+")   #! open data file
f_out_S = "E:\\1\\S_L_C.txt" #! Von-Neuman entropy address file for output
f2 = open(f_out_S,"w+") #! open output file
m = 10000 # m = Number of divisions between x1 and x2 in the fortran program
ti = 0
tf = 10
dt = tf/m
t = np.arange(ti,tf,dt)

s_list = [] # for plot s on based t
l_list = [] # for plot l on based t
c_list = [] # for plot c on based t
for i in range(m):
    a = t[i]
    # print(a)
    #? ----------------------------------------------------------------------------
    x = f1.readline().split("     ") #! read data
    z  = []
    for i in x: #! delete space in number
        x = re.sub(" ", "", i)
        z.append(x)
    #? ----------------------------------------------------------------------------
    rho = list(z) #! change to list
    rho.pop(0) #! delete first column
    rho.pop(32)#! delete last column
    rho_f = []
    for i in rho:
        rho_f.append(float(i))
    # print(rho_f)
    rho = rho_f
    #? ----------------------------------------------------------------------------
    i = 0
    c_rho = []
    while i < 32: #! convert to complex number
        a = float(rho[i])
        b = float(rho[i+1])
        d = complex(a,b)
        c_rho = c_rho + [d]
        i = i + 2
    # print(c_rho)
    # print(rho1)
    #? ----------------------------------------------------------------------------
    # convert to 4*4
    rho_arr = np.array(c_rho) # convert to array 1D
    rho = rho_arr.reshape(4,4) # convert to 4*4
    rho_2 = np.dot(rho,rho)
    #? ----------------------------------------------------------------------------
    #! Von-Neuman Entropy 
    rho_val,rho_vec = linalg.eig(rho) 
    # print(rho_val)
    rho_real = np.real(rho_val) 
    # print(rho_real)
    rho_list = list(rho_real)
    # # print(rho_list)
    for i in rho_list:
        if i==0:
            rho_list.remove(i)
    rho_real = np.array(rho_list)
    # print(rho_real)
    for i in rho_real:
        if i==0:
            print(rho_real)
    S = 0
    for i in rho_real:
        S = S-(i * np.log2(i))
    s_list.append(S)
    S = str(S)
    f2.write(S)  #! write in output file
    f2.write("\n")
    #? ----------------------------------------------------------------------------
#     # ! Linear Entropy and Concurrence
    L = 2 * (1 - (np.trace(rho_2)))
    l_list.append(L)

#     # Concurrence = np.sqrt(L)
#     # c_list.append(Concurrence)

#     L = str(L)
#     # Concurrence = str(Concurrence)

#     f2.write(S)
#     f2.write("     ")
#     # f2.write(str(landa_1))
#     # f2.write("     ")
#     # f2.write(str(landa_2))
#     # f2.write("     ")
#     f2.write(L)
#     # f2.write("     ")
#     # f2.write(Concurrence)
#     f2.write("\n")
#     #? ----------------------------------------------------------------------------  
#     #? ----------------------------------------------------------------------------
  
# # print(s_list)
plt.plot(t,s_list,label = "Von Neuman Entropy")
# plt.plot(t,l_list,label = "Linear Entropy")
# plt.plot(t,c_list,label = "Concurrence")
plt.legend(title = "g = 0.2",loc = "upper right")
plt.xlabel("T")
plt.ylabel("Entanglement Measure")
plt.show()

    #? ----------------------------------------------------------------------------







