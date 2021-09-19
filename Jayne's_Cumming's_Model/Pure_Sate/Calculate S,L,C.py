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
    c_y = []
    while i < 32: #! convert to complex number
        a = float(rho[i])
        b = float(rho[i+1])
        d = complex(a,b)
        c_y = c_y + [d]
        i = i + 2
    # print(c_y)
    rho_con = np.array(c_y)
    rho1 = rho_con.reshape(4,4)
    # print(rho1)
    #? ----------------------------------------------------------------------------
    #! Von-Neuman Entropy
    c_1 = rho1[0,0] + rho1[1,1]
    # print(rho[5]+rho[16])
    c_2 = rho1[0,2] + rho1[1,3]
    c_3 = rho1[2,0] + rho1[3,1]
    c_4 = rho1[2,2] + rho1[3,3]
    c_1 = np.real(c_1)
    c_2 = np.real(c_2)
    c_3 = np.real(c_3)
    c_4 = np.real(c_4)
    # print("c_1 = ",c_1,"\nc_2 = ",c_2,"\nc_3 = ",c_3,"\nc_4 = ",c_4,)  
    landa_1 = 0.5*(c_1 + c_4 + np.sqrt((c_1 * c_1) + (4*c_2*c_3) - (2*c_1*c_4) + (c_4 * c_4)))
    landa_2 = 0.5*(c_1 + c_4 - np.sqrt((c_1 * c_1) + (4*c_2*c_3) - (2*c_1*c_4) + (c_4 * c_4))) 
    S = -(landa_1 * np.log2(landa_1)) - (landa_2 * np.log2(landa_2)) 
    s_list.append(S)
    S = str(S)
    # f2.write(S)  #! write in output file
    # f2.write("\n")
    #? ----------------------------------------------------------------------------
    #! Linear Entropy and Concurrence
    L = 2 * (1 - ((c_1*c_1) + (2*c_2*c_3) + (c_4 * c_4)))
    l_list.append(L)
    # Concurrence = np.sqrt(L)
    # c_list.append(Concurrence)
    L = str(L)
    # Concurrence = str(Concurrence)
    f2.write(S)
    f2.write("     ")
    f2.write(str(landa_1))
    f2.write("     ")
    f2.write(str(landa_2))
    f2.write("     ")
    f2.write(L)
    # f2.write("     ")
    # f2.write(Concurrence)
    f2.write("\n")
    #? ----------------------------------------------------------------------------  
    #? ----------------------------------------------------------------------------
  
# print(s_list)
plt.plot(t,s_list,label = "Von Neuman Entropy")
plt.plot(t,l_list,label = "Linear Entropy")
# plt.plot(t,c_list,label = "Concurrence")
plt.legend()
plt.xlabel("T")
plt.ylabel("Entanglement Measure")
plt.show()

    #? ----------------------------------------------------------------------------







