from os import write
from scipy import linalg
import numpy as np
import re
import matplotlib.pyplot as plt
f_in = "E:\\1\\JC_Model_Rho_Mixed.txt" #! address file for input
f1 = open(f_in,"r+")   #! open data file
f_out_S = "E:\\1\\S_L_C.txt" #! Von-Neuman entropy address file for output
f2 = open(f_out_S,"w+") #! open output file
m = 10000 #! m = Number of divisions between x1 and x2 in the fortran program
t = np.arange(0,30,0.003)
s_list = [] # for plot s on based t
l_list = [] # for plot l on based t
c_list = [] # for plot c on based t
for i in range(m):
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
    # print(type(rho_f[0]))
    rho = rho_f
    #? ----------------------------------------------------------------------------
    i = 0
    c_y = []
    while i < 32: #! convert to complex number
        a = float(0)
        b = float(rho[i+1])
        d = complex(a,b)
        c_y = c_y + [d]
        i = i + 2
    # print(c_y)
    #? ----------------------------------------------------------------------------
    #! Von-Neuman Entropy
    c_i = []
    c_i.append((rho[0]+rho[10]))
    c_i.append((rho[1]+rho[11]))
    c_i.append((rho[4]+rho[14]))
    c_i.append((rho[5]+rho[15]))
    c_i.append((rho[16]+rho[26]))
    c_i.append((rho[17]+rho[27]))
    c_i.append((rho[20]+rho[30]))
    c_i.append((rho[21]+rho[31]))
    # print("c_i = ",c_i)
#     c_y = []
#     while i < 9: #! convert to complex number
#         a = float(0)
#         b = float(c_i[i+1])
#         d = complex(a,b)
#         c_y = c_y + [d]
#         i = i + 2
#     landa_1 = 0.5*(c_1 + c_4 + np.sqrt((c_1 * c_1) + (4*c_2*c_3) - (2*c_1*c_4) + (c_4 * c_4))) 
#     print("landa_1 = ",landa_1)
#     landa_2 = 0.5*(c_1 + c_4 - np.sqrt((c_1 * c_1) + (4*c_2*c_3) - (2*c_1*c_4) + (c_4 * c_4))) 
#     print("landa_2 = ",landa_2)
#     S = -(landa_1 * np.log2(landa_1)) - (landa_2 * np.log2(landa_2)) 
#     s_list.append(S)
#     S = str(S)
#     # f2.write(S)  #! write in output file
#     # f2.write("\n")
#     #? ----------------------------------------------------------------------------
#     #! Linear Entropy and Concurrence
#     L = 2 * (1 - ((c_1*c_1) + (2*c_2*c_3) + (c_4 * c_4)))
#     l_list.append(L)
#     Concurrence = np.sqrt(L)
#     c_list.append(Concurrence)
#     L = str(L)
#     Concurrence = str(Concurrence)
#     f2.write(S)
#     f2.write("     ")
#     f2.write(L)
#     f2.write("     ")
#     f2.write(Concurrence)
#     f2.write("\n")
#     #? ----------------------------------------------------------------------------  
#     #? ----------------------------------------------------------------------------
  
# print(s_list)
plt.plot(t,s_list,label = "Von Neuman Entropy")
plt.plot(t,l_list,label = "Linear Entropy")
plt.plot(t,c_list,label = "Concurrence")
plt.legend()
plt.xlabel("T")
plt.ylabel("Entanglement Measure")
plt.show()

    #? ----------------------------------------------------------------------------







