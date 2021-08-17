from scipy import linalg
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import re
#? --------------------------------------------------------------------------------
f_in = "E:\\1\\Rabi_Model_for rho_Mixed.txt" # address file for input
f1 = open(f_in,"r+")   # open data file
f_out = "E:\\1\\Concurrence.txt" # address file for output
f2 = open(f_out,"w+") # open output file
t = np.arange(0,30,0.003)
C_list = []
#? --------------------------------------------------------------------------------
m = 10000 # m = Number of divisions between x1 and x2 in the fortran program
for i in range(m):
    # ----------------------------------------------------------------------------
    # read data
    x = f1.readline().split("     ") # read data
    rho  = []
    for i in x: # delete space in number
        x = re.sub(" ", "", i)
        rho.append(x)
    rho.pop(0) # delete first column
    rho.pop(32)# delete last column
    # print("rho = ", rho)
    #? ----------------------------------------------------------------------------
    # convert to complex number
    i = 0
    c_rho = []
    while i < 32: 
        # print("i = ",i)
        a = float(rho[i])
        b = float(rho[i+1])
        d = complex(a,b)
        c_rho = c_rho + [d]
        i = i + 2
    #? ----------------------------------------------------------------------------
    # convert to 4*4
    rho_arr = np.array(c_rho) # convert to array 1D
    rho = rho_arr.reshape(4,4) # convert to 4*4
    # y_array_reshape_real = np.real(y_array_reshape) # convert to real matrix
    rho_cmpconj = np.conj(rho)
    # print ("y_cmpconj =  ", y_cmpconj)
    #? ----------------------------------------------------------------------------
    # calculate Xi
    sigma_y = np.array([[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]])
    _1_ = np.dot(rho,sigma_y) # rho*. sigma_y
    # print("_1_ = ",_1_)
    _2_ = np.dot(_1_,rho_cmpconj) # sigma_y. rho*. sigma_y
    # print("_2_ = ",_2_)
    Xi = np.dot(_2_,sigma_y)# rho.sigma_y. rho*. sigma_y
    # print("_3_ = ",_3_)
    #? ----------------------------------------------------------------------------
    # calculate eigenvalue and eigenvector
    rho_val,rho_vec = linalg.eig(Xi) 
    # print("y Eigenvalue = ",y_value)
    # print ("y Eigenvector = ",y_vector)
    #? ----------------------------------------------------------------------------
    # convert to list and real number
    rho_val_list = list(rho_val) # convert array to list
    rho_val_new = np.real(rho_val_list) # convert to real numbers
    # print("rho_val_kist = ",rho_val_list)
    rho_val_real = [] # round number 4
    for i in rho_val_new:
        v = round(i,4)
        rho_val_real.append(v)
    rho_val_list = rho_val_real
    # print("rho_value_list = ",rho_value_list)
    # print("rho_value_real = ",rho_value_real)
    #? ----------------------------------------------------------------------------
    # sort of max to min
    rho_val_list.sort(reverse=True)   
    # print("y_value_real = ",rho_value_list)
    #? ----------------------------------------------------------------------------
    # calculate c(rho)
    # print("np.sqrt(rho_val_list[0]) =",(rho_val_list[0]),"np.sqrt(rho_val_list[1]) =",(rho_val_list[1])
    #     ,"np.sqrt(rho_val_list[2]) =",(rho_val_list[2]),"np.sqrt(rho_val_list[3]) =",(rho_val_list[3]))
    c_rho = sqrt(rho_val_list[0])-sqrt(rho_val_list[1])-sqrt(rho_val_list[2])-sqrt(rho_val_list[3])
    # print("np.sqrt(rho_val_list[0]) =",sqrt(rho_val_list[0]),"np.sqrt(rho_val_list[1]) =",sqrt(rho_val_list[1])
    #     ,"np.sqrt(rho_val_list[2]) =",sqrt(rho_val_list[2]),"np.sqrt(rho_val_list[3]) =",sqrt(rho_val_list[3]),"c_rho = ",c_rho)
    c_rho = round(c_rho,5)
    if c_rho < 0:
        rho_val_list.append(0)
        C_list.append(0)
    else:
        rho_val_list.append(c_rho)
        C_list.append(c_rho)
    #? ----------------------------------------------------------------------------
    # convert to str
    rho_val_list_str = [] 
    for i in rho_val_list: 
        rho_val_list_str.append(str(i))
    # print(y_value_real_str)
    #? ---------------------------------------------------------------------------
    # write in output file
    for i in rho_val_list_str: 
        f2.write(i)
        f2.write("      ")
    f2.write("\n")
    #? ----------------------------------------------------------------------------
f1.close()
f2.close() 
plt.plot(t,C_list)
plt.xlabel("T")
plt.ylabel("Concurence")
plt.title("Concurence")
plt.show()