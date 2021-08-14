from functools import partial
from scipy import linalg
import numpy as np
import re
def par_trs_1(y): # input is a 4*4 matrix

    y_1_1 = y[0:2,0:2] # broken to 2*2 matrix
    y_1_2 = y[0:2,2:4]
    y_2_1 = y[2:4,0:2]
    y_2_2 = y[2:4,2:4]
    y_swap = y_2_1 # transpose type 1
    y_2_1 = y_1_2
    y_1_2 = y_swap
    
    _1 = np.concatenate((y_1_1,y_1_2),axis = 1) #  accession,axis = 1 --->
    _2 = np.concatenate((y_2_1,y_2_2),axis = 1) 

    _all_1 = np.concatenate((_1,_2)) # accession,axis = 0 ||
    return _all_1
def par_trs_2(y): # input is a 4*4 matrix

    y_1_1 = y[0:2,0:2] # broken to 2*2 matrix
    y_1_2 = y[0:2,2:4]
    y_2_1 = y[2:4,0:2]
    y_2_2 = y[2:4,2:4]
    
    y_11_t = np.transpose(y_1_1)
    y_12_t = np.transpose(y_1_2)
    y_21_t = np.transpose(y_2_1)
    y_22_t = np.transpose(y_2_2)
    
    _1 = np.concatenate((y_11_t,y_12_t),axis = 1)
    _2 = np.concatenate((y_21_t,y_22_t),axis = 1) 
    _all_2 = np.concatenate((_1,_2))
    return _all_2
#? --------------------------------------------------------------------------------

f_in = "E:\\1\\JC_Model_Heisenberg_Mixed.txt" # address file for input
f1 = open(f_in,"r+")   # open data file
f_out_1 = "E:\\1\\partial_transpose_type_1 eigenvalues.txt" # address file for partial_transpose_type_1 eigenvalues
f2 = open(f_out_1,"w+") # open output file
f_out_2 = "E:\\1\\partial_transpose_type_2 eigenvalues.txt" # address file for partial_transpose_type_2 eigenvalues
f3 = open(f_out_2,"w+") # open output file
#? --------------------------------------------------------------------------------
m = 10000 # m = Number of divisions between x1 and x2 in the fortran program
for i in range(m):
    #? ----------------------------------------------------------------------------
    # read data
    x = f1.readline().split("     ") # read data
    rho  = []
    for i in x: # delete space in number
        x = re.sub(" ", "", i)
        rho.append(x)
    rho.pop(0) # delete first column
    rho.pop(32)# delete last column
    #? ----------------------------------------------------------------------------
    # convert to complex number
    i = 0
    c_rho = []
    while i < 32: 
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
    #? ----------------------------------------------------------------------------
    # partial transpose
    partial_transpose_type_1 = par_trs_1(rho)
    partial_transpose_type_2 = par_trs_2(rho)
    #? ----------------------------------------------------------------------------
    # calculate eigenvalue and eigenvector
    rho_val_pats_1,rho_vec_pats_1 = linalg.eigh(partial_transpose_type_1)
    rho_val_pats_2,rho_vec_pats_2 = linalg.eigh(partial_transpose_type_2)
    # print("partial_transpose_type_1 = ",rho_val_pats_1)
    # print ("partial_transpose_type_2 = ",rho_val_pats_2)
    #? ----------------------------------------------------------------------------
    # convert to list and real number
    rho_val_1 = list(rho_val_pats_1) # convert array to list
    rho_val_2 = list(rho_val_pats_2) # convert array to list
    rho_val_real1 = np.real(rho_val_1) # convert to real numbers
    rho_val_real2 = np.real(rho_val_2) # convert to real numbers
    # print("rho_val_1 = ",rho_val_1)
    # print("rho_val_2 = ",rho_val_2)
    # print("rho_val_new1 = ",rho_val_real1)
    # print("rho_val_new2 = ",rho_val_real1)
    rho_val_re_n1 = [] # round number 5
    rho_val_re_n2 = []
    for i in rho_val_real1:# convert partial_transpose_type_1 to real number
        v = round(i,4)
        rho_val_re_n1.append(v)
    rho_val_1 = rho_val_re_n1
    for i in rho_val_real2:# convert partial_transpose_type_2 to real number
        v = round(i,4)
        rho_val_re_n2.append(v)
    rho_val_2 = rho_val_re_n2
    # print("rho_value_list = ",rho_value_list)
    # print("rho_value_real = ",rho_value_real)
    # #? ----------------------------------------------------------------------------
    # sort of min to max
    rho_val_1.sort()   
    rho_val_2.sort()  
    # print("rho_val_2 = ",rho_val_2)
    # #? ----------------------------------------------------------------------------

    # #? ----------------------------------------------------------------------------
    # convert to str
    rho_val_1_str = [] 
    rho_val_2_str = [] 
    for i in rho_val_1: 
        rho_val_1_str.append(str(i))
    for i in rho_val_2: 
        rho_val_2_str.append(str(i))
    # print(y_value_real_str)
    # #? ---------------------------------------------------------------------------
    # write in output file for partial_transpose_type_1 eigenvalue
    for i in rho_val_1_str: 
        f2.write(i)
        f2.write("     ")
    f2.write("\n")
    # write in output file for partial_transpose_type_2 eigenvalue
    for i in rho_val_2_str: 
        f3.write(i)
        f3.write("     ")
    f3.write("\n")
    #? ----------------------------------------------------------------------------
f1.close()
f2.close()
f3.close() 