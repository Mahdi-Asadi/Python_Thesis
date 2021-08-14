# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 18:21:22 2021

@author: Mahdi
"""
# calculate Trace rho^2
import numpy as np
rho = np.array([[0.08,0.06,0.05,0.03],[0.06,0.29,0.13,0.16]
               ,[0.05,0.13,0.21,0.07],[0.03,0.16,0.07,0.42]])
rho_2 = np.dot(rho,rho)
tr_rho_2 = np.trace(rho_2)
print("rho_2 = \n",rho_2)
print("tr_rho_2 = \n", tr_rho_2)