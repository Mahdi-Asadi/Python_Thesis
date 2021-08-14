# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 13:13:45 2021

@author: Mahdi
"""
x = []
t = 0
import random
for i in range(16):
     rand_num = round(random.random(),2)# between 0 to 1
     x = x + [rand_num]
b = [0, 5, 10,15]
n = x[0]+x[5]+x[10]+x[15]
for i in b:
   x[i] = round(x[i]/n,2)
a = x[0]+x[5]+x[10]+x[15]
z = [3,7,11,15]
j = 0
for i in z:
    h = i + 1
    print(x[j:h])
    j = i + 1
print("a = ",a)