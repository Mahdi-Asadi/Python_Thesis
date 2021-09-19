import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
f_out = "E:\\1\\P_rk4.txt" # address file for output
f2 = open(f_out,"w+")
def dy_dx(x,y):
        wa=1      # atomic frequency  
        wp=0.6    # field frequency
        g=0.6    # coupling strength 
        n_p = 1     # number of photons
        C = np.sqrt(n_p+1)    
        dydx_0= 0
        dydx_1= 0
        dydx_2= -C*y[5]-wp*y[3]
        dydx_3= -C*y[4]-wp*y[2]
        dydx_4= -C*y[3]-wa*y[5]
        dydx_5= +C*y[2]+wa*[4]
        dydx_6= -(wa+wp)*y[7] 
        dydx_7= +(wa+wp)*y[6]
        dydx_8= +C*y[17]+wp*y[9]
        dydx_9= -C*y[16]-wp*y[8]
        dydx_10= +C*(-y[13]+y[19])
        dydx_11= +C*(y[12]-y[18])
        dydx_12= +C*(y[21]-y[11])-(wa-wp)*y[13]
        dydx_13= -C*(y[20]-y[10])+(wa-wp)*y[12]
        dydx_14= +C*y[23]-wa*y[15]
        dydx_15= -C*y[22]+wa*y[14]
        dydx_16= +C*y[9]+wa*y[17]
        dydx_17= -C*y[8]-wa*y[16]
        dydx_18= +C*(y[11]-y[21])-(-wa+wp)*y[19]
        dydx_19= -C*(y[10]-y[20])+(-wa+wp)*y[18]
        dydx_20= +C*(y[13]-y[19])
        dydx_21= -C*(y[12]-y[18])
        dydx_22= +C*y[15]-wp*y[23]
        dydx_23= -C*y[14]+wp*y[22]
        dydx_24= +(wa+wp)*y[25]
        dydx_25= -(wa+wp)*y[24]
        dydx_26= -C*y[29]+wa*y[27]
        dydx_27= +C*y[28]-wa*y[26]
        dydx_28= -C*y[27]+wp*y[29]
        dydx_29= +C*y[26]-wp*y[28]
        dydx_30=  0
        dydx_31=  0
        return [dydx_0,dydx_1,dydx_2,dydx_3,dydx_4,dydx_5,dydx_6,\
                dydx_7,dydx_8,dydx_9,dydx_10,dydx_11,dydx_12,\
                dydx_13,dydx_14,dydx_15,dydx_16,dydx_17,\
                dydx_18,dydx_19,dydx_20,dydx_21,dydx_22,\
                dydx_23,dydx_24,dydx_25,dydx_26,dydx_27,\
                dydx_28,dydx_29,dydx_30,dydx_31]

y_0 = [0.08,0,0.06,0,0.5,0,0.03,0,0.06,0,0.29,0,0.12,0,\
        0.16,0,0.05,0,0.13,0,0.21,0,0.07,0,0.03,0,0.16,0,0.04,0,0.42,0]# initial value


m = 1000
ti = 0
tf = 30
h = tf/m
tspan = np.arange(ti,tf,h)
v = solve_ivp(dy_dx,[0,30],y_0,"RK45",t_eval=tspan) # 4 answer of dydx_1,...,dydx_4
plt.plot(v.t, v.y[0,:],'-', label='r(t)') 
plt.xlabel("x")
plt.ylabel("y")
plt.show()