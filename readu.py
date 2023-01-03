# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 16:33:39 2023

@author: Administrator
"""

import numpy as np
import math
import matplotlib.pyplot as plt

t = 1 
#e = 2.718281828459045
u_FTCS = np.loadtxt("1d_FTCS_u.dat")
u_BTCS = np.loadtxt("1d_BTCS_u.dat")
x = np.linspace(-1,1,num=81)

u_ana = np.linspace(-1,1,num=81)
u_diff_FTCS = np.linspace(-1,1,num=81)
u_diff_BTCS = np.linspace(-1,1,num=81)
for i in range(81):
    u_ana[i] = -math.pow(math.e, -t)*math.sin(math.pi*x[i])
    u_diff_FTCS[i] = abs(u_FTCS[i]-u_ana[i])
    u_diff_BTCS[i] = abs(u_BTCS[i]-u_ana[i])

plt.figure(figsize=(12, 5),dpi=300)
plt.subplot(121)
#plt.plot(x,u,'b--',label="FTCS simulation") #simulation
plt.plot(x,u_FTCS,'k--',label="FTCS simulation") #simulation
plt.plot(x,u_BTCS,'b--',label="BTCS simulation") #simulation
plt.plot(x,u_ana,'g-',label="analytical solution")
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.title("Solution Field")

plt.subplot(122)
plt.plot(x,u_diff_FTCS,label="FTCS simulation")
plt.plot(x,u_diff_BTCS,label="BTCS simulation")
plt.xlabel("x")
plt.ylabel("Difference $\epsilon$")
plt.title("Discretization error")
plt.legend()

plt.tight_layout()