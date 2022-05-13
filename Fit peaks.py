# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 17:46:17 2022

@author: ismae
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import cmath

def distr(x,Gam, A, a):
    #return A/((x-a)+(Gam/2)*1j) #BW
    return A*np.exp(-(x-a)**2/Gam)
# def distr(x,Gam, A, a):
#     return A/((1+((x-a)**2/Gam)))
def func(x,Gam0,Gam1,Gam2, A0, A1, A2, a0,a1,a2):
    #return abs(distr(x,Gam0, A0, a0) + distr(x,Gam1, A1, a1) + distr(x,Gam2, A2, a2) + distr(x,Gam1, A1, -a1) + distr(x,Gam2, A2, -a2))**2
    return distr(x,Gam0, A0, a0) + distr(x,Gam1, A1, a1) + distr(x,Gam2, A2, a2) + distr(x,Gam1, A1, -a1) + distr(x,Gam2, A2, -a2)
P0 = [1.,1.,1., 12., 3., 0.5, 0., -6.,-13.]

f="C:/Users/ismae/OneDrive/Desktop/Physics/Thesis/Script/Plot Values.csv"
data = np.loadtxt(f, delimiter=",", skiprows=1)
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
data[:,0] = data[:,0]-xmax
p,cov=fit(func,data[:,0],data[:,1],p0=P0)
fig = plt.figure(figsize=(15,15))
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
ax = fig.add_subplot(111)
ax.plot(data[:,0],data[:,1], "k--")
ax.plot(xplt, func(xplt,*p), "b-")
for i in range(0,3):
    ax.plot(xplt, distr(xplt, p[i],p[i+3],p[i+6]), "g-")
    if i>0:
        ax.plot(xplt, distr(xplt, p[i],p[i+3],-1*p[i+6]), "y-")
print(p)