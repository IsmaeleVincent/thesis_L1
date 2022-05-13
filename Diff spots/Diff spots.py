#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:40:43 2022

@author: aaa
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 17:46:17 2022

@author: ismae
"""
from scipy.special import erfc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import cmath
from scipy.stats import chisquare as cs
def gauss(x, A, x0,sx):
    return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))
def distr(x,A0,A1,A2,x0,x1,x2,s0,s1,s2):
    return gauss(x,A0,x0,s0)+gauss(x, A1, x1,s1)+gauss(x, A2, x2,s2)
f="diff spots.csv"
data = np.loadtxt(f, delimiter=",", skiprows=1)
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
ymax = np.amax(data[:,1])
ymin = np.amin(data[:,1])
data[:,0] = (data[:,0]-xmax)
data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
P0 = [1.,1.,0.5, 0., 5., 10., 0.5, 1.5, 2.]
bound = ([0,0,0,-3.,3.,3.,0.1,0.,0.],[1.1,1.5,1.5,3.,50.,50., 50.,50.,50.])
p,cov=fit(distr,data[:,0],data[:,1], p0=P0, bounds = bound)
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
ax.plot(data[:,0],data[:,1], "ko")
ax.plot(xplt, distr(xplt,*p), "b--")
color=["r-","g-","k-"]
for i in range(3):
    ax.plot(xplt, gauss(xplt, p[i], p[i+3], p[i+6]), color[i%3])

#ax.set_xlim([6,16])
#chi=cs(data[:,1], distr(data[:,0],*p))
print(p)
print(np.diag(cov)**0.5)