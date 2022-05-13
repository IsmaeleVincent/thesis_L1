# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 21:34:28 2022

@author: ismae
"""

from scipy.special import erfc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import cmath
from scipy.stats import chisquare as cs

def distr(x,A,x0, s0, tg_th):
    sx = s0+tg_th*(x-x0)
    return A/(sx)*np.exp(-0.5*(x-x0)**2/sx**2)
P0=[ 0.196986777,  10.097921, 0.79071826, 0.]
#P0=[0.86315847, 2e-3*9.57562695, 0.06995962, 2e-3*0.70166142]
f="Horiz Values 0deg.csv"
data = np.loadtxt(f, skiprows=1, delimiter=",")
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
ymax = np.amax(data[:,1])
ymin = np.amin(data[:,1])
data[:,0] = (data[:,0]-xmax+10)
#data[:,0] = (data[:,0]-xmax+10)*2e-3
#data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
p,cov=fit(distr,data[:,0],data[:,1],p0=P0)
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
ax.plot(data[:,0],data[:,1], "ro")
ax.plot(xplt, distr(xplt,*p), "b-")
print("p=",p)
print("cov=", np.diag(cov)**0.5)
f="Horiz Values offBragg.csv"
data = np.loadtxt(f, skiprows=1, delimiter=",")
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
ymax = np.amax(data[:,1])
ymin = np.amin(data[:,1])
data[:,0] = (data[:,0]-xmax+10)
#data[:,0] = (data[:,0]-xmax+10)*2e-3
#data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
p,cov=fit(distr,data[:,0],data[:,1],p0=P0)
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
ax.plot(data[:,0],data[:,1], "go")
ax.plot(xplt, distr(xplt,*p), "y-")
print("p=",p)
print("cov=", np.diag(cov)**0.5)
f="Horiz Values leer.csv"
data = np.loadtxt(f, skiprows=1, delimiter=",")
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
ymax = np.amax(data[:,1])
ymin = np.amin(data[:,1])
data[:,0] = (data[:,0]-xmax+10)
#data[:,0] = (data[:,0]-xmax+10)*2e-3
#data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
p,cov=fit(distr,data[:,0],data[:,1],p0=P0)
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
ax.plot(data[:,0],data[:,1], "ko")
ax.plot(xplt, distr(xplt,*p), "r-")
print("p=",p)
print("cov=", np.diag(cov)**0.5)
print(p[2]*2e-3/2.05)