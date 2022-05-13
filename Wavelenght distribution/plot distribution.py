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

# def distr(x, A, x0,sx,tg):
#     return A/(sx+tg*(x-x0))*np.exp(-(x-x0)**2/(2*(sx+tg*(x-x0))**2))
    #return A*np.exp(-((b*x-a))**2*Gam)
    #return A/((1+((x-a)**2*Gam)))
    #return A*(b/(x-a))**2*np.exp(-(b/(x-a))**2*Gam)
def func(l,A,mu,sig):
    return A/(2.)*np.exp(A/(2.)*(2.*mu+A*sig**2-2*l))
def distr(l,A,mu,sig,a,b):
          # return func(b/(l-a),A,mu,sig)*erfc((mu+A*sig**2-b/(l-a))/(np.sqrt(2)*sig))
          return func(b*np.sqrt(abs(l-a)),A,mu,sig)*erfc((mu+A*sig**2-b*np.sqrt(abs(l-a)))/(np.sqrt(2)*sig))
          #return func((l-a),A,mu,sig)*erfc((mu+A*sig**2-(l-a))/(np.sqrt(2)*sig))
# #def distr(x,Gam, A, a):
#     return A/((1+((x-a)**2/Gam)))
# def func(x,Gam0,Gam1,Gam2, A0, A1, A2, a0,a1,a2):
#     #return abs(distr(x,Gam0, A0, a0) + distr(x,Gam1, A1, a1) + distr(x,Gam2, A2, a2) + distr(x,Gam1, A1, -a1) + distr(x,Gam2, A2, -a2))**2
#     return distr(x,Gam0, A0, a0) + distr(x,Gam1, A1, a1) + distr(x,Gam2, A2, a2) + distr(x,Gam1, A1, -a1) + distr(x,Gam2, A2, -a2)
#P0 = [13.396, -2.298,  0.383,  6.000,  -8.557]
#P0 = [1.66, 9.29, 0.95, 0.44, 3.29]
P0 = [87.6770335,   5.5756957,   0.3880741 ,  2.42898191,  2.00215552]
#P01 = [78.1141,  5.53583488,  0.582523, 2.50755852,  2.010684 ]
#P0 = [2.70, 4.52, 0.28, 5.46]
#P0 = [0.26,  0.48, 5.27, 0.9]
#P0 = [21.59, 347.1, 1.100, 0.3216, 621.2814]
#P0 = [0.73974169, 10.00004964,  0.73957559,  0.17626706,0.1]
f="Vert Values leer simu.csv"
f1="Vert Values leer.csv"
data = np.loadtxt(f, delimiter=",", skiprows=1)
xmax = data[:,0][data[:,1]==np.amax(data[:,1])]
ymax = np.amax(data[:,1])
ymin = np.amin(data[:,1])
data[:,0] = (data[:,0]-xmax+10)
data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
p,cov=fit(distr,data[:,0],data[:,1], p0=P0)
xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
data1 = np.loadtxt(f1, delimiter=",", skiprows=1)
xmax1 = data1[:,0][data1[:,1]==np.amax(data1[:,1])]
ymax1 = np.amax(data1[:,1])
ymin1 = np.amin(data1[:,1])
data1[:,0] = data1[:,0]-xmax1+10
data1[:,1] = (data1[:,1]-ymin1)/ymax1
p1,cov1=fit(distr,data1[:,0],data1[:,1])
xplt1=np.linspace(data1[:, 0][0], data1[:, 0][-1], 1000)
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
ax.plot(data1[:,0],data1[:,1], "ro")
ax.plot(xplt1, distr(xplt1,*p1), "g-")
ax.plot(data[:,0],data[:,1], "ko")
ax.plot(xplt, distr(xplt,*p), "b-")
#ax.set_xlim([6,16])
#chi=cs(data[:,1], distr(data[:,0],*p))
print(p)
print(np.diag(cov)**0.5)
print(p1)
print(np.diag(cov1)**0.5)

# c = 299792458 #m/s
# eps= 8.8541878128e-12 #F/m
# e=1.60217662e-19 #C
# u=1.66053906660e-27 #Kg
# m_e= 9.10938356e-31 #kg
# KB = 1.38064852e-23 #J/K 
# h = 6.626176e-34 #J s 
# hbar = 1.0545718e-34
# m_n = 1.00866491588

# b = h/(m_n*u*1.9)*(4*10**-3/9.81)**0.5
# print(b)
# print(h/(m_n*u*(5e-9)))
