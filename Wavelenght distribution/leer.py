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
from PIL import Image as im

c = 299792458 #m/s
eps= 8.8541878128e-12 #F/m
e=1.60217662e-19 #C
u=1.66053906660e-27 #Kg
m_e= 9.10938356e-31 #kg
KB = 1.38064852e-23 #J/K 
h = 6.626176e-34 #J s 
hbar = 1.0545718e-34
m_n = 1.00866491588
g=9.80665
lambda_par=1595.7292122046995	#+/-	147.471394720765
mu=0.004632543663155012	#+/-	5.46776175965519e-05
sigma=0.0006873016655595522	
l_max = 3.7
Dl=l_max*0.4/2
z0=2.15
def v(l):
    return h/(m_n*u*l)*1e9
# def rho(l):
#     return np.exp(-0.5*(l-l_max)**2/Dl**2)/(Dl*(2*np.pi)**0.5)
def func(l,A,mu,sig):
    return A/(2.)*np.exp(A/(2.)*(2.*mu+A*sig**2-2*l))
def rho(l,A,mu,sig):
    return func(l,A,mu,sig)*erfc((mu+A*sig**2-l)/(np.sqrt(2)*sig))
def distr(x,y,l,A,x0,y0,th):
    y0 = y0-1/2*g*(z0)**2/v(l)**2
    phi = np.arctan(y0/z0)
    sx =th*z0
    sy=th*((z0**2+y0**2)**0.5+phi*(y-y0))
    return A/(sx*sy)*np.exp(-0.5*((x-x0)**2/sx**2+(y-y0)**2/sy**2))
    # return A/((1+((x-a)**2*Gam)))
    # return A*(b/(x-a))**2*np.exp(-(b/(x-a))**2*Gam)
# def func(l,A,mu,sig):
#     return A/(2.)*np.exp(A/(2.)*(2.*mu+A*sig**2-2*l))
# def distr(l,A,mu,sig,a,b):
#           #return func(b/(l-a),A,mu,sig)*erfc((mu+A*sig**2-b/(l-a))/(np.sqrt(2)*sig))
#           return func(b*np.sqrt(abs(l-a)),A,mu,sig)*erfc((mu+A*sig**2-b*np.sqrt(abs(l-a)))/(np.sqrt(2)*sig))
#           #return func((l-a),A,mu,sig)*erfc((mu+A*sig**2-(l-a))/(np.sqrt(2)*sig))
#p= [3.07584609, 3# print(p1)
# print(np.diag(cov1)**0.5).71207594, 9.84870124, 0.09066904, 10.]
P0=[1., 0., 0., 0.001]
xy=np.linspace(-128e-3, 128e-3, 1280)
# X, Y = np.meshgrid(xy, xy)
# Z = distr(X,Y,*P0)
matrix = np.zeros((1280,1280))
lam=np.linspace(1., 15., 100)
for l in range(len(lam)):
    a=rho(lam[l]/1e3,lambda_par, mu, sigma)
    for i in range(600,800):
        for j in range(400,800):
            if abs(abs(np.arctan((xy[j])/z0))-lam[l]*1e-3)<2.e-3:
                matrix[i][j+int(np.sign(xy[j])*lam[l]*z0*5)] += a*distr(xy[j],-xy[i],lam[l],*P0)
            else:
                matrix[i][j] += a*distr(xy[j],-xy[i],lam[l],*P0)
m_max = np.amax(matrix)
matrix= matrix/m_max
# with open("leer simu.mpa", 'w') as f:
#         np.savetxt(f, matrix)
# for i in range(500):
#         for j in range(500):
#             matrix[i][j] += distr(xy[j],xy[i],5,*P0)  
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
image = im.fromarray((matrix))
ax.set_xlim([550,750])
ax.set_ylim([800,600])
im=ax.imshow(matrix,cmap='plasma')
image.save("leer diff simu.tiff")
