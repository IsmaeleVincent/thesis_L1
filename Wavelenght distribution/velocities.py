# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:09:13 2022

@author: ismae
"""
import numpy as np
import matplotlib.pyplot as plt

c = 299792458 #m/s
eps= 8.8541878128e-12 #F/m
e=1.60217662e-19 #C
u=1.66053906660e-27 #Kg
m_e= 9.10938356e-31 #kg
KB = 1.38064852e-23 #J/K 
h = 6.626176e-34 #J s 
hbar = 1.0545718e-34
m_n = 1.00866491588

def v(l):
    return h/(m_n*u*l)*1e9
def Dy(l):
    return 1/2*9.8*(1.9)**2/v(l)**2

v_max = h/(m_n*u*2)*1e9
xplt=np.linspace(1, 6, 1000)
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
#ax.plot(xplt, h/(m_n*u*xplt)*1e9, "ro")
ax.plot(xplt,Dy(xplt), "k-")
ax.plot(xplt,Dy(1)+0.002+0*xplt, "b-")
ax.plot(xplt,Dy(1)+0.001+0*xplt, "r-")