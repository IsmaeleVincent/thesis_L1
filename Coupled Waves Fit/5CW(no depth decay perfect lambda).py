 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 15:09:28 2022

@author: aaa
"""
"""
This module defines the vector field for 5 coupled wave equations
(without a decay, Uchelnik) and first and second harmonics in the modulation; phase: 0 or pi (sign of n2).
Fit parameters are: n1,n2, d, and wavelength; 
Fit 5/(5) orders!
!!!Data: X,order,INTENSITIES
Fit  background for second orders , first and subtract it for zero orders (background fixed)
"""
from scipy.integrate import ode
from scipy import integrate
import numpy as np
from numpy.linalg import eig,solve
import inspect,os,time
from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.special import erfc
import matplotlib.pyplot as plt
import matplotlib as mpl
import socket
import shutil
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
import math
pi=np.pi
rad=pi/180

sorted_fold_path="/home/aaa/Desktop/Thesis/Script/Trial/Sorted data/" #insert folder of sorted meausements files
allmeasurements = sorted_fold_path+"All measurements/"
allrenamed = allmeasurements +"All renamed/"
allmatrixes = allmeasurements + "All matrixes/"
allpictures = allmeasurements + "All pictures/"
allrawpictures = allmeasurements + "All raw pictures/"
alldata_analysis = allmeasurements + "All Data Analysis/"
allcropped_pictures = alldata_analysis + "All Cropped Pictures/"
allcontrolplots = alldata_analysis + "All Control plots/"
allcontrolfits = alldata_analysis + "All Control Fits/"
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg")
foldername.append("79-4U_77c88deg")
foldername.append("79-8U_76c76deg")
foldername.append("79-12U_75c64deg")
foldername.append("79-16U_74c52deg")
tilt=[0,40,48,61,69,71,79,80,81,79,79,79,79]
n_theta=[26,46,28,17,16,20,21,20,19,48,43,59,24]  #number of measurements files for each folder (no flat, no compromised data)
n_pixel = 16384 #number of pixels in one measurement
"""
This block fits the diffraction efficiencies n(x)= n_0 + n_1 cos(Gx)
"""

n_diff= 2 #number of peaks for each side, for example: n=2 for 5 diffracted waves
lam= 3.5e-3 #incoming wavelenght in micrometers
LAM= 0.5 #grating constant in micrometers
b=2*pi/lam #beta value 
G=2*pi/LAM
bcr1=7.0#scattering lenght x density
bcr2=0.
n_0 =1.00
n_1 = bcr1*2*pi/b**2 
#print(n_1)

def k_jz(theta, j, G):
    k_jz=b*(1-(np.sin(theta)-j*G/b)**2)**0.5
    return k_jz
def dq_j (theta, j, G):
    return b*np.cos(theta) - k_jz(theta, j, G)
for k in range(1,2):#len(foldername)):
    d=78/np.cos(tilt[k]*rad)
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff.mpa',skiprows=1)    
    th=np.linspace(diff_eff[0,0],diff_eff[-1,0], 1000)*rad
    S=np.zeros((2*n_diff+1,len(th)),dtype=np.complex)
    sum_diff = np.zeros(len(th)) 
    for t in range(len(th)):
        A = np.zeros((2*n_diff+1,2*n_diff+1), dtype=np.complex)
        for i in range(len(A[0])):
            A[i][i]=b**2*(n_0**2-1)/(2*k_jz(th[t],i-n_diff,G))-dq_j(th[t],i-n_diff,G)
            if(i+1<len(A[0])):
                A[i][i+1]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G))
                A[i+1][i]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G))
        A=-1j*A
        # print(A)
        w,v = np.linalg.eig(A)
        # for i in range(len(w)):
        #     print("AV=wV", np.allclose(np.dot(A,v[:,i]),w[i]*v[:,i]))
        # print(, )
        # print(w,v)
        v0=np.zeros(2*n_diff+1)
        v0[n_diff]=1
        #print(v0)
        c = np.linalg.solve(v,v0)
        for i in range(len(w)):
            v[:,i]=v[:,i]*c[i]*np.exp(w[i]*d)
        for i in range(len(S[:,0])):
            S[i,t] = sum(v[i,:])
        #print(np.allclose(np.dot(v, c), v0))
    eta = S.copy().real
    for t in range(len(th)):
        for i in range(2*n_diff+1):
            eta[i,t] = abs(S[i,t])**2*k_jz(th[t],i-n_diff,G)/(b*np.cos(th[t]))
        sum_diff[t]= sum(eta[:,t])
    # print(sum_diff)
    # for i in range(len(S[0,:])): 
    #         S[:,i]=(abs(S[:,i])/sum(abs(S[:,i])**2)**0.5)
    for i in range(len(diff_eff[:,0])):
            s=sum(diff_eff[i,1:])
            diff_eff[i,1:]=diff_eff[i,1:]/s
    # fig = plt.figure(figsize=(15,15))
    # ax = fig.add_subplot(111)
    fig, ax = plt.subplots(n_diff+2,figsize=(10,10))
    ax[0].set_title(foldername[k]) 
    ax[0].plot(th,eta[n_diff,:])  
    ax[0].plot(diff_eff[:,0]*rad,diff_eff[:,2*2+2],'o')
    for i in range(1,n_diff+1):
        ax[i].plot(th,eta[n_diff-i,:])
        ax[i].plot(th,eta[n_diff+i,:])   
        if i<3:
            ax[i].plot(diff_eff[:,0]*rad,diff_eff[:,6-2*i],'o')
            ax[i].plot(diff_eff[:,0]*rad,diff_eff[:,6+2*i],'o')
    ax[n_diff+1].plot(th, sum_diff)
    ax[n_diff+1].set_ylim([0.5,1.5])
    #   plt.errorbar(diff_eff[:,0],diff_eff[:,2*j+2],yerr=diff_eff[:,2*j+1],capsize=1)