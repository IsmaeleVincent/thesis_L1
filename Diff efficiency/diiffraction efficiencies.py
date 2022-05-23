#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 15:07:46 2022

@author: aaa
"""
import os
import shutil
import numpy as np
from PIL import Image as im
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
import math
from scipy.interpolate import interp1d

sorted_fold_path="/home/exp-k03/Desktop/thesis/Sorted data/" #insert folder of sorted meausements files
allmeasurements = sorted_fold_path+"All measurements/"
allrenamed = allmeasurements +"All renamed/"
allmatrixes = allmeasurements + "All matrixes/"
allpictures = allmeasurements + "All pictures/"
allrawpictures = allmeasurements + "All raw pictures/"
alldata_analysis = allmeasurements + "All Data Analysis/"
allcropped_pictures = alldata_analysis + "All Cropped Pictures/"
allcontrolplots = alldata_analysis + "All Control plots/"
allcorrectedfits = alldata_analysis + "All Corrected Fits/"
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg")
foldername.append("79-4U_77c88deg")
foldername.append("79-8U_76c76deg")
foldername.append("79-12U_75c64deg")
foldername.append("79-16U_74c52deg")
n_theta=[26,46,28,17,16,20,21,20,19,48,43,59,24]  #number of measurements files for each folder (no flat, no compromised data)
step_theta=[0.03,0.03,0.05,0.05,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03]
n_pixel = 16384 #number of pixels in one measurement


"""
This block calculates the diffraction efficiencies
"""
plot=0
def gauss(x, A, x0,sx):
      return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))

for k in range(len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
    stack = np.stack(matrixes,axis=2)
    xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
    yabsmax = xyzabsmax[0][0]
    xabsmax = xyzabsmax[1][0]
    zabsmax = xyzabsmax[2][0]
    roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
    data_and_fit  =  np.loadtxt(data_analysis+foldername[k]+'_fit+data.mpa',skiprows=1)
    diff_eff = np.zeros((len(stack[0,0,:]),12))
    #print(foldername[k])
    for z in range(len(stack[0,0,:])):
        zprofile0 = np.zeros(len(stack[0,0,:]))
        zprofile0 += stack[yabsmax,xabsmax,:]
        # for i in range(3):
        #     for j in range(3):
        #         zprofile0 += stack[yabsmax+i-1,xabsmax+j-1,:].copy()/6
    zmin1=roi[:,7][roi[:,0]==yabsmax]
    zmin2=roi[:,8][roi[:,0]==yabsmax]
    f2 = interp1d(np.where(zprofile0)[0], zprofile0, kind='cubic')
    zplt=np.linspace(0,len(zprofile0)-1, 10000)
    # if (k==6):
    #     zplt=np.linspace(4,len(zprofile0)-1, 10000)
    zplt1=np.linspace(zmin1,zmin2, 1000)
    zmax= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
    if(zmin1>0 and zmin2<len(stack[0,0,:])-1):
        zplt1=np.linspace(0,zmax, 1000)
        z1=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        zplt1=np.linspace(zmax,len(stack[0,0,:])-1, 1000)
        z2=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        c=(z1+z2)*0.5
    else:
        c= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
        z1=zmin1
        z2=zmin2
    if (k==0):
        c=19
    if (k==4):
        c=12
    if (k==6):
        c=19.1
    if (k==7):
        c=16
    if (k==8):
        c=18
    if (k==12):
        c=23
    if(plot):    
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111)
        ax.set_title(foldername[k])
        ax.axvline(zmax, color="b")
        ax.plot(np.where(zprofile0)[0],zprofile0, "ko")
        ax.plot(zplt,f2(zplt), "b-")
        ax.axvline(z1, color="r")
        ax.axvline(z2, color="g")
        ax.axvline(c, color="k")
    P0m = np.zeros(9)
    P0p = np.zeros(9)
    #print(foldername[k])
    for y in range(len(roi[:,0])):
        for z in range(len(stack[0,0,:])):
            diff_eff[z][0]=(z-c)*step_theta[k]
            if data_and_fit[z*len(roi[:,0])+y][1]>0:
                bckg = 0#data_and_fit[z*len(roi[:,0])+y][2]
                for j in range(len(P0m)):
                    P0m[j]=data_and_fit[z*len(roi[:,0])+y][j+3] 
                    P0p[j]=data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)]
                data=np.zeros((roi[y][2]-roi[y][1]+1,2))
                data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
                data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
                data[:,0] = (data[:,0]-xabsmax)
                xplt=data[:, 0]
                color=["r-","g-","k-"]
                for j in range(3):
                    if (P0m[2-j]>0):
                        diff_eff[z][2+j*2]+=sum(gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))+bckg*len(xplt)
                        diff_eff[z][3+j*2]= diff_eff[z][2+j*2]**0.5 #sum((bckg+gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))**0.5)
                        if (j>0 and P0p[j]>0):
                            diff_eff[z][6+j*2]+=sum(bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))
                            diff_eff[z][7+j*2]=diff_eff[z][6+j*2]**0.5#sum((bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))**0.5)
            if(k==6 and (z==3 or z==4)):
                diff_eff[z][:]*=0
    # if (k==6):
    #     diff_eff=np.delete(diff_eff, [3,4])
    with open(data_analysis+foldername[k]+'_diff_eff.mpa', 'w') as f:
        np.savetxt(f,diff_eff, header="theta err counts-2 err counts-1 err counts-0 err counts1 err counts1 err", fmt="%.6f")

"""
This block calculates the diffraction efficiencies for 3 lines
"""
plot=0
def gauss(x, A, x0,sx):
      return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))

for k in range(len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
    stack = np.stack(matrixes,axis=2)
    xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
    yabsmax = xyzabsmax[0][0]
    xabsmax = xyzabsmax[1][0]
    zabsmax = xyzabsmax[2][0]
    roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
    data_and_fit  =  np.loadtxt(data_analysis+foldername[k]+'_fit+data.mpa',skiprows=1)
    diff_eff = np.zeros((len(stack[0,0,:]),12))
    #print(foldername[k])
    for z in range(len(stack[0,0,:])):
        zprofile0 = np.zeros(len(stack[0,0,:]))
        zprofile0 += stack[yabsmax,xabsmax,:]
        # for i in range(3):
        #     for j in range(3):
        #         zprofile0 += stack[yabsmax+i-1,xabsmax+j-1,:].copy()/6
    zmin1=roi[:,7][roi[:,0]==yabsmax]
    zmin2=roi[:,8][roi[:,0]==yabsmax]
    f2 = interp1d(np.where(zprofile0)[0], zprofile0, kind='cubic')
    zplt=np.linspace(0,len(zprofile0)-1, 10000)
    # if (k==6):
    #     zplt=np.linspace(4,len(zprofile0)-1, 10000)
    zplt1=np.linspace(zmin1,zmin2, 1000)
    zmax= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
    if(zmin1>0 and zmin2<len(stack[0,0,:])-1):
        zplt1=np.linspace(0,zmax, 1000)
        z1=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        zplt1=np.linspace(zmax,len(stack[0,0,:])-1, 1000)
        z2=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        c=(z1+z2)*0.5
    else:
        c= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
        z1=zmin1
        z2=zmin2
    if (k==0):
        c=19
    if (k==4):
        c=12
    if (k==6):
        c=19.1
    if (k==7):
        c=16
    if (k==8):
        c=18
    if (k==12):
        c=23
    if(plot):    
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111)
        ax.set_title(foldername[k])
        ax.axvline(zmax, color="b")
        ax.plot(np.where(zprofile0)[0],zprofile0, "ko")
        ax.plot(zplt,f2(zplt), "b-")
        ax.axvline(z1, color="r")
        ax.axvline(z2, color="g")
        ax.axvline(c, color="k")
    P0m = np.zeros(9)
    P0p = np.zeros(9)
    #print(foldername[k])
    y0= np.where(roi[:,0]==yabsmax)[0][0]
    print("yabsmax=",yabsmax)
    for y in range(y0-1,y0+2):#len(roi[:,0])):
        print(roi[y,0])
        for z in range(len(stack[0,0,:])):
            diff_eff[z][0]=(z-c)*step_theta[k]
            if data_and_fit[z*len(roi[:,0])+y][1]>0:
                bckg = 0#data_and_fit[z*len(roi[:,0])+y][2]
                for j in range(len(P0m)):
                    P0m[j]=data_and_fit[z*len(roi[:,0])+y][j+3] 
                    P0p[j]=data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)]
                data=np.zeros((roi[y][2]-roi[y][1]+1,2))
                data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
                data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
                data[:,0] = (data[:,0]-xabsmax)
                xplt=data[:, 0]
                color=["r-","g-","k-"]
                for j in range(3):
                    if (P0m[2-j]>0):
                        diff_eff[z][2+j*2]+=sum(gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))+bckg*len(xplt)
                        diff_eff[z][3+j*2]= diff_eff[z][2+j*2]**0.5 #sum((bckg+gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))**0.5)
                        if (j>0 and P0p[j]>0):
                            diff_eff[z][6+j*2]+=sum(bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))
                            diff_eff[z][7+j*2]=diff_eff[z][6+j*2]**0.5#sum((bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))**0.5)
            if(k==6 and (z==3 or z==4)):
                diff_eff[z][:]*=0
    # if (k==6):
    #     diff_eff=np.delete(diff_eff, [3,4])
    with open(data_analysis+foldername[k]+'_diff_eff_3lines.mpa', 'w') as f:
        np.savetxt(f,diff_eff, header="theta err counts-2 err counts-1 err counts-0 err counts1 err counts1 err", fmt="%.6f")

"""
This block calculates the diffraction efficiencies for 1 lines
"""
plot=0
def gauss(x, A, x0,sx):
      return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))

for k in range(len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
    stack = np.stack(matrixes,axis=2)
    xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
    yabsmax = xyzabsmax[0][0]
    xabsmax = xyzabsmax[1][0]
    zabsmax = xyzabsmax[2][0]
    roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
    data_and_fit  =  np.loadtxt(data_analysis+foldername[k]+'_fit+data.mpa',skiprows=1)
    diff_eff = np.zeros((len(stack[0,0,:]),12))
    #print(foldername[k])
    for z in range(len(stack[0,0,:])):
        zprofile0 = np.zeros(len(stack[0,0,:]))
        zprofile0 += stack[yabsmax,xabsmax,:]
        # for i in range(3):
        #     for j in range(3):
        #         zprofile0 += stack[yabsmax+i-1,xabsmax+j-1,:].copy()/6
    zmin1=roi[:,7][roi[:,0]==yabsmax]
    zmin2=roi[:,8][roi[:,0]==yabsmax]
    f2 = interp1d(np.where(zprofile0)[0], zprofile0, kind='cubic')
    zplt=np.linspace(0,len(zprofile0)-1, 10000)
    # if (k==6):
    #     zplt=np.linspace(4,len(zprofile0)-1, 10000)
    zplt1=np.linspace(zmin1,zmin2, 1000)
    zmax= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
    if(zmin1>0 and zmin2<len(stack[0,0,:])-1):
        zplt1=np.linspace(0,zmax, 1000)
        z1=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        zplt1=np.linspace(zmax,len(stack[0,0,:])-1, 1000)
        z2=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
        c=(z1+z2)*0.5
    else:
        c= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
        z1=zmin1
        z2=zmin2
    if (k==0):
        c=19
    if (k==4):
        c=12
    if (k==6):
        c=19.1
    if (k==7):
        c=16
    if (k==8):
        c=18
    if (k==12):
        c=23
    if(plot):    
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111)
        ax.set_title(foldername[k])
        ax.axvline(zmax, color="b")
        ax.plot(np.where(zprofile0)[0],zprofile0, "ko")
        ax.plot(zplt,f2(zplt), "b-")
        ax.axvline(z1, color="r")
        ax.axvline(z2, color="g")
        ax.axvline(c, color="k")
    P0m = np.zeros(9)
    P0p = np.zeros(9)
    #print(foldername[k])
    print("yabsmax=",yabsmax)
    y0= np.where(roi[:,0]==yabsmax)[0][0]
    for y in range(y0,y0+1):#len(roi[:,0])):
        print(roi[y,0])
        for z in range(len(stack[0,0,:])):
            diff_eff[z][0]=(z-c)*step_theta[k]
            if data_and_fit[z*len(roi[:,0])+y][1]>0:
                bckg = 0#data_and_fit[z*len(roi[:,0])+y][2]
                for j in range(len(P0m)):
                    P0m[j]=data_and_fit[z*len(roi[:,0])+y][j+3] 
                    P0p[j]=data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)]
                data=np.zeros((roi[y][2]-roi[y][1]+1,2))
                data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
                data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
                data[:,0] = (data[:,0]-xabsmax)
                xplt=data[:, 0]
                color=["r-","g-","k-"]
                for j in range(3):
                    if (P0m[2-j]>0):
                        diff_eff[z][2+j*2]+=sum(gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))+bckg*len(xplt)
                        diff_eff[z][3+j*2]= diff_eff[z][2+j*2]**0.5 #sum((bckg+gauss(xplt, P0m[2-j], -P0m[2-j+3], P0m[2-j+6]))**0.5)
                        if (j>0 and P0p[j]>0):
                            diff_eff[z][6+j*2]+=sum(bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))
                            diff_eff[z][7+j*2]=diff_eff[z][6+j*2]**0.5#sum((bckg+gauss(xplt, P0p[j], P0p[j+3], P0p[j+6]))**0.5)
            if(k==6 and (z==3 or z==4)):
                diff_eff[z][:]*=0
    # if (k==6):
    #     diff_eff=np.delete(diff_eff, [3,4])
    with open(data_analysis+foldername[k]+'_diff_eff_1line.mpa', 'w') as f:
        np.savetxt(f,diff_eff, header="theta err counts-2 err counts-1 err counts-0 err counts1 err counts1 err", fmt="%.6f")



"""
This shows something inteesting but I'm not sure what to do with it yet
"""

# for k in range(len(foldername)):
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#     yabsmax = xyzabsmax[0][0]
#     xabsmax = xyzabsmax[1][0]
#     zabsmax = xyzabsmax[2][0]
#     roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
#     print(foldername[k])
#     fig, ax = plt.subplots(3,figsize=(10,10))
#     for j in range(3):
#         for z in range(len(stack[0,0,:])):
#             zprofile0 = np.zeros(len(stack[0,0,:]))
#             zprofile0 += stack[yabsmax+j,xabsmax,:].copy()
#             # for i in range(3):
#             #     for j in range(3):
#             #         zprofile0 += stack[yabsmax+i-1,xabsmax+j-1,:].copy()/6
#         zmin1=roi[:,7][roi[:,0]==yabsmax]
#         zmin2=roi[:,8][roi[:,0]==yabsmax]
#         print(zmin1, zmin2)
#         #fig = plt.figure(figsize=(15,15))
#         #ax = fig.add_subplot(111)
#         #ax.set_title(foldername[k])
#         f2 = interp1d(np.where(zprofile0)[0], zprofile0, kind='cubic')
#         zplt=np.linspace(0,len(zprofile0)-1, 10000)
#         if (k==6):
#             zplt=np.linspace(4,len(zprofile0)-1, 10000)
#         zplt1=np.linspace(zmin1,zmin2, 1000)
#         zmax= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
#         if(zmin1>0 and zmin2<len(stack[0,0,:])-1):
#             zplt1=np.linspace(0,zmax, 1000)
#             z1=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
#             zplt1=np.linspace(zmax,len(stack[0,0,:])-1, 1000)
#             z2=zplt1[f2(zplt1)==np.amin(f2(zplt1))]
#             c=(z1+z2)*0.5
#         else:
#             c= zplt1[f2(zplt1)==np.amax(f2(zplt1))]
#             z1=zmin1
#             z2=zmin2
#         ax[j].axvline(zmax, color="b")
#         ax[j].plot(np.where(zprofile0)[0],zprofile0, "ko")
#         ax[j].plot(zplt,f2(zplt), "b-")
#         ax[j].axvline(z1, color="r")
#         ax[j].axvline(z2, color="g")
#         ax[j].axvline(c, color="k")
"""
# This block plots the diffraction efficiencies
"""

for k in range(len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff.mpa',skiprows=1)
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)
    ax.set_title(foldername[k])
    for j in range(5):
        ax.plot(diff_eff[:,0],diff_eff[:,2*j+2],'o')
        ax.errorbar(diff_eff[:,0],diff_eff[:,2*j+2],yerr=diff_eff[:,2*j+1],capsize=1)
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff_3lines.mpa',skiprows=1)
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)
    ax.set_title(foldername[k]+"_diff_eff_3lines")
    for j in range(5):
        ax.plot(diff_eff[:,0],diff_eff[:,2*j+2],'o')
        ax.errorbar(diff_eff[:,0],diff_eff[:,2*j+2],yerr=diff_eff[:,2*j+1],capsize=1)
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff_1line.mpa',skiprows=1)
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)
    ax.set_title(foldername[k]+"_diff_eff_1line")
    for j in range(5):
        ax.plot(diff_eff[:,0],diff_eff[:,2*j+2],'o')
        ax.errorbar(diff_eff[:,0],diff_eff[:,2*j+2],yerr=diff_eff[:,2*j+1],capsize=1)


