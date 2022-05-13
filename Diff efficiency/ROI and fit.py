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
n_theta=[26,46,28,18,16,20,21,20,19,48,43,59,24]  #number of measurements files for each folder (no flat, no compromised data)
n_pixel = 16384 #number of pixels in one measurement

"""
This block creates a 3D matrix (stack) with measurments along the x,y plane and theta along z, 
and uses it to calculate the region of interest (ROI)
"""
# for k in range(len(foldername)):
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     xyzmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#     ymax = xyzmax[0][0]
#     xmax = xyzmax[1][0]
#     zmax = xyzmax[2][0]
#     y1=0
#     a=0.
#     b=0.
#     s=0.
#     while y1<ymax and abs(a)<3:#(abs(a)<2 or abs(b)<2):
#         s+=stack[y1,xmax,0]
#         avg=s/(y1+1)
#         if(y1>65):
#             a=stack[y1+3,xmax,0]/avg
#             b=stack[y1+4,xmax,0]/avg
#         y1+=1
#     #y1-=1
#     y2=127
#     a=0.
#     b=0.
#     s=0.
#     while y2>ymax and abs(a)<3:#(abs(a)<2 or abs(b)<2):
#         s+=stack[y2,xmax,0]
#         avg=s/(128-y2)
#         if(128-y2>48):
#             a=stack[y2-3,xmax,0]/avg
#             b=stack[y1-4,xmax,0]/avg
#         y2-=1
#     #y2+=1
#     x1=np.zeros(y2-y1,dtype=int)
#     x2=np.zeros(y2-y1,dtype=int)+127
#     for i in range(y2-y1):    
#         a=0.
#         b=0.
#         s=0.
#         while x1[i]<xmax and abs(a)<3:#(abs(a)<2 or abs(b)<2):
#             s+= np.amax(stack[y1+i,x1[i],:])
#             avg=s/(x1[i]+1)
#             if(x1[i]>20):
#                 a=np.amax(stack[y1+i,x1[i]+3,:])/avg
#                 b=np.amax(stack[y1+i,x1[i]+4,:])/avg
#             x1[i]+=1
#         #x1[i]-=1
#         a=0.
#         b=0.
#         s=0.
#         while x2[i]>xmax and abs(a)<3:#(abs(a)<2 or abs(b)<2):
#             s+=np.amax(stack[y1+i,x2[i],:])
#             avg=s/(128-x2[i])
#             if(128-x2[i]>20):
#                 a=np.amax(stack[y1+i,x2[i]-3,:])/avg
#                 b=np.amax(stack[y1+i,x1[i]-4,:])/avg
#             x2[i]-=1
#         #x2[i]+=1
#         roi= np.zeros((y2-y1,3))
#         for i in range(y2-y1):
#             roi[i][0] = y1+i
#             roi[i][1]= x1[i]
#             roi[i][2]= x2[i]
#     with open(data_analysis+foldername[k]+'_ROI.mpa', 'w') as f:
#                 np.savetxt(f,roi, header="line x1 x2", fmt="%1.0f")

"""
This block checks the number of peaks of each line and puts it in the ROI file
the plots the max points along z for each line to check if the ROI is reasonable
"""
# for k in range(len(foldername)):
#     print(foldername[k])
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     controlplots = data_analysis + "Control Plots/" 
#     if os.path.exists(controlplots):
#         shutil.rmtree(controlplots)
#     os.makedirs(controlplots)
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI.mpa',skiprows=1).astype(int)    
#     y=np.zeros(128)
#     npeaks=np.zeros((len(roi[:,0]),2)).astype(int)
#     xpeaks=np.zeros((len(roi[:,0]),2)).astype(int)
#     zmins=np.zeros((len(roi[:,0]),2)).astype(int)
#     zmin1=0
#     zmin2=0
#     for i in range(len(roi[:,0])):
#         for j in range(128):
#             y[j]=np.amax(stack[roi[i][0],j,:])
#         xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#         xabsmax = xyzabsmax[1][0]
#         if(roi[i][1]==roi[i][2]):
#             xmax=roi[i][1]
#         else:
#             xmax= np.where(y[roi[i][1]:(roi[i][2]+1)]==np.max(y[roi[i][1]:(roi[i][2]+1)]))[0][0]+roi[i][1]
#         if(abs(xmax-roi[i][1])<5):
#             npeaks[i][0]=0
#             a=roi[i][1]
#         else:
#             a=np.where(y[roi[i][1]:(xmax-1)]==np.amax(y[roi[i][1]:(xmax-1)]))[0][0] + roi[i][1]
#             if(abs((roi[i][1])-a)<10):
#                 npeaks[i][0]=1
#             else:
#                 npeaks[i][0]=2
#         xpeaks[i][0]=a
#         if(abs(roi[i][2]-xmax)<5):
#             npeaks[i][1]=0
#             b=roi[i][2]
#         else:
#             b=np.where(y[(xmax+2):(roi[i][2]+1)]==np.amax(y[(xmax+2):(roi[i][2]+1)]))[0][0]+xmax+2
#             if(abs((roi[i][2])-b)<10):
#                 npeaks[i][1]=1
#             else:
#                 npeaks[i][1]=2
#         xpeaks[i][1]=b
#     for i in range(len(roi[:,0])):
#         for j in range(128):
#             y[j]=np.amax(stack[roi[i][0],j,:])
#         xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#         xabsmax = xyzabsmax[1][0]
#         fig, axs = plt.subplots(4,figsize=(10,10))
#         fig.suptitle(foldername[k] +'-Line ' +str("%0d"%(roi[0][0]+i))+"\nleft peaks="+str(npeaks[i][0])+"\nright peaks="+str(npeaks[i][1]))
#         axs[0].plot(y)
#         axs[0].set_xlim(0,roi[i][1]+3)
#         axs[0].set_ylim(0,y[int(roi[i][1]+3)]+5)
#         axs[0].axvline(roi[i][1],color='r')
#         axs[1].plot(y)
#         axs[1].set_xlim(roi[i][2]-3,128)
#         axs[1].set_ylim(0,y[int(roi[i][2])]+5)
#         axs[1].axvline(roi[i][2],color='r')
#         axs[2].plot(y)
#         axs[2].axvline(roi[i][1],color='r')
#         axs[2].axvline(roi[i][2],color='r')
#         axs[2].axvline(xpeaks[i][0],color='y', ls = '--')
#         axs[2].axvline(xpeaks[i][1],color='k', ls = '--')
#         if (npeaks[i][0]>0 or npeaks[i][1]>0):
#             aus = (stack[roi[i][0],xabsmax,:].copy()+stack[roi[i][0],xabsmax-1,:].copy()+stack[roi[i][0],xabsmax+1,:].copy())/3
#             if (np.amax(aus)>150):
#                 zmin1 = np.where(aus==np.amin(aus))[0][0]
#                 if (np.amax(npeaks[:,1])<2):
#                     zmin1 = len(aus)-1
#                 for l in range(5):
#                     if (zmin1+l<len(aus)):
#                         aus[zmin1+l] += 1000
#                     if (zmin1-l>-1):
#                         aus[zmin1-l] += 1000
#                 zmin2=np.where(aus==np.amin(aus))[0][0]
#                 if (zmin1>zmin2):
#                     zmin1,zmin2 = zmin2,zmin1
#             else:
#                 if (zmin1>0):
#                     zmin1+=-1
#                 if (zmin2<len(aus)-1):
#                     zmin2+=1
#         zmins[i][0]=zmin1
#         zmins[i][1]=zmin2
#         axs[3].plot((stack[roi[i][0],xabsmax,:].copy()+stack[roi[i][0],xabsmax-1,:].copy()+stack[roi[i][0],xabsmax+1,:].copy())/3)
#         axs[3].axvline(zmins[i][0],color='g')
#         axs[3].axvline(zmins[i][1],color='k')
#         plt.savefig(controlplots+foldername[k] +'_line' +str("%0d"%(roi[0][0]+i))+'.png')
#         plt.close(fig)
#     roipeaks= np.concatenate((roi,npeaks,xpeaks,zmins),axis=1)
#     with open(data_analysis+foldername[k]+'_ROI+Peaks.mpa', 'w') as f:
#         np.savetxt(f,roipeaks, header="line x1 x2 left_peaks right_peaks xmax-1 xmax+1 zmin-1 zmin+1", fmt="%1.0f")

"""
This block creates cropped pictures of the ROI, to check that everything is ok
"""
# for k in range(len(foldername)):
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     cropped_pictures = data_analysis+"Cropped Pictures/"
#     if os.path.exists(cropped_pictures):
#         shutil.rmtree(cropped_pictures)
#     os.makedirs(cropped_pictures)
#     for j in range (1,n_theta[k]+1):
#         f = sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa"
#         data_matrix = np.loadtxt(f)
#         roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI.mpa',skiprows=1)
#         ymin= int(np.amin(roi[:,0]))
#         ymax= int(np.amax(roi[:,0])) 
#         xmin= int(np.amin(roi[:,1]))
#         xmax= int(np.amax(roi[:,2]))
#         image = im.fromarray((data_matrix[ymin:ymax,xmin:xmax]))
#         image.save(cropped_pictures+foldername[k]+"_cropped_"+str("%03d" % (j,))+'.tiff')
"""
This block creates a folder with the plot of the max points along z for each 
line to check if the ROI is reasonable
"""
# for k in range(len(foldername)):
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     controlplots = data_analysis + "Control Plots/" 
#     if os.path.exists(controlplots):
#         shutil.rmtree(controlplots)
#     os.makedirs(controlplots)
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI.mpa',skiprows=1)    
#     y=np.zeros(128)
#     for i in range(len(roi[:,0])):
#         for j in range(128):
#             y[j]=np.amax(stack[int(roi[0][0])+i,j,:])
#         fig, axs = plt.subplots(3,figsize=(10,10))
#         fig.suptitle(foldername[k] +'-Line ' +str("%0d"%(roi[0][0]+i)))
#         axs[0].plot(y)
#         axs[0].set_xlim(0,roi[i][1]+3)
#         axs[0].set_ylim(0,y[int(roi[i][1]+3)]+5)
#         axs[0].axvline(roi[i][1],color='r')
#         axs[1].plot(y)
#         axs[1].set_xlim(roi[i][2]-3,128)
#         axs[1].set_ylim(0,y[int(roi[i][2])]+5)
#         axs[1].axvline(roi[i][2],color='r')
#         axs[2].plot(y)
#         axs[2].axvline(roi[i][1],color='r')
#         axs[2].axvline(roi[i][2],color='r')
#         plt.savefig(controlplots+foldername[k] +'_line' +str("%0d"%(roi[0][0]+i))+'.png')
#         plt.close(fig)
"""
This block copies the files in a common folder just 
for the sake of simplicity
"""
# if os.path.exists(alldata_analysis):
#     shutil.rmtree(alldata_analysis)
# os.makedirs(alldata_analysis)
# if os.path.exists(allcropped_pictures):
#     shutil.rmtree(allcropped_pictures)
# os.makedirs(allcropped_pictures)
# if os.path.exists(allcontrolplots):
#     shutil.rmtree(allcontrolplots)
# os.makedirs(allcontrolplots)
# for k in range(len(foldername)):
#     folder = sorted_fold_path+foldername[k]
#     for j in range (1,n_theta[k]+1):
#         croppicname = foldername[k]+"_cropped_"+str("%03d" % (j,))+".tiff"
#         shutil.copy(folder+"/Data Analysis/Cropped Pictures/"+croppicname, allcropped_pictures+croppicname)
#     roi =  np.loadtxt(folder+"/Data Analysis/"+foldername[k]+'_ROI.mpa',skiprows=1)    
#     for j in range(len(roi[:,0])):
#         contplotname = foldername[k] +'_line' +str("%0d"%(roi[0][0]+j))+'.png'
#         shutil.copy(folder+"/Data Analysis/Control Plots/"+contplotname, allcontrolplots+contplotname)
"""
This block calculates the fits
"""

# plot=0
# def gauss(x, A, x0,sx):
#     return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))
# def distrm3(x,A0,A1,A2,x0,x1,x2,s0,s1,s2):
#     return bg1+gauss(x,A0,-x0,s0)+gauss(x, A1,-x1,s1)+gauss(x, A2, -x2,s2)
# def distrp3(x,A0,A1,A2,x0,x1,x2,s0,s1,s2):
#     return bg1+gauss(x,A0,x0,s0)+gauss(x, A1, x1,s1)+gauss(x, A2, x2,s2)
# def distr1(x, A, x0,sx):
#     return bg1+A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))

# for k in range(0,1):#len(foldername)):
#     bckg=0.
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#     yabsmax = xyzabsmax[0][0]
#     xabsmax = xyzabsmax[1][0]
#     zabsmax = xyzabsmax[2][0]
#     roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
#     print(foldername[k])
#     data_and_fit = np.zeros((len(stack[0,0,:])*len(roi[:,0]), 21))
#     for y in range(len(roi[:,0])):
#         if (roi[y][2] == roi[y][1]):
#             bckg = sum(stack[roi[y][0],xabsmax-20:xabsmax-3,0])/len(stack[roi[y][0],xabsmax-20:xabsmax-3,0])
#         if ((roi[y][3]==0 or roi[y][4]==0) and roi[y][2]-roi[y][1]>5):
#             P0 = [1.,0.,0.5]
#             bound = [[0.2,-3.,0.01,],[1.5,3.,1.]]
#             for z in range(len(stack[0,0,:])):
#                 data=np.zeros((roi[y][2]-roi[y][1]+1,2))
#                 data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
#                 data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
#                 ymax = np.amax(data[:,1])
#                 ymin = np.amin(data[:,1])
#                 bg1 = bckg/ymax
#                 data[:,0] = (data[:,0]-xabsmax)
#                 data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
#                 p,cov=fit(distr1,data[:,0][data[:,0]>-3],data[:,1][data[:,0]>-3], p0=P0, bounds = bound)
#                 P0=p.copy()
#                 if plot:
#                     fig = plt.figure(figsize=(15,15))
#                     ax = fig.add_subplot(111)
#                     ax.set_title(foldername[k] +'-Line ' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z)))
#                     ax.plot(data[:,0],data[:,1], "k-")
#                     xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
#                     ax.plot(xplt, distr1(xplt,*p), "r-")
#                 data_and_fit[z*len(roi[:,0])+y][0] = z
#                 data_and_fit[z*len(roi[:,0])+y][1] = roi[y][0]
#                 data_and_fit[z*len(roi[:,0])+y][2] = bckg + ymin
#                 for j in range(len(P0)):
#                     data_and_fit[z*len(roi[:,0])+y][j*3+3] = P0[j]
#                 data_and_fit[z*len(roi[:,0])+y][3]*=(ymax-ymin)
#         if (roi[y][3]>1 or roi[y][4]>1):
#             if(roi[y][3]==1):
#                 P0m = [1.,0.3,0.1, 0., abs(roi[y][5]-xabsmax)-2, abs(roi[y][5]-xabsmax)+8, 0.5, 0.5, 0.5]
#                 boundm = [[0.2,0.,0.,-2.,abs(roi[y][5]-xabsmax)-4, abs(roi[y][5]-xabsmax)+2,0.01,0.1,0.1],[1.5,1.5,0.4,2.,abs(roi[y][1]-xabsmax),abs(roi[y][1]-xabsmax), 1.,1.5,1.5]]
#             if(roi[y][3]>1):    
#                 P0m = [1.,0.3,0.1, 0., abs(roi[y][5]-xabsmax)-2, abs(roi[y][5]-xabsmax)+8, 0.5, 0.5, 0.5]
#                 boundm = [[0.2,0.,0.,-2.,abs(roi[y][5]-xabsmax)-4, abs(roi[y][5]-xabsmax)+5,0.01,0.1,0.1],[1.5,1.5,0.4,2.,abs(roi[y][1]-xabsmax),abs(roi[y][1]-xabsmax), 1.,1.5,1.5]]
#             if(roi[y][4]==1):
#                 P0p = [1.,0.,0., 0., abs(roi[y][6]-xabsmax)-2, abs(roi[y][6]-xabsmax)+3, 0.5, 1., 1.5]
#                 boundp = [[0.2,0,0,-2.,abs(roi[y][6]-xabsmax)-4, abs(roi[y][6]-xabsmax)+1,0.01,0.1,0.1],[2.,0.4,1e-8,2.,abs(roi[y][2]-xabsmax),abs(roi[y][2]-xabsmax), 1,1.5,1.5]]
#             if(roi[y][4]>1):    
#                 P0p = [1.,0.,0., 0., abs(roi[y][6]-xabsmax)-2, abs(roi[y][6]-xabsmax)+3, 0.5, 1., 1.5]
#                 boundp = [[0.2,0,0,-2.,abs(roi[y][6]-xabsmax)-4, abs(roi[y][6]-xabsmax)+2,0.01,0.1,0.1],[2.,0.4,1e-8,2.,abs(roi[y][2]-xabsmax),abs(roi[y][2]-xabsmax), 1,1.5,1.5]]
#             P0maus=P0m.copy()
#             P0paus=P0p.copy()
#             boundmaus=boundm.copy()
#             boundpaus=boundp.copy()
#             zmin1=roi[y][7]
#             zmin2=roi[y][8]
#             for z in range(len(stack[0,0,:])):
#                 data=np.zeros((roi[y][2]-roi[y][1]+1,2))
#                 data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
#                 data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
#                 ymax = np.amax(data[:,1])
#                 ymin = np.amin(data[:,1])
#                 bg1 = bckg/ymax
#                 data[:,0] = (data[:,0]-xabsmax)
#                 data[:,1] = (data[:,1]-ymin)/(ymax-ymin)
#                 if plot:
#                     fig = plt.figure(figsize=(15,15))
#                     ax = fig.add_subplot(111)
#                     ax.set_title(foldername[k] +'-Line ' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z)))
#                     ax.plot(data[:,0],data[:,1], "ko")
#                 if(roi[y][3]==0):
#                     boundm[1][1] = 1e-8
#                     P0m[1]=0.
#                     boundm[1][2] = 1e-8
#                     P0m[2]=0.
#                 if(z>zmin1-2):
#                     if(z>zmin1+2):
#                         boundm[1][2] = 1e-8
#                         P0m[2]=0.
#                         boundp[1][1]= 1.5
#                     else:
#                         boundm[1][2] = 0.1
#                         boundm[1][8] = 1.
#                         P0m[2]=0.09
#                 if(roi[y][4]==0):
#                     boundp[1][1] = 1e-8
#                     P0p[1]=0.
#                     boundp[1][2] = 1e-8
#                     P0p[2]=0.
#                 if(z>zmin2-4):
#                     if(z>zmin2+2):
#                         boundp[1][2] = 0.3
#                         boundp[1][8] = 1.5
#                     else:
#                         boundp[1][2] = 0.1
#                         boundp[1][8] = 1.
#                 for j in range(len(P0p)):
#                     if (P0m[j]<=boundm[0][j] or P0m[j]>=boundm[1][j]):
#                         P0m[j]=(boundm[1][j]-boundm[0][j])/3+boundm[0][j]
#                     if (P0p[j]<=boundp[0][j] or P0p[j]>=boundp[1][j]):
#                         P0p[j]=(boundp[1][j]-boundp[0][j])/3+boundp[0][j]
#                     # if (boundm[0][j]>=boundm[1][j]):
#                     #     print(j,"m",boundm[0][j],boundm[1][j])
#                     #     boundm[0][j]=boundmaus[0][j]
#                     #     boundm[1][j]=boundmaus[1][j]
#                     # if (boundp[0][j]>=boundp[1][j]):
#                     #     print(j,"p",boundp[0][j],boundp[1][j])
#                     #     boundp[0][j]=boundpaus[0][j]
#                     #     boundp[1][j]=boundpaus[1][j]
#                 #print(P0m, P0p)
#                 # for j in range(len(P0p)):
#                 #     if (P0p[j]<boundp[0][j] or P0p[j]>boundp[1][j]):
#                 #         print(P0p[j],"p",boundp[0][j],boundp[1][j])
#                 #     if (P0m[j]<boundm[0][j] or P0m[j]>boundm[1][j]):
#                 #         print(j,P0m[j],"m",boundm[0][j],boundm[1][j])
#                 boundmt = tuple(boundm.copy())
#                 boundpt = tuple(boundp.copy())
#                 try:
#                     pm,covm=fit(distrm3,data[:,0][data[:,0]<3],data[:,1][data[:,0]<3], p0=P0m, bounds = boundmt)
#                     pp,covp=fit(distrp3,data[:,0][data[:,0]>-3],data[:,1][data[:,0]>-3], p0=P0p, bounds = boundpt)
#                 except RuntimeError or ValueError:
#                     try:
#                         P0m=P0maus.copy()
#                         P0p=P0paus.copy()
#                         boundm=boundmaus.copy()
#                         boundp=boundpaus.copy()
#                         if(roi[y][3]==0):
#                             boundm[1][1] = 1e-8
#                             P0m[1]=0.
#                             boundm[1][2] = 1e-8
#                             P0m[2]=0.
#                         if(z>zmin1-2):
#                             if(z>zmin1+2):
#                                 boundm[1][2] = 1e-8
#                                 P0m[2]=0.
#                                 boundp[1][1]= 1.5
#                             else:
#                                 boundm[1][2] = 0.1
#                                 boundm[1][8] = 1.
#                                 P0m[2]=0.09
#                         if(roi[y][4]==0):
#                             boundp[1][1] = 1e-8
#                             P0p[1]=0.
#                             boundp[1][2] = 1e-8
#                             P0p[2]=0.
#                         if(z>zmin2-4):
#                             if(z>zmin2+3):
#                                 boundp[1][2] = 0.3
#                                 boundp[1][8] = 1.5
#                             else:
#                                 boundp[1][2] = 0.1
#                                 boundp[1][8] = 1.
#                         for j in range(len(P0p)):
#                             if (P0m[j]<=boundm[0][j] or P0m[j]>=boundm[1][j]):
#                                 P0m[j]=(boundm[1][j]-boundm[0][j])/3+boundm[0][j]
#                             if (P0p[j]<=boundp[0][j] or P0p[j]>=boundp[1][j]):
#                                 P0p[j]=(boundp[1][j]-boundp[0][j])/3+boundp[0][j]
#                         boundmt = tuple(boundm.copy())
#                         boundpt = tuple(boundp.copy())
#                         pm,covm=fit(distrm3,data[:,0][data[:,0]<3],data[:,1][data[:,0]<3], p0=P0m, bounds = boundmt)
#                         pp,covp=fit(distrp3,data[:,0][data[:,0]>-3],data[:,1][data[:,0]>-3], p0=P0p, bounds = boundpt)
#                     except RuntimeError or ValueError:
#                         print('error'+foldername[k] +' Line ' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z)))
#                         print(P0m, P0p)
#                         print(boundmt)
#                         print(boundpt)
#                 data_and_fit[z*len(roi[:,0])+y][0] = z
#                 data_and_fit[z*len(roi[:,0])+y][1] = roi[y][0]
#                 data_and_fit[z*len(roi[:,0])+y][2] = bckg + ymin
#                 P0m=pm.copy()
#                 P0p=pp.copy()
#                 for j in range(len(P0m)):
#                     data_and_fit[z*len(roi[:,0])+y][j+3] = P0m[j]
#                     data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)] = P0p[j]
#                 for j in range(3):
#                     data_and_fit[z*len(roi[:,0])+y][j+3]*=(ymax-ymin)
#                     data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)]*=(ymax-ymin)
#                 if plot:
#                     xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
#                     ax.plot(xplt, distrm3(xplt,*pm), "b--")
#                     ax.plot(xplt, distrp3(xplt,*pp), "b--")
#                     color=["r-","g-","k-"]
#                     for i in range(3):
#                         ax.plot(xplt,(bg1+gauss(xplt, pm[i], -pm[i+3], pm[i+6])), color[i%3])
#                         ax.plot(xplt, (bg1+gauss(xplt, pp[i], pp[i+3], pp[i+6])), color[i%3])
#                 if (P0m[4]+4<abs(roi[y][1]-xabsmax)):
#                     boundm[0][5]=P0m[4]+4
#                 else:
#                     boundm[1][5]=P0m[4]+6
#                     boundm[0][5]=P0m[4]+4
#                 boundm[0][7]=P0m[6]
#                 boundm[0][8]=P0m[6]
#                 boundp[0][8]=P0p[6]
#                 boundp[0][7]=P0p[6]
#                 if (P0p[4]+4<abs(roi[y][2]-xabsmax)):
#                     boundp[0][5]=P0p[4]+4
#                 else:
#                     boundp[1][5]=P0p[4]+6
#                     boundp[0][5]=P0p[4]+4
#     with open(data_analysis+foldername[k]+'_fit+data.mpa', 'w') as f:
#         np.savetxt(f,data_and_fit, header="theta line bckg A0m A1m A2m x0m x1m x2m s0m s1m s2m A0p A1p A2p x0p x1p x2p s0p s1p s2p", fmt="%i %i "+"%.18e "*19)
"""
This block creates plots of the fits from the fit+data file to check them
"""
# def gauss(x, A, x0,sx):
#     return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))
# def distrm3(x,A0,A1,A2,x0,x1,x2,s0,s1,s2):
#     return gauss(x,A0,-x0,s0)+gauss(x, A1,-x1,s1)+gauss(x, A2, -x2,s2)
# def distrp3(x,A0,A1,A2,x0,x1,x2,s0,s1,s2):
#     return gauss(x,A0,x0,s0)+gauss(x, A1, x1,s1)+gauss(x, A2, x2,s2)
# def distr1(x, A, x0,sx):
#     return A/sx*np.exp(-(x-x0)**2/(2*(sx)**2))

# for k in range(0,7):#len(foldername)):
#     data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
#     controlfits = data_analysis + "Control Fits/" 
#     if os.path.exists(controlfits):
#         shutil.rmtree(controlfits)
#     os.makedirs(controlfits)
#     matrixes = [np.loadtxt(sorted_fold_path+foldername[k]+"/Matrixes/"+foldername[k]+"_"+str("%03d" % (j,))+".mpa") for j in range (1,n_theta[k]+1)]
#     stack = np.stack(matrixes,axis=2)
#     xyzabsmax = np.where(stack[:,:,:]==np.amax(stack[:,:,:]))
#     yabsmax = xyzabsmax[0][0]
#     xabsmax = xyzabsmax[1][0]
#     zabsmax = xyzabsmax[2][0]
#     roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
#     data_and_fit  =  np.loadtxt(data_analysis+foldername[k]+'_fit+data.mpa',skiprows=1)
#     P0m = np.zeros(9)
#     P0p = np.zeros(9)
#     print(foldername[k])
#     for y in range(len(roi[:,0])):
#         for z in range(len(stack[0,0,:])):
#             if data_and_fit[z*len(roi[:,0])+y][1]>0:
#                 bckg = data_and_fit[z*len(roi[:,0])+y][2]
#                 for j in range(len(P0m)):
#                     P0m[j]=data_and_fit[z*len(roi[:,0])+y][j+3] 
#                     P0p[j]=data_and_fit[z*len(roi[:,0])+y][j+3+len(P0m)]
#                 data=np.zeros((roi[y][2]-roi[y][1]+1,2))
#                 data[:,0] =  np.arange(roi[y][1],roi[y][2]+1)
#                 data[:,1] =  stack[roi[y][0],roi[y][1]:(roi[y][2]+1),z]
#                 data[:,0] = (data[:,0]-xabsmax)
                
#                 fig = plt.figure(figsize=(15,15))
#                 ax = fig.add_subplot(111)
#                 ax.set_title(foldername[k] +'-Line ' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z)))
#                 if(P0m[2]>0 or P0m[2]>0):
#                     ax.plot(data[:,0],data[:,1], "ko")
#                     xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
#                     ax.plot(xplt,bckg + distrm3(xplt,*P0m), "b--")
#                     ax.plot(xplt,bckg + distrp3(xplt,*P0p), "b--")
#                     color=["r-","g-","k-"]
#                     for i in range(3):
#                         ax.plot(xplt,(bckg+gauss(xplt, P0m[i], -P0m[i+3], P0m[i+6])), color[i%3])
#                         ax.plot(xplt, (bckg+gauss(xplt, P0p[i], P0p[i+3], P0p[i+6])), color[i%3])
#                 else:
#                     if(P0m[0]>0):
#                           ax.plot(data[:,0],data[:,1], "k-")
#                           xplt=np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
#                           ax.plot(xplt,(bckg+gauss(xplt, P0m[0], -P0m[3], P0m[6])), color[i%3])
#                 plt.savefig(controlfits+foldername[k] +'_line_' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z))+'_fit.png')
#                 plt.close(fig)

"""
This block copies the fits in a common folder just 
for the sake of simplicity
"""
if os.path.exists(allcontrolfits):
    shutil.rmtree(allcontrolfits)
os.makedirs(allcontrolfits)
for k in range(len(foldername)):
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    folder = sorted_fold_path+foldername[k]
    roi =  np.loadtxt(data_analysis+foldername[k]+'_ROI+Peaks.mpa',skiprows=1).astype(int)
    for y in range(len(roi[:,0])):
        for z in range(1,n_theta[k]+1):
            contfitname = foldername[k] +'_line_' +str("%0d"%(roi[0][0]+y))+'_theta'+str("%0d"%(z))+'_fit.png'
            try:
                shutil.copy(folder+"/Data Analysis/Control Fits/"+contfitname, allcontrolfits+contfitname)        
            except FileNotFoundError:
                a=0
#                 print("not there")