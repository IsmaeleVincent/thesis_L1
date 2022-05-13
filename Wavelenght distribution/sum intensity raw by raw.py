# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:42:45 2022

@author: ismae
"""
"""
ATTENTION: All the blocks of this script completely re-write folders if they already exist
"""
import os
import shutil
import numpy as np
from PIL import Image as im
import matplotlib.pyplot as plt

sorted_fold_path="/home/aaa/Desktop/Thesis/Script/Trial/Sorted data/" #insert folder of sorted meausements files
allmeasurements = sorted_fold_path+"All measurements/"
allrenamed = allmeasurements +"All renamed/"
allmatrixes = allmeasurements + "All matrixes/"
allpictures = allmeasurements + "All pictures/"
allrawpictures = allmeasurements + "All raw pictures/"
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg/")
n_theta=[25,45,27,17,15,19,20,19,18] #number of measurements along theta
n_pixel = 16384 #number of pixels in one measurement
x_ROI = [38,89]
y_ROI = [66,76]
vert_profile = np.zeros((y_ROI[1]-y_ROI[0]+1,2))
matrix_aus= np.zeros((x_ROI[1]-x_ROI[0]+1,y_ROI[1]-y_ROI[0]+1))
for k in range(0,1):#(len(tiltangles)):
    matrixes = sorted_fold_path+foldername[k]+"Matrixes/"
    for j in range (1,2):#(n_theta[k]+1):
        with open(matrixes+str(tiltangles[k])+"deg_"+str("%03d" % (j,))+".mpa",'r') as f:
            data_matrix = np.loadtxt(f).astype(float)
            #data_matrix.transpose()
            for y in range(y_ROI[0],y_ROI[1]+1):
                vert_profile[y-y_ROI[0]][0] = y
                for x in range(x_ROI[0],x_ROI[1]+1):
                    vert_profile[y-y_ROI[0]][1] += data_matrix[128-y][x]
                    matrix_aus[x-x_ROI[0]][y-y_ROI[0]]= data_matrix[128-y][x]
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
ax.plot(vert_profile[:,0],vert_profile[:,1],"b-")
#image = im.fromarray((matrix_aus))
#im=ax.imshow(matrix_aus,cmap='plasma')
# print(vert_profile)
# with open("C:/Users/ismae/OneDrive/Desktop/Physics/Thesis/Script/Trial/sum.txt", 'w') as f:
#     np.savetxt(f, vert_profile, fmt="%1.0f")

