# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:41:38 2022

@author: ismae
"""

import os
import shutil
import numpy as np
from PIL import Image as im

sorted_fold_path="C:/Users/ismae/OneDrive/Desktop/Physics/Thesis/Script/Trial/Sorted data/" #insert folder of sorted meausements files
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
"""
This block renames the data files so they respect a specific labeling pattern: 
"Tilt angle_number of measurement along theta", for example "0deg_001"
It has to be done folder by folder
"""
# orig_fold_path="C:/Users/ismae/OneDrive/Desktop/Physics/Thesis/Script/meaningfulDataFolder/0deg/" 
# renamed = sorted_fold_path+"0deg/Renamed/"
# # if os.path.exists(renamed):
# #     shutil.rmtree(renamed)
# # os.makedirs(renamed)
# shutil.copy(orig_fold_path+"#1AR_@0_offBrag@min10deg"+str("%03d" % (i,))+".mpa", renamed+"0deg_offBrag.mpa")

"""
This block creates for each measurement a file with the matrix version of the 
128x128 measurements
"""
data=np.zeros(n_pixel);
data_aus=data
matrixes = sorted_fold_path+foldername[0]+"Matrixes/"
f = open(sorted_fold_path+foldername[0]+"Renamed/leer.mpa",'r')
for i,line in enumerate(f):
    if i>124 and i<(126+n_pixel):
        data[i-125]=float(line)
        data_aus=data
f.close()
data_matrix = data_aus.reshape(128, 128)
with open(matrixes+str(tiltangles[0])+"leer.mpa", 'w') as f:
        np.savetxt(f, data_matrix, fmt="%1.0f")
"""
This creates for each matrix file a corresponding picture (flipped of 180 because im.fromarray gives
a rotated pic for some reason)
"""
max_count = 0.
for k in range(len(tiltangles)):
    for j in range (1,n_theta[k]+1):
        f = sorted_fold_path+foldername[k]+"Matrixes/"+str(tiltangles[k])+"deg_"+str("%03d" % (j,))+".mpa"
        data_matrix = np.loadtxt(f)
        if np.amax(data_matrix)>max_count:
            max_count=np.amax(data_matrix)
k=0
pictures = sorted_fold_path+foldername[k]+"Pictures/"
rawpictures = sorted_fold_path+foldername[k]+"Raw pictures/"
f = sorted_fold_path+foldername[k]+"Matrixes/"+str(tiltangles[k])+"leer.mpa"
data_matrix = np.loadtxt(f).astype(float)/max_count
data_raw = np.loadtxt(f).astype(float)
image = im.fromarray((data_matrix))
image=image.rotate(180)
image.save(pictures+str(tiltangles[k])+"leer.tiff")
imageraw = im.fromarray(data_raw)
imageraw = imageraw.rotate(180)
imageraw.save(rawpictures+str(tiltangles[k])+"leer.tiff")

"""
This creates common folders "Renamed", "Matrixes and "Pictures" for all measurements at all tilt angles 
"""

# if os.path.exists(allmeasurements):
#         shutil.rmtree(allmeasurements)
# os.makedirs(allmeasurements)
# if os.path.exists(allpictures):
#         shutil.rmtree(allpictures)
# os.makedirs(allpictures)
# if os.path.exists(allrawpictures):
#         shutil.rmtree(allrawpictures)
# os.makedirs(allrawpictures)
# if os.path.exists(allmatrixes):
#     shutil.rmtree(allmatrixes)
# os.makedirs(allmatrixes)
# if os.path.exists(allrenamed):
#     shutil.rmtree(allrenamed)
# os.makedirs(allrenamed)

# for k in range(len(tiltangles)):
#     for j in range (1,n_theta[k]+1):
#         filename = str(tiltangles[k])+"deg_"+str("%03d" % (j,))+".mpa"
#         picname = str(tiltangles[k])+"deg_"+str("%03d" % (j,))+".tiff"
#         folder = sorted_fold_path+foldername[k]
#         shutil.copy(folder+"Renamed/"+filename, allrenamed+filename)
#         shutil.copy(folder+"Matrixes/"+filename, allmatrixes+filename)
#         shutil.copy(folder+"Pictures/"+picname, allpictures+picname)
#         shutil.copy(folder+"Raw pictures/"+picname, allrawpictures+picname)
