 # -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 17:00:42 2022

@author: ismae
"""
from PIL import Image as im
import matplotlib.pyplot as plt
import numpy as np
import time
n_pixel = 16384
data=np.zeros(n_pixel);
data_aus=data
data_line = np.zeros(128)
folderspath="C:/Users/ismae/OneDrive/Desktop/Physics/Thesis/Script/meaningfulDataFolder/" #insert folder of meausements files
tiltangles=[0,40,48,61,69,71,79,80,81]
foldername=[]
n_line = 0 #number of line of interest
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg/Renamed/")
n_theta=[26,46,28,18,16,20,21,20,19] #number of measurements along theta
totalcounts = np.zeros(len(tiltangles))
tot=0
max_count = 0.
for k in range(1,2):#range(len(tiltangles)):
    for j in range (1,n_theta[k]):
        f = open(folderspath+foldername[k]+str(tiltangles[k])+"deg_"+str("%03d" % (j,))+".mpa",'r')
        for i,line in enumerate(f):
            # if i==0:
            #     for j,line in enumerate(f):
            #         if j>124 and j<(126+n_pixel):
            #             totalcounts[k]+=float(line)
            #             tot+=float(line)
            #             data[i-125]+=float(line)
            #             data_aus=data
            if i>124 and i<(126+n_pixel):
                totalcounts[k]+=float(line)
                data[i-125]=float(line)
                data_aus=data
                max_aus = np.amax(data)
                if max_aus > max_count:
                    max_count=max_aus
        # f.close()
        # data_aus_matrix = data_aus.reshape(128, 128)
        # fig = plt.figure(figsize=(15,15))
        # ax = fig.add_subplot(111)
        # ax.plot(data_line, "k--")
        # im=ax.imshow(data_aus_matrix,cmap='plasma')
        # ax.set_xlim([40,85])
        # ax.set_ylim([40,70])
        # plt.colorbar(im)
        # plt.show()
        # time.sleep(0.005)

for j in range (1,n_theta[0]):
    f = open(folderspath+foldername[0]+str(tiltangles[0])+"deg_"+str("%03d" % (j,))+".mpa",'r')
    data=np.zeros(n_pixel);
    for i,line in enumerate(f):
        if i>124 and i<(126+n_pixel):
            data[i-125]=float(line)/(max_count)
    data_matrix = data.reshape(128, 128)
    image = im.fromarray((data_matrix*255).astype('uint8'))
    image.save(str(j)+'.png')
# for n_line in range(55,65):
#     for i in range (128):
#         data_line[i] = data_matrix[n_line][i]
#     fig = plt.figure(figsize=(15,15))
#     ax = fig.add_subplot(111)
#     ax.plot(data_line, "k--")
#     time.sleep(5.0)
#     #ax.set_ylim([0,5000])
# # print(lines)
# # print(data)
# # print(data_matrix)
# fig = plt.figure(figsize=(15,15))
# ax = fig.add_subplot(111)
# im=ax.imshow(data_matrix,cmap='plasma')
# #ax.plot(tiltangles,totalcounts, "k--")
# #plt.legend(loc=1, prop={'size': 14}, frameon=True)
# ax.plot(data_line, "k--")
# #ax.yaxis.tick_left()
# #ax.set_xlabel('$m$',FontSize=20)
# #ax.set_ylabel('$k_B T$',FontSize=20)
# ax.tick_params(axis='both', which='major',length=8, labelsize=16)
# ax.tick_params(axis='both', which='minor',length=5)
#ax.set_xlim([40,80])
# ax.set_ylim([55,60])
# plt.colorbar(im,cmap='plasma')
# plt.show()

