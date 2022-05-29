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
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib as mpl
import socket
import shutil
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
import scipy.integrate as integrate
import math
from datetime import datetime
from scipy.interpolate import interp1d
pi=np.pi
rad=pi/180

sorted_fold_path="/home/exp-k03/Desktop/thesis/Sorted data/" #insert folder of sorted meausements files
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
##############################################################################
"""
Wavelenght distribution: Exponentially Modified Gaussian
"""
def func(l,A,mu,sig):
    return A/(2.)*np.exp(A/(2.)*(2.*mu+A*sig**2-2*l))
def rho(l,A,mu,sig):
    return func(l,A,mu,sig)*erfc((mu+A*sig**2-l)/(np.sqrt(2)*sig))

# const=0.0052592164045898665
# 	#+/-	147.471394720765
mu=3e-3 #0.004632543663155012	#+/-	5.46776175965519e-05
#lambda_par= 1/(const-mu)
lambda_par=500
sigma=0.005
sigma1= (sigma**2+1/lambda_par**2)**0.5

##############################################################################

"""
Angular distribution: Gaussian
"""
div=0.0006/2
def ang_gauss(x,x0):
    sig=div
    return 1/((2*pi)**0.5*sig)*np.exp(-(x-x0)**2/(2*sig**2))


##############################################################################
k=2
n_diff= 4 #number of peaks for each side, for example: n=2 for 5 diffracted waves

LAM= 0.5 #grating constant in micrometers
G=2*pi/LAM
bcr1=5.0 #scattering lenght x density
bcr2=-2
bcr3=0
n_0 =1.
phi=-pi
wl=np.linspace(mu-2.5*sigma,mu+1/lambda_par+3*sigma1, 10000)
a = rho(wl,lambda_par, mu, sigma)/sum(rho(wl,lambda_par, mu, sigma))
from scipy.interpolate import UnivariateSpline
spl = UnivariateSpline(wl, a, k=3, s=0)
d=spl.antiderivative()(wl)
s=50
y=np.arange(d[d==np.amin(d)],d[d==np.amax(d)]+1e-8,  2e-8)
x=np.zeros(len(y))
guess=0
print(len(y))
for i in range(len(y)):
    aus =abs(spl.antiderivative()(wl)-y[i])
    x[i]=wl[aus==np.amin(aus)]
plt.plot(wl,d/np.amax(d))
plt.plot(wl,a/np.amax(a),"b")
plt.plot(x,x*0,"k.")
a=rho(x,lambda_par, mu, sigma)/sum(rho(x,lambda_par, mu, sigma))
plt.plot(x,a/np.amax(a),"g.")
print(a[0],a[-1])
# a = rho1(wl,lambda_par, mu, sigma)/sum(rho1(wl,lambda_par, mu, sigma))
# from scipy.interpolate import UnivariateSpline
# spl = UnivariateSpline(wl, a, k=3, s=0)
# d=spl.antiderivative()(wl)
# s=100
# y=np.linspace(d[d==np.amin(d)],d[d==np.amax(d)],  s)
# x=np.zeros(s)
# for i in range(s):
#     aus =abs(spl.antiderivative()(wl)-y[i])
#     x[i]=wl[aus==np.amin(aus)]
# plt.plot(wl,d/np.amax(d))
# plt.plot(wl,a/np.amax(a))
# plt.plot(x,x*0,"k.")
# a=rho1(x,lambda_par, mu, sigma)/sum(rho1(x,lambda_par, mu, sigma))
# plt.plot(x,a/np.amax(a),"g.")

# print(a[0],a[-1])
# data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
# diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff.mpa',skiprows=1)
# x1=diff_eff[:,0]*rad
# th=np.linspace(x1[0]-3*div,x1[-1]+3*div,10000)
# print(x1)
# asd=np.zeros(10000)
# for i in range(len(x1)):
#     asd += ang_gauss(th,x1[i])
# spl=UnivariateSpline(th, asd, k=3, s=0)
# d=spl.antiderivative()(th)
# plt.plot(th, asd/np.amax(asd))
# plt.plot(th, d/np.amax(d))
# s=len(x1)*100
# y=np.linspace(d[d==np.amin(d)],d[d==np.amax(d)],  s)
# x=np.zeros(s)
# for i in range(s):
#     aus =abs(spl.antiderivative()(th)-y[i])
#     x[i]=th[aus==np.amin(aus)]
# plt.plot(x,x*0,"k.")
# plt.plot(x1,x1*0,"r.")


# plt.plot(qwe) 
# print(sum(qwe)/len(qwe))