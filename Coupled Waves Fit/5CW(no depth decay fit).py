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
import scipy.integrate as integrate
import math
from scipy.interpolate import UnivariateSpline
from datetime import datetime
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
lambda_par=1595.7292122046995	#+/-	147.471394720765
mu=3.2e-3#0.004632543663155012	#+/-	5.46776175965519e-05
sigma=0.0006873016655595522	
##############################################################################

"""
Angular distribution: Gaussian
"""
div=0.0006
def ang_gauss(x,x0):
    sig=div
    return 1/((2*pi)**0.5*sig)*np.exp(-(x-x0)**2/(2*sig**2))


##############################################################################

n_diff= 3 #number of peaks for each side, for example: n=2 for 5 diffracted waves

LAM= 0.5 #grating constant in micrometers
G=2*pi/LAM
bcr1=5.0 #scattering lenght x density
bcr2=2
bcr3=0
n_0 =1.
phi=pi
d=78
#print(n_1)

def k_jz(theta, j, G,b):
    k_jz=b*(1-(np.sin(theta)-j*G/b)**2)**0.5
    return k_jz
def dq_j (theta, j, G,b):
    return b*np.cos(theta) - k_jz(theta, j, G, b)
fitting=0
plotting=1
for k in range(0,1):#len(foldername)):
    wlpoints=10
    print(foldername[k])
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+'_diff_eff.mpa',skiprows=1)
    fit_res =  np.loadtxt(data_analysis+foldername[k]+'_fit_results.mpa',skiprows=1)
    for i in range(len(diff_eff[:,0])): 
        s=sum(diff_eff[i,1:])
        diff_eff[i,1:]=diff_eff[i,1:]/s
    diff_eff_fit=np.zeros((5, len(diff_eff[:,5])))
    diff_eff_fit[2,:]=diff_eff[:,2*2+2].copy()
    for i in range(1,3):
        diff_eff_fit[2-i,:]=diff_eff[:,6-2*i].copy()
        diff_eff_fit[2+i,:]=diff_eff[:,6+2*i].copy()
    def fit_func(x, bcr1, bcr2, mu1, phi,d):
        d=d/np.cos(tilt[k]*rad)
        wl=np.linspace(mu1-3.5*sigma, mu1+7*sigma, 10000)
        a = rho(wl,lambda_par, mu1, sigma)/sum(rho(wl,lambda_par, mu1, sigma))
        spl = UnivariateSpline(wl, a, k=4, s=0)
        I=spl.antiderivative()(wl)
        y=np.linspace(I[I==np.amin(I)],I[I==np.amax(I)],  wlpoints)
        xp=np.zeros(wlpoints)
        for i in range(wlpoints):
            aus =abs(spl.antiderivative()(wl)-y[i])
            xp[i]=wl[aus==np.amin(aus)]
        wl=xp.copy()
        a=rho(xp,lambda_par, mu1, sigma)/sum(rho(xp,lambda_par, mu1, sigma))
        th=np.linspace(diff_eff[0,0]*rad-3*div,diff_eff[-1,0]*rad+3*div, len(diff_eff[:,0])*wlpoints//10+10)
        tx=np.zeros(len(diff_eff[:,0]),dtype=int)
        for i in range(len(diff_eff[:,0])):
            for j in range(1,len(th)-1):
                if(th[j-1]<=x[i]*rad and th[j+1]>=x[i]*rad):
                    tx[i]=j
        for i in range(len(tx)):
            th[tx[i]]=x[i]*rad
        
        # plt.plot(x[0:len(tx)]*rad,x[0:len(tx)]*0,"k.")
        # plt.plot(th[tx],th[tx]*0 +1,"b.")
        # print("here")
        S=np.zeros((2*n_diff+1,len(th)),dtype=complex)
        eta=S.copy().real
        eta_aus=eta.copy()
        sum_diff = np.zeros(len(th))
        for l in range(len(wl)):
            lam=wl[l] #single wavelenght in micrometers
            b=2*pi/lam #beta value 
            n_1 = bcr1*2*pi/b**2
            n_2 = bcr2*2*pi/b**2
            # n_3 = bcr3*2*pi/b**2
            for t in range(len(th)):
                A = np.zeros((2*n_diff+1,2*n_diff+1), dtype=complex)
                for i in range(len(A[0])):
                    A[i][i]=b**2*(n_0**2-1)/(2*k_jz(th[t],i-n_diff,G,b))-dq_j(th[t],i-n_diff,G,b)
                    if(i+1<len(A[0])):
                        A[i][i+1]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G,b))
                        A[i+1][i]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G,b))
                    if(i+2<len(A[0]) and bcr2!=0):
                        A[i][i+2]=b**2*n_0*n_2*np.exp(1j*phi)/(2*k_jz(th[t],i-n_diff,G,b))
                        A[i+2][i]=b**2*n_0*n_2*np.exp(-1j*phi)/(2*k_jz(th[t],i-n_diff,G,b))
                    # if(i+3<len(A[0]) and bcr3!=0):
                    #     A[i][i+3]=b**2*n_0*n_3/(2*k_jz(th[t],i-n_diff,G,b))
                    #     A[i+3][i]=b**2*n_0*n_3/(2*k_jz(th[t],i-n_diff,G,b))
                A=-1j*A
                w,v = np.linalg.eig(A)
                v0=np.zeros(2*n_diff+1)
                v0[n_diff]=1
                c = np.linalg.solve(v,v0)
                for i in range(len(w)):
                    v[:,i]=v[:,i]*c[i]*np.exp(w[i]*d)
                for i in range(len(S[:,0])):
                    S[i,t] = sum(v[i,:])
            for t in range(len(th)):
                for i in range(2*n_diff+1):
                    eta_aus[i,t] = abs(S[i,t])**2*k_jz(th[t],i-n_diff,G,b)/(b*np.cos(th[t]))
                sum_diff[t] = sum(eta[:,t])
            eta+=eta_aus*a[l]
        eta_ang = eta.copy()
        for i in range(len(eta[:,0])):
            for j in range(len(tx)):
                eta_ang[i,tx[j]] = sum(ang_gauss(th,th[tx[j]])*eta[i,:])/sum(ang_gauss(th,th[tx[j]]))
        eta_fit=np.zeros((2*n_diff+1,len(diff_eff[:,0])))
        for i in range(len(tx)):
            eta_fit[:,i]=eta_ang[:,tx[i]]
        aaa=eta_fit[n_diff-2:n_diff+3].ravel()
        #plt.plot(aaa)
        return aaa
    now=datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    P0= fit_res[0]
    if (fitting):
        B=([3, 0, 2e-3,-2*pi,65],[8, 5, 5e-3,2*pi,90])
        ff=diff_eff_fit.ravel()
        xx=np.zeros(len(diff_eff[:,0])*5)
        xx[0:len(diff_eff[:,0])]=diff_eff[:,0]
       # plt.plot(ff,"k")
        try:
            for i in range(3):
                p,cov=fit(fit_func,xx,ff, p0=P0, bounds=B)
                P0=p
                wlpoints=10**i*wlpoints
                print(p)
        except RuntimeError:
            print("Error: fit not found")
        wlpoints=wlpoints//10
        print(p)
        print(np.diag(cov)**0.5)
        now1=datetime.now()
        print("fit time=",now1-now)
        # with open(data_analysis+foldername[k]+'_fit_results.mpa', 'w') as f:
        #     np.savetxt(f,(p,np.diag(cov)**0.5), header="bcr1 bcr2 mu phi thickness", fmt="%.6f")
    print(wlpoints)


    def plot_func(th, bcr1, bcr2, mu1,phi,d):
        d=d/np.cos(tilt[k]*rad)
        wl=np.linspace(mu1-3.5*sigma, mu1+7*sigma, 10000)
        a = rho(wl,lambda_par, mu1, sigma)/sum(rho(wl,lambda_par, mu1, sigma))
        # plt.plot(wl,a/np.amax(a))
        spl = UnivariateSpline(wl, a, k=4, s=0)
        I=spl.antiderivative()(wl)
        y=np.linspace(I[I==np.amin(I)],I[I==np.amax(I)],  wlpoints)
        xp=np.zeros(wlpoints)
        for i in range(wlpoints):
            aus =abs(spl.antiderivative()(wl)-y[i])
            xp[i]=wl[aus==np.amin(aus)]
        wl=xp.copy()
        a=rho(xp,lambda_par, mu1, sigma)/sum(rho(xp,lambda_par, mu1, sigma))
        # plt.plot(xp,xp*0,"k.")
        # print(sum(a))
        # plt.plot(xp,a/np.amax(a),".")
        S=np.zeros((2*n_diff+1,len(th)),dtype=complex)
        eta=S.copy().real
        eta_aus=eta.copy()
        sum_diff = np.zeros(len(th))
        for l in range(len(wl)):
            lam=wl[l] #single wavelenght in micrometers
            b=2*pi/lam #beta value 
            n_1 = bcr1*2*pi/b**2
            n_2 = bcr2*2*pi/b**2
            n_3 = bcr3*2*pi/b**2
            for t in range(len(th)):
                A = np.zeros((2*n_diff+1,2*n_diff+1), dtype=complex)
                for i in range(len(A[0])):
                    A[i][i]=b**2*(n_0**2-1)/(2*k_jz(th[t],i-n_diff,G,b))-dq_j(th[t],i-n_diff,G,b)
                    if(i+1<len(A[0])):
                        A[i][i+1]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G,b))
                        A[i+1][i]=b**2*n_0*n_1/(2*k_jz(th[t],i-n_diff,G,b))
                    if(i+2<len(A[0]) and bcr2!=0):
                        A[i][i+2]=b**2*n_0*n_2*np.exp(1j*phi)/(2*k_jz(th[t],i-n_diff,G,b))
                        A[i+2][i]=b**2*n_0*n_2*np.exp(-1j*phi)/(2*k_jz(th[t],i-n_diff,G,b))
                    if(i+3<len(A[0]) and bcr3!=0):
                        A[i][i+3]=-b**2*n_0*n_3/(2*k_jz(th[t],i-n_diff,G,b))
                        A[i+3][i]=-b**2*n_0*n_3/(2*k_jz(th[t],i-n_diff,G,b))
                A=-1j*A
                w,v = np.linalg.eig(A)
                v0=np.zeros(2*n_diff+1)
                v0[n_diff]=1
                c = np.linalg.solve(v,v0)
                for i in range(len(w)):
                    v[:,i]=v[:,i]*c[i]*np.exp(w[i]*d)
                for i in range(len(S[:,0])):
                    S[i,t] = sum(v[i,:])
            for t in range(len(th)):
                for i in range(2*n_diff+1):
                    eta_aus[i,t] = abs(S[i,t])**2*k_jz(th[t],i-n_diff,G,b)/(b*np.cos(th[t]))
                sum_diff[t] = sum(eta[:,t])
            eta+=eta_aus*a[l]
        eta_ang = eta.copy()
        for i in range(len(eta[:,0])):
            for j in range(len(eta[0,:])):
                eta_ang[i,j] = sum(ang_gauss(th,th[j])*eta[i,:])/sum(ang_gauss(th,th[j]))
        return eta_ang
    
    if (plotting):
        if (not fitting):
            p=P0
            print(p)
            now1=datetime.now()
            wlpoints=500
        thx=np.linspace(diff_eff[0,0]*rad,diff_eff[-1,0]*rad, 100)
        eta=plot_func(thx, *p)
        #print("here")
        # bbb=eta.ravel()
        # plt.plot(bbb)
        fig, ax = plt.subplots(n_diff+2,figsize=(10,10))
        ax[0].set_title(foldername[k])
        ax[0].plot(thx,eta[n_diff,:])
        ax[0].plot(diff_eff[:,0]*rad,diff_eff_fit[2,:],'o')
        for i in range(1,n_diff+1):
            ax[i].plot(thx,eta[n_diff-i,:])
            ax[i].plot(thx,eta[n_diff+i,:])   
            if i<3:
                ax[i].plot(diff_eff[:,0]*rad,diff_eff[:,6-2*i],'o')
                ax[i].plot(diff_eff[:,0]*rad,diff_eff[:,6+2*i],'o')
        # ax[n_diff+1].plot(th, sum_diff)
        # ax[n_diff+1].set_ylim([0.5,1.5])
        #   plt.errorbar(diff_eff[:,0],diff_eff[:,2*j+2],yerr=diff_eff[:,2*j+1],capsize=1)
        now2=datetime.now()
        print("plot time=",now2-now1)
duration = 0.2  # seconds

freq = 440  # Hz
for i in range (6):
    os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq+i%3*62))
    if i%3==2:
        os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))
for i in range (2):
    os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq+2*62))
    os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq+2*62+31))
    os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq+3*62+31))
    time.sleep(0.2)