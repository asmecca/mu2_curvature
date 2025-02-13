#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import CubicSpline
from scipy.interpolate import PPoly

plt.rc('font',**{'size':11})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


Tpc_latt = 182

plotpath = '/Users/antoniosmecca/Documents/Physics/pdoc_Swansea/mu2/Code'

mu=float(sys.argv[1])
a_inv_gev = 5.63

################### BIS
if mu == 0.0 :
    filename= "interpolate/boot_zero.dat"
    f=open(filename,'r')
else:
    filename= f"interpolate/boot_full_mu_{mu}.dat"
    f=open(filename,'r')

T,R,R_err=[],[],[]
Nt_array=[]
for line in f.readlines():
    if line[0] != '#':
        x=line.split()
        T += [int(x[0])]
        R += [float(x[1])]
        len_Nt = int(x[2])
        boot_samples = int(x[3])
T = np.asarray(T)
R = np.asarray(R)
#R_err = np.asarray(R_err)
f.close()

R_boot = np.zeros((len_Nt,boot_samples),dtype=float)
T_boot = np.zeros((len_Nt,boot_samples),dtype=float)
for i in range(len_Nt):
    for b in range(0,boot_samples):
        R_boot[i][b] = R[i*boot_samples+b]
        T_boot[i][b] = T[i*boot_samples+b]

for i in range(len_Nt):
    for b in range(0,boot_samples):
        T_boot[i][b] = (1/T_boot[i][b])*a_inv_gev
T_boot_plot = np.flip(T_boot)
R_boot_plot = np.flip(R_boot)

# Removing T=169n MeV point from interpolation
print('length: ',len(T_boot))
T_boot = np.delete(T_boot_plot,1,axis=0)
R_boot = np.delete(R_boot_plot,1,axis=0)
len_Nt = len_Nt-1
print('length: ',len(T_boot))

T_mean = np.zeros(len_Nt,dtype=float)
T_mean_plot = np.zeros(len_Nt,dtype=float)
R_mean = np.zeros(len_Nt,dtype=float)
R_mean_plot = np.zeros(len_Nt,dtype=float)
R_stdev = np.zeros(len_Nt,dtype=float)
R_stdev_plot = np.zeros(len_Nt,dtype=float)
for i in range(len_Nt):
    T_mean[i] = np.mean(T_boot[i])
    T_mean_plot[i] = np.mean(T_boot_plot[i])
    R_mean[i] = np.mean(R_boot[i])
    R_stdev[i] = np.std(R_boot[i])
    R_mean_plot[i] = np.mean(R_boot_plot[i])
    R_stdev_plot[i] = np.std(R_boot_plot[i])

x=np.arange((1/41)*a_inv_gev,(1/25)*a_inv_gev,0.0001)
x_plot=np.arange((1/41)*a_inv_gev,(1/15)*a_inv_gev,0.0001)
cs = CubicSpline(T_mean, R_mean)
cs_boot = np.zeros(boot_samples,dtype=PPoly)

R_t = np.transpose(R_boot)
for i in range(0,boot_samples):
    cs_boot[i] = CubicSpline(T_mean,R_t[i])

cs_plus = CubicSpline(T_mean, (R_mean+R_stdev))
cs_minus = CubicSpline(T_mean, (R_mean-R_stdev))

i=len(x)-1
while cs(x[i]) > 0 and i >0:
    Tpc = x[i]
    i -= 1


Tpc_boot = np.zeros(boot_samples,dtype=float)
for b in range(0,boot_samples):
    i=len(x)-1
    while cs_boot[b](x[i]) > 0 and i > 0:
        Tpc_boot[b] = x[i]
        i -= 1

Tpc_err = np.std(Tpc_boot)

mu_mev = int(float(mu)*1000)
plt.xlabel(r'$T/T_{c}$')
plt.ylabel(r'$R(\mu_q;T)$')
plt.title(r'$\mu_q = $ '+str(mu_mev)+' $\mathrm{MeV}$')

for i in range(len_Nt):    
    plt.errorbar(x=T_mean_plot[i]*1000/Tpc_latt,y=R_mean_plot[i],yerr=R_stdev_plot[i],marker='o',color='blue')
plt.plot(x_plot*1000/Tpc_latt,cs(x_plot),color='darkorange')
plt.plot(x_plot*1000/Tpc_latt,cs_plus(x_plot),linestyle='--',color='darkorange')
plt.plot(x_plot*1000/Tpc_latt,cs_minus(x_plot),linestyle='--',color='darkorange')
plt.axhline(y=0.0,linestyle='--',color='darkgray')
plt.axvline(x=1,linestyle='--',color='darkgray')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/boot_interpolate_full_mu_method1_bis_mu_{mu}.png",dpi=300)
plt.show()

g=open("interpolate/boot_pseudo_T.dat",'a')
for b in range(0,boot_samples):
    g.write(str(mu)+' '+str(Tpc_boot[b])+' '+str(boot_samples)+'\n') 
g.close()
