#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import CubicSpline
from scipy.interpolate import PPoly
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('font',**{'size':18})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


Tpc_latt_2L = 167
Tpc_latt = 182

plotpath = '/Users/antoniosmecca/Documents/Physics/pdoc_Swansea/mu2/Code'

mu=float(sys.argv[1])

a_inv_gev_2L = 6.079
a_inv_gev = 5.63

################### BIS
if mu == 0.0 :
    filename= "gen2/interpolate/boot_zero.dat"
    filename_2L= "Gen2L/interpolate/boot_gen2l_zero.dat"
else:
    filename= f"gen2/interpolate/boot_full_mu_0.056.dat"
    filename_2L= f"Gen2L/interpolate/boot_gen2l_full_mu_{mu}.dat"

f=open(filename,'r')

T,R,R_err=[],[],[]
T_2L,R_2L,R_err_2L=[],[],[]
Nt_array=[]
Nt_array_2L=[]
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

f=open(filename_2L,'r')
for line in f.readlines():
    if line[0] != '#':
        x=line.split()
        T_2L += [int(x[0])]
        R_2L += [float(x[1])]
        len_Nt_2L = int(x[2])
        #R_err += [float(x[2])]
T_2L = np.asarray(T_2L)
R_2L = np.asarray(R_2L)
#R_err = np.asarray(R_err)
f.close()

R_boot = np.zeros((len_Nt,boot_samples),dtype=float)
T_boot = np.zeros((len_Nt,boot_samples),dtype=float)
for i in range(len_Nt):
    for b in range(0,boot_samples):
        R_boot[i][b] = R[i*boot_samples+b]
        T_boot[i][b] = T[i*boot_samples+b]
R_boot_2L = np.zeros((len_Nt_2L,boot_samples),dtype=float)
T_boot_2L = np.zeros((len_Nt_2L,boot_samples),dtype=float)
for i in range(len_Nt_2L):
    for b in range(0,boot_samples):
        R_boot_2L[i][b] = R_2L[i*boot_samples+b]
        T_boot_2L[i][b] = T_2L[i*boot_samples+b]

for i in range(len_Nt):
    for b in range(0,boot_samples):
        T_boot[i][b] = (1/T_boot[i][b])*a_inv_gev
T_boot = np.flip(T_boot)
R_boot = np.flip(R_boot)
for i in range(len_Nt_2L):
    for b in range(0,boot_samples):
        T_boot_2L[i][b] = (1/T_boot_2L[i][b])*a_inv_gev_2L
T_boot_2L_plot = np.flip(T_boot_2L)
R_boot_2L_plot = np.flip(R_boot_2L)

T_boot_2L = T_boot_2L_plot
R_boot_2L = R_boot_2L_plot

# Removing T = 169 from Gen2L data
if mu != 0.0:
    print('length: ',len(T_boot_2L_plot))
    T_boot_2L = np.delete(T_boot_2L_plot,1,axis=0)
    R_boot_2L = np.delete(R_boot_2L_plot,1,axis=0)
    len_Nt_2L = len_Nt-1
    print('length: ',len(T_boot_2L))


T_mean = np.zeros(len_Nt,dtype=float)
R_mean = np.zeros(len_Nt,dtype=float)
R_stdev = np.zeros(len_Nt,dtype=float)
for i in range(len_Nt):
    T_mean[i] = np.mean(T_boot[i])
    R_mean[i] = np.mean(R_boot[i])
    R_stdev[i] = np.std(R_boot[i])
T_mean_2L = np.zeros(len_Nt_2L,dtype=float)
R_mean_2L = np.zeros(len_Nt_2L,dtype=float)
R_stdev_2L = np.zeros(len_Nt_2L,dtype=float)
T_mean_2L_plot = np.zeros(len_Nt_2L+1,dtype=float)
R_mean_2L_plot = np.zeros(len_Nt_2L+1,dtype=float)
R_stdev_2L_plot = np.zeros(len_Nt_2L+1,dtype=float)
for i in range(len_Nt_2L):
    T_mean_2L[i] = np.mean(T_boot_2L[i])
    R_mean_2L[i] = np.mean(R_boot_2L[i])
    R_stdev_2L[i] = np.std(R_boot_2L[i])
for i in range(len_Nt_2L): #+1
    T_mean_2L_plot[i] = np.mean(T_boot_2L_plot[i])
    R_mean_2L_plot[i] = np.mean(R_boot_2L_plot[i])
    R_stdev_2L_plot[i] = np.std(R_boot_2L_plot[i])

x=np.arange((1/41)*a_inv_gev,(1/25)*a_inv_gev,0.0001)
x_plot=np.arange((1/41)*a_inv_gev,(1/15)*a_inv_gev,0.0001)
cs = CubicSpline(T_mean, R_mean)
cs_boot = np.zeros(boot_samples,dtype=PPoly)
x_2L=np.arange((1/41)*a_inv_gev_2L,(1/25)*a_inv_gev_2L,0.0001)
x_plot_2L=np.arange((1/41)*a_inv_gev_2L,(1/15)*a_inv_gev_2L,0.0001)
cs_2L = CubicSpline(T_mean_2L, R_mean_2L)
cs_boot_2L = np.zeros(boot_samples,dtype=PPoly)

R_t = np.transpose(R_boot)
for i in range(0,boot_samples):
    cs_boot[i] = CubicSpline(T_mean,R_t[i])
R_t_2L = np.transpose(R_boot_2L)
for i in range(0,boot_samples):
    cs_boot_2L[i] = CubicSpline(T_mean_2L,R_t_2L[i])

cs_plus_2L = CubicSpline(T_mean_2L, (R_mean_2L+R_stdev_2L))
cs_minus_2L = CubicSpline(T_mean_2L, (R_mean_2L-R_stdev_2L))
cs_plus = CubicSpline(T_mean, (R_mean+R_stdev))
cs_minus = CubicSpline(T_mean, (R_mean-R_stdev))

i=len(x)-1
while cs(x[i]) > 0 and i >0:
    Tpc = x[i]
    i -= 1
i=len(x)-1
while cs_2L(x_2L[i]) > 0 and i >0:
    Tpc_2L = x_2L[i]
    i -= 1


Tpc_boot = np.zeros(boot_samples,dtype=float)
for b in range(0,boot_samples):
    i=len(x)-1
    while cs_boot[b](x[i]) > 0 and i > 0:
        Tpc_boot[b] = x[i]
        i -= 1
Tpc_boot_2L = np.zeros(boot_samples,dtype=float)
for b in range(0,boot_samples):
    i=len(x)-1
    while cs_boot_2L[b](x_2L[i]) > 0 and i > 0:
        Tpc_boot_2L[b] = x_2L[i]
        i -= 1

Tpc_err = np.std(Tpc_boot)
Tpc_err_2L = np.std(Tpc_boot_2L)


mu_mev = int(float(mu)*1000)

plt.figure(figsize=(8.5, 5.7))
#plt.title(r'$\mu_q = $ '+str(mu_mev)+' $\mathrm{MeV}$')
plt.text(185,-0.15,r'$\mu_q = $ '+str(mu_mev)+' $\mathrm{MeV}$',fontsize=20)
plt.xlabel(r'$T$ $[\mathrm{MeV}]$',fontsize=18)
plt.ylabel(r'$\overline{R}(\mu_q,T)$',fontsize=18)
if mu != 0.084:
    plt.ylim(-0.25,0.15)
plt.xlim(140,203)
for i in range(len_Nt):    
    if i == 0:
        if mu == 0.0:
            plt.errorbar(x=T_mean[i]*1000,y=R_mean[i],yerr=R_stdev[i],marker='o',color='#1f77b4',label='Generation 2')
        else:
            plt.errorbar(x=T_mean[i]*1000,y=R_mean[i],yerr=R_stdev[i],marker='o',color='#1f77b4',label=r'Generation 2') 
    else:
        if mu == 0.0:
            plt.errorbar(x=T_mean[i]*1000,y=R_mean[i],yerr=R_stdev[i],marker='o',color='#1f77b4')
        else:
            plt.errorbar(x=T_mean[i]*1000,y=R_mean[i],yerr=R_stdev[i],marker='o',color='#1f77b4')
for i in range(len_Nt_2L+1):    
    if i == 0.0:
        if mu == 0.0:
            plt.errorbar(x=T_mean_2L_plot[i]*1000,y=R_mean_2L_plot[i],yerr=R_stdev_2L_plot[i],marker='s',color='#ff7f0e',label=r'Generation 2L')
        else:
            plt.errorbar(x=T_mean_2L_plot[i]*1000,y=R_mean_2L_plot[i],yerr=R_stdev_2L_plot[i],marker='s',color='#ff7f0e',label=r'Generation 2L') 
    else:
        if mu == 0.0:
            plt.errorbar(x=T_mean_2L_plot[i]*1000,y=R_mean_2L_plot[i],yerr=R_stdev_2L_plot[i],marker='s',color='#ff7f0e')
        else:
            plt.errorbar(x=T_mean_2L_plot[i]*1000,y=R_mean_2L_plot[i],yerr=R_stdev_2L_plot[i],marker='s',color='#ff7f0e')
plt.plot(x_plot*1000,cs(x_plot),color='#1f77b4')
plt.fill_between(x_plot*1000,cs_plus(x_plot),cs_minus(x_plot),alpha=0.2,color='#1f77b4')
plt.plot(x_plot_2L*1000,cs_2L(x_plot_2L),color='#ff7f0e')
plt.fill_between(x_plot_2L*1000,cs_plus_2L(x_plot_2L),cs_minus_2L(x_plot_2L),alpha=0.2,color='#ff7f0e')
plt.axhline(y=0.0,linestyle='--',color='darkgray')
plt.axvline(x=Tpc_latt,linestyle='--',color='black')
plt.errorbar(x=182,xerr=2,y=0,color='black',capsize=4)

plt.errorbar(x=167,xerr=3,y=0,color='black',capsize=4)

plt.axvline(x=Tpc_latt_2L,linestyle='dotted',color='black')
plt.legend(loc='best')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/2L2_interpolate_full_mu_mu_{mu}.pdf",dpi=300)
plt.show()

