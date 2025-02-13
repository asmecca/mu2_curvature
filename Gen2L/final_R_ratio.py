#!/usr/bin/env python3

# The iteration is: for each s_id there are all the timeslices, so for bin=1 we have from t=1 to t=Nt

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit
from utils import correlated_chi_squared
from utils import get_correlators

from cycler import cycler


def symmetrise(corr,Nt):
    for t in range(1,int(Nt/2)):
        corr[Nt-t] = corr[t]


plt.rc('font',**{'size':14})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

default_cycler = (cycler(marker=['v','s','^','*','o','<','p'])+cycler(color=['darkred','red','orange','gold','lightblue','#1f99b4','blue']))
plt.rc('axes', prop_cycle=default_cycler)
plt.rc('lines', linestyle='')

marker_list=['v','<','^','>','o','s','p']
color_list=['darkred','red','orange','gold','#71daeb','#0fb7f5','blue']

a_inv_gev = 6.079

T_array = [0.380, 0.304, 0.253, 0.217, 0.190, 0.169, 0.152] 
tminarray = [3, 4, 4, 5, 6, 7, 8]

Nt_array =[16, 20, 24, 28, 32, 36, 40]

boot_samples = 2000

mu = float(sys.argv[1])


path_to_corr = "DATA/Gen2L" # Path to directory with name NtxNs containing the correlators

#########################
f = open(path_to_corr+f"/16x32/analysis_mu_{mu}/boot/res.vector.dat",'r')
C=[]
for l in f.readlines():
    x=l.split()
    C += [float(x[0])]
boot_check = int(len(C)/16)
print('boot samples is: ',boot_check)
########################

corr_t,corr_g_zero,corr_g_zero_A,corr_total,corr_total_A = get_correlators(path_to_corr,boot_samples,Nt_array,T_array,a_inv_gev,mu)

plotpath = '..'
if not os.path.exists(plotpath+'/plots'):
        os.makedirs(plotpath+'/plots')


b_vec=[]
b_axial=[]

for i in range(0,len(Nt_array)):
    b_vec += [np.asarray(corr_g_zero[i])]
    b_axial += [np.asarray(corr_g_zero_A[i])]
    b_vec[i] = np.asarray(b_vec[i])
    b_axial[i] = np.asarray(b_axial[i])
    tmp1 = np.zeros((Nt_array[i],boot_samples),dtype=float)
    tmp2 = np.zeros((Nt_array[i],boot_samples),dtype=float)
    mean_a = np.zeros(Nt_array[i],dtype=float)
    stdev_a = np.zeros(Nt_array[i],dtype=float)
    mean_v = np.zeros(Nt_array[i],dtype=float)
    stdev_v = np.zeros(Nt_array[i],dtype=float)
    for t in range(0,Nt_array[i]):
        for b in range(0,boot_samples):
            tmp1[t][b] = float(b_axial[i][b*Nt_array[i]+t])
            tmp2[t][b] = float(b_vec[i][b*Nt_array[i]+t])
        mean_a[t] = np.mean(tmp1[t])
        stdev_a[t] = np.std(tmp1[t])
        mean_v[t] = np.mean(tmp2[t])
        stdev_v[t] = np.std(tmp2[t])


################ Order O(1) ratio ############

b_R=[]
b_R_num=[]
b_R_den=[]

for i in range(0,len(Nt_array)):
    R_num_tmp=[]
    R_den_tmp=[]
    for b in range(0,boot_samples):
        for t in range(0,Nt_array[i]):
            R_num_tmp += [b_vec[i][b*Nt_array[i]+t]/b_vec[i][b*Nt_array[i]+int(Nt_array[i]/2)] - b_axial[i][b*Nt_array[i]+t]/b_axial[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
            R_den_tmp += [b_axial[i][b*Nt_array[i]+t]/b_axial[i][b*Nt_array[i]+int(Nt_array[i]/2)] + b_vec[i][b*Nt_array[i]+t]/b_vec[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
    R_num_tmp = np.asarray(R_num_tmp)
    R_den_tmp = np.asarray(R_den_tmp)
    b_R_num += [R_num_tmp]
    b_R_den += [R_den_tmp]        


for i in range(0,len(Nt_array)):
    b_R += [b_R_num[i]/b_R_den[i]]

R_var_tot = [] #np.zeros((len(Nt_array),[]),dtype=float)
for i in range(0,len(Nt_array)):
    tmin = 1 
    plt.xlabel(r"$\tau T$")
    plt.ylabel(r"$R(\tau;\mu_q)$")
    R_mean = np.zeros(Nt_array[i],dtype=float)
    R_stdev = np.zeros(Nt_array[i],dtype=float)
    R_var = np.zeros(Nt_array[i],dtype=float)
    tmp = np.zeros((Nt_array[i],boot_samples),dtype=float)
    for t in range(0,Nt_array[i]):
        for b in range(0,boot_samples):
            tmp[t][b] = float(b_R[i][b*Nt_array[i]+t])
        R_mean[t] = np.mean(tmp[t])
        R_stdev[t] = np.std(tmp[t])
        R_var[t] = np.var(tmp[t])
    R_var_tot += [R_var]
    plt.errorbar(x=np.asarray(corr_t[i][tmin:Nt_array[i]])*T_lat[i],y=R_mean[tmin:Nt_array[i]],yerr=R_stdev[tmin:Nt_array[i]],label=r"$T=$"+str(int(T_array[i]*1000))+"$\mathrm{MeV}$",color=color_list[i],capsize=2)
plt.legend(ncol=2,loc='best')
plt.text(0.5, 0.07, r'$\mu_q = 0$', fontsize=18, ha='center')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/gen2l_R_tau_ratio_plot.pdf",dpi=300)
plt.clf()

    
R_av=[]

for i in range(0,len(Nt_array)):
    R_tmp=[]
    tmin = tminarray[i]-1
    for b in range(0,boot_samples):
        tmp_1=0
        tmp_2=0
        for j in range(tmin,int(Nt_array[i]/2)-1):
            tmp_1 = tmp_1 + b_R[i][b*Nt_array[i]+j]/((R_var_tot[i][j]))
            tmp_2 = tmp_2 + 1/((R_var_tot[i][j]))
        R_tmp += [tmp_1/tmp_2]
    R_av += [np.asarray(R_tmp)]

plt.xlabel(r"$N_{\tau}$")
plt.ylabel(r"$R$ $O(1)$")

R_av_mean = np.zeros(len(Nt_array),dtype=float)
R_av_stdev = np.zeros(len(Nt_array),dtype=float)
for i in range(0,len(Nt_array)):    
    R_av_mean[i] = np.mean(R_av[i])
    R_av_stdev[i] = np.std(R_av[i])
    plt.errorbar(x=Nt_array[i],y=R_av_mean[i],yerr=R_av_stdev[i],color='red',marker='.')
#plt.legend()
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/gen2l_R_ratio_plot.png",dpi=300)
#plt.show()
plt.clf()

##################### O(mu^2) ratio #############################

b_vec_mu=[]
b_axial_mu=[]

for i in range(0,len(Nt_array)):
    b_vec_mu += [np.asarray(corr_total[i])]    
    b_axial_mu += [np.asarray(corr_total_A[i])]
    b_vec_mu[i] = np.asarray(b_vec_mu[i])
    b_axial_mu[i] = np.asarray(b_axial_mu[i])
    tmp1 = np.zeros((Nt_array[i],boot_samples),dtype=float)
    tmp2 = np.zeros((Nt_array[i],boot_samples),dtype=float)
    mean_a = np.zeros(Nt_array[i],dtype=float)
    stdev_a = np.zeros(Nt_array[i],dtype=float)
    mean_v = np.zeros(Nt_array[i],dtype=float)
    stdev_v = np.zeros(Nt_array[i],dtype=float)
    for t in range(0,Nt_array[i]):
        for b in range(0,boot_samples):
            tmp1[t][b] = float(b_axial_mu[i][b*Nt_array[i]+t])
            tmp2[t][b] = float(b_vec_mu[i][b*Nt_array[i]+t])
        mean_a[t] = np.mean(tmp1[t])
        stdev_a[t] = np.std(tmp1[t])
        mean_v[t] = np.mean(tmp2[t])
        stdev_v[t] = np.std(tmp2[t])


b_R_mu=[]
b_R_num_mu=[]
b_R_den_mu=[]

for i in range(0,len(Nt_array)):
        R_num_tmp=[]
        R_den_tmp=[]
        for b in range(0,boot_samples):
            for t in range(0,Nt_array[i]):
                R_num_tmp += [b_vec_mu[i][b*Nt_array[i]+t]/b_vec_mu[i][b*Nt_array[i]+int(Nt_array[i]/2)] - b_axial_mu[i][b*Nt_array[i]+t]/b_axial_mu[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
                R_den_tmp += [b_axial_mu[i][b*Nt_array[i]+t]/b_axial_mu[i][b*Nt_array[i]+int(Nt_array[i]/2)] + b_vec_mu[i][b*Nt_array[i]+t]/b_vec_mu[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
        R_num_tmp = np.asarray(R_num_tmp)
        R_den_tmp = np.asarray(R_den_tmp)
        b_R_num_mu += [R_num_tmp]
        b_R_den_mu += [R_den_tmp]        

#R_num_mu = np.asarray(R_num_mu)
#R_den_mu = np.asarray(R_den_mu)

for i in range(0,len(Nt_array)):    
    b_R_mu += [b_R_num_mu[i]/b_R_den_mu[i]]

mu_mev = mu*1000
R_mu_var_total=np.zeros((len(Nt_array),64),dtype=float)
for i in range(0,len(Nt_array)):
        tmin = 1
        plt.xlabel(r"$\tau T$")
        plt.ylabel(r"$R(\tau,\mu_q)$")
        plt.title('Gen2L')
        R_mu_mean = np.zeros(Nt_array[i],dtype=float)
        R_mu_stdev = np.zeros(Nt_array[i],dtype=float)
        R_mu_var = np.zeros(Nt_array[i],dtype=float)
        tmp = np.zeros((Nt_array[i],boot_samples),dtype=float)
        for t in range(0,Nt_array[i]):
            for b in range(0,boot_samples):
                tmp[t][b] = float(b_R_mu[i][b*Nt_array[i]+t])
            R_mu_mean[t] = np.mean(tmp[t])
            R_mu_stdev[t] = np.std(tmp[t])
        symmetrise(R_mu_mean,Nt_array[i])
        symmetrise(R_mu_stdev,Nt_array[i])
        R_mu_var_total[i] = R_mu_var
        plt.errorbar(x=np.asarray(corr_t[i][tmin:Nt_array[i]])*T_lat[i],y=R_mu_mean[tmin:Nt_array[i]],yerr=R_mu_stdev[tmin:Nt_array[i]],label=r"$T=$"+str(int(T_array[i]*1000))+"$\mathrm{MeV}$",color=color_list[i],capsize=2)
plt.text(0.5, 0.21, r'$\mu_q =$'+str(mu_mev)+' $\mathrm{MeV}$', fontsize=18, ha='center')    
plt.legend(ncol=2,loc='best')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/gen2l_R_mu_tau_ratio_plot_mu_{mu}.pdf",dpi=300)        
plt.clf()

R_av_mu=[]
for i in range(0,len(Nt_array)):
    R_tmp=[]
    tmin = tminarray[i]-1    
    for b in range(0,boot_samples):
        tmp_1=0
        tmp_2=0            
        for j in range(tmin, int(len(R_mu[i])/2)-1):             
            tmp_1 = tmp_1 + b_R_mu[i][b*Nt_array[i]+j]/R_mu[i][j].var
            tmp_2 = tmp_2 + 1/R_mu[i][j].var
        R_tmp += [tmp_1/tmp_2]
    R_av_mu += [np.asarray(R_tmp)]

plt.xlabel(r"$N_{\tau}$")
plt.ylabel(r"$R$ $O(\mu^2)$")
R_av_mu_mean = np.zeros(len(Nt_array),dtype=float)
R_av_mu_stdev = np.zeros(len(Nt_array),dtype=float)
for i in range(0,len(Nt_array)):
    R_av_mu_mean[i] = np.mean(R_av_mu[i])
    R_av_mu_stdev[i] = np.std(R_av_mu[i])
    plt.errorbar(x=Nt_array[i],y=R_av_mu_mean[i],yerr=R_av_mu_stdev[i],color='red',marker='o')
#plt.legend()
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/gen2l_R_mu_ratio_plot_mu_{mu}.png",dpi=300)
#plt.show()
plt.clf()

######################### COMPARISON ######################
plt.xlabel(r"$N_{\tau}$")
plt.ylabel(r"$R$")
R_mean = np.zeros(len(Nt_array),dtype=float)
R_stdev = np.zeros(len(Nt_array),dtype=float)
R_mu_mean = np.zeros(len(Nt_array),dtype=float)
R_mu_stdev = np.zeros(len(Nt_array),dtype=float)
for i in range(0,len(Nt_array)):
    R_mean[i] = np.mean(R_av[i])
    R_stdev[i] = np.std(R_av[i])
    R_mu_mean[i] = np.mean(R_av_mu[i])
    R_mu_stdev[i] = np.std(R_av_mu[i])
    if i==0:
        plt.errorbar(x=Nt_array[i],y=R_mean[i],yerr=R_stdev[i],color='blue',marker='o',label=r'$O(1)$')
        plt.errorbar(x=Nt_array[i],y=R_mu_mean[i],yerr=R_mu_stdev[i],color='red',marker='s',label=r"$O(\mu^2)$")
    else:
        plt.errorbar(x=Nt_array[i],y=R_mean[i],yerr=R_stdev[i],color='blue',marker='o')
        plt.errorbar(x=Nt_array[i],y=R_mu_mean[i],yerr=R_mu_stdev[i],color='red',marker='s')        
plt.legend(loc='best')
plt.axhline(y=0.0,linestyle='--',color='black')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/gen2l_comp_R_plot_mu_{mu}.png",dpi=300)
plt.clf()

datapath = "./interpolate"

if not os.path.exists(datapath):
        os.makedirs(datapath)

g=open("interpolate/boot_gen2l_zero.dat","w")
for i in range(0,len(Nt_array)):
    for b in range(0,boot_samples):
        g.write(str(Nt_array[i])+' '+str(R_av[i][b])+' '+str(len(Nt_array))+'\n')


g=open(f"interpolate/boot_gen2l_full_mu_{mu}.dat","w")
for i in range(0,len(Nt_array)):
    for b in range(0,boot_samples):
        g.write(str(Nt_array[i])+' '+str(R_av_mu[i][b])+' '+str(len(Nt_array))+'\n')

    
