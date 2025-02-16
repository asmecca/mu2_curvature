#!/usr/bin/env python3

# The iteration is: for each s_id there are all the timeslices, so for bin=1 we have from t=1 to t=Nt

import numpy as np
import pyerrors as pe
import matplotlib.pyplot as plt
import os
import gvar as gv
import sys
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utils import get_correlators_gen2

from cycler import cycler

def get_corr(filename,tag,Nt):
    lookup = tag #' G0.S0'
    numb_lines=[]
    corr=[]
    with open(filename) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                numb_lines += [int(num)]
    end_num=numb_lines[0]+Nt
    f=open(filename,'r')
    content=f.readlines()
    G0=content[(numb_lines[0]+1):end_num]
    corr=np.zeros(Nt-1,dtype=float)
    time=np.zeros(Nt-1,dtype=float)    
    for i in range(0,len(G0)):        
        time[i] = float(G0[i].split()[0])
        corr[i] = float(G0[i].split()[1])
    return time,corr

plt.rc('font',**{'size':14})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

#default_cycler = (cycler(marker=['v','<','^','>','o','s','p'])+cycler(color=['darkred','red','orange','gold','lightblue','cyan','blue']))
default_cycler = (cycler(marker=['v','s','^','*','o','<','p'])+cycler(color=['darkred','red','orange','gold','lightblue','cyan','blue']))
plt.rc('axes', prop_cycle=default_cycler)
plt.rc('lines', linestyle='')

marker_list=['v','<','^','>','o','s','p']
color_list=['darkred','red','orange','gold','#71daeb','cyan','blue']
a_inv_gev = 5.63
T_array = [0.352, 0.281, 0.235, 0.201, 0.176, 0.156, 0.141]
tminarray = [3, 4, 4, 5, 6, 7, 8]#, 8]
Nt_array =[16, 20, 24, 28, 32, 36, 40]

mu = 0.028 #float(sys.argv[1])

freefile = "./plot-GV-GA-24x36-32-28-24-20-16_7.agr"

path_to_corr = "DATA/gen2" # Path to directory with name NtxNs containing the correlators

# Gets the correlators for the different contributions and channels
# corr_g_zero -> O(1) connected correlator, vector channel
# corr_g_zero_A -> O(1) connected correlator, axial-vector channel
# corr_total -> Correlator containing the sum of O(1) and O(mu^2) contributions at a certain value of mu, vector channel
# corr_total_A -> Correlator containing the sum of O(1) and O(mu^2) contributions at a certain value of mu, axial-vector channel
# boot_samples -> number of bootstrap samples used in the analysis

corr_t,corr_g_zero,corr_g_zero_A,corr_total,corr_total_A,boot_samples = get_correlators_gen2(path_to_corr,Nt_array,T_array,a_inv_gev,mu)

plotpath = '..'
if not os.path.exists(plotpath+'/plots'):
        os.makedirs(plotpath+'/plots')

gv_vec=[]
gv_axial=[]

b_vec=[]
b_axial=[]

for i in range(0,len(Nt_array)):
        b_vec += [np.asarray(corr_g_zero[i])]
        b_axial += [np.asarray(corr_g_zero_A[i])]
        b_vec[i] = np.asarray(b_vec[i])
        b_axial[i] = np.asarray(b_axial[i])
        tmp1 = np.zeros((Nt_array[i],2000),dtype=object)
        tmp2 = np.zeros((Nt_array[i],2000),dtype=object)
        mean_a = np.zeros(Nt_array[i],dtype=object)
        stdev_a = np.zeros(Nt_array[i],dtype=object)
        mean_v = np.zeros(Nt_array[i],dtype=object)
        stdev_v = np.zeros(Nt_array[i],dtype=object)
        for t in range(0,Nt_array[i]):
            for b in range(0,2000):
                tmp1[t][b] = float(b_axial[i][b*Nt_array[i]+t])
                tmp2[t][b] = float(b_vec[i][b*Nt_array[i]+t])
            mean_a[t] = np.mean(tmp1[t])
            stdev_a[t] = np.std(tmp1[t])
            mean_v[t] = np.mean(tmp2[t])
            stdev_v[t] = np.std(tmp2[t])
        gv_vec += [gv.gvar(mean_v,stdev_v)]
        gv_axial += [gv.gvar(mean_a,stdev_a)]
        gv_vec[i] = np.asarray(gv_vec[i])
        gv_axial[i] = np.asarray(gv_axial[i])        


# Getting free correlator
tag_arr = [' G0.S0',' G0.S1',' G0.S2', 'G0.S3', 'G0.S4', 'G0.S5']
time=np.zeros(len(tag_arr),dtype=object)
corr=np.zeros(len(tag_arr),dtype=object)
for i in range(0,len(Nt_array)):
    if i < len(tag_arr):        
        tmp_time,tmp_corr = get_corr(freefile,tag_arr[i],Nt_array[i])
        time[i] = tmp_time
        corr[i] = tmp_corr

gv_R=[]
gv_R_num=[]
gv_R_den=[]

for i in range(0,len(Nt_array)):
        gv_R_num += [gv_vec[i]/gv_vec[i][int(Nt_array[i]/2)] - gv_axial[i]/gv_axial[i][int(Nt_array[i]/2)] ]
        gv_R_den += [gv_axial[i]/gv_axial[i][int(Nt_array[i]/2)] + gv_vec[i]/gv_vec[i][int(Nt_array[i]/2)] ]

for i in range(0,len(Nt_array)):
    gv_R += [gv_R_num[i]/gv_R_den[i]]

b_R_naive=[]
b_R_num_naive=[]
b_R_den_naive=[]

for i in range(0,len(Nt_array)):
        R_num_tmp=[]
        R_den_tmp=[]
        for b in range(0,2000):
            for t in range(0,Nt_array[i]):
                R_num_tmp += [b_vec[i][b*Nt_array[i]+t]/b_vec[i][b*Nt_array[i]+int(Nt_array[i]/2)] - b_axial[i][b*Nt_array[i]+t]/b_axial[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
                R_den_tmp += [b_axial[i][b*Nt_array[i]+t]/b_axial[i][b*Nt_array[i]+int(Nt_array[i]/2)] + b_vec[i][b*Nt_array[i]+t]/b_vec[i][b*Nt_array[i]+int(Nt_array[i]/2)] ]
        R_num_tmp = np.asarray(R_num_tmp)
        R_den_tmp = np.asarray(R_den_tmp)
        b_R_num_naive += [R_num_tmp]
        b_R_den_naive += [R_den_tmp]        


for i in range(0,len(Nt_array)):    
    b_R_naive += [b_R_num_naive[i]/b_R_den_naive[i]]

for i in range(0,len(Nt_array)):
    if i < len(tag_arr)-2:
        tmin = 1
        tmax = tmaxarray[i]        
        plt.xlabel(r"$\tau T$",fontsize=16)
        plt.ylabel(r"$R(\tau,0)$",fontsize=16)
        R_naive_mean = np.zeros(Nt_array[i],dtype=float)
        R_naive_stdev = np.zeros(Nt_array[i],dtype=float)
        tmp = np.zeros((Nt_array[i],2000),dtype=float)
        for t in range(0,Nt_array[i]):
            for b in range(0,2000):
                tmp[t][b] = float(b_R_naive[i][b*Nt_array[i]+t])
            R_naive_mean[t] = np.mean(tmp[t])
            R_naive_stdev[t] = np.std(tmp[t])
        plt.errorbar(x=np.asarray(corr_t[i][tmin:Nt_array[i]])*T_lat[i],y=R_naive_mean[tmin:Nt_array[i]],yerr=R_naive_stdev[tmin:Nt_array[i]],label=r"$T=$"+str(int(T_array[i]*1000))+"$\mathrm{MeV}$",color=color_list[i],capsize=3,markersize=8)

for i in range(len(Nt_array)):
    if i < len(tag_arr)-2:
        plt.plot(time[i],corr[i],linestyle='--',marker='.',markersize=1,color=color_list[i])
plt.text(0.4,0.2,r'Generation 2',fontsize=16)
plt.legend(ncol=2,loc='best')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/free_R_tau_ratio.pdf",dpi=300)
#pdf.savefig( fig )
#plt.show()
plt.clf()

