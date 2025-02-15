#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.optimize import curve_fit
from scipy import special

from utils import correlated_chi2
from utils import bootstrap_normal
from utils import compute_covariance_matrix
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('font',**{'size':14})
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

def calculate_correlation_matrix(data):
    """
    Calculate the correlation matrix for a given dataset.
    
    Parameters:
        data (numpy.ndarray): A 2D array where each row represents a dataset 
                              corresponding to a different lattice spacing, 
                              and each column is a specific measurement.
                              
    Returns:
        numpy.ndarray: The correlation matrix.
    """
    return np.corrcoef(data)

def check_significant_correlations(corr_matrix, threshold=0.5):
    """
    Identify pairs of variables with significant correlations.
    
    Parameters:
        corr_matrix (numpy.ndarray): The correlation matrix.
        threshold (float): The correlation coefficient threshold for significance.
        
    Returns:
        list of tuple: Pairs of indices with significant correlations.
    """
    significant_pairs = []
    n = corr_matrix.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(corr_matrix[i, j]) > threshold:
                significant_pairs.append((i, j, corr_matrix[i, j]))
    return significant_pairs

def plot_correlation_matrix(corr_matrix, labels=None):
    """
    Plot the correlation matrix as a heatmap.
    
    Parameters:
        corr_matrix (numpy.ndarray): The correlation matrix to plot.
        labels (list of str, optional): Labels for the axes. Default is None.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(corr_matrix, cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label="Correlation Coefficient")
    plt.title("Correlation Matrix Heatmap")
    
    if labels:
        plt.xticks(range(len(labels)), labels, rotation=45)
        plt.yticks(range(len(labels)), labels)
    
    plt.tight_layout()
    plt.show()

def fit_func(x,k):
    return T0*(1-k*(x/(T0**2)))

def fit_func_2L(x,k):
    return T0_2L*(1-k*(x/(T0_2L**2)))

# Other groups results
Budapest = [0.0153,0.0018] #0.0153(18)
Pisa = [0.0145,0.0025] #0.0145(25)
Pisa_old = [0.013,0.002] #0.013(2)(1)
Bazavov = [0.012, 0.004] #0.012(4)
Bellwied = [0.0149,0.0021]

# Pseudocritical temperatures obtained using the chiral condensate for each ensemble
T_chiral = [182,2] #182(2)
T_chiral_2L = [167,3] #167(3)

# Getting bootstrap data for averaged R(tau,\mu^2) for each ensemble
filename= f"gen2/interpolate/boot_pseudo_T.dat"
filename_2L= f"Gen2L/interpolate/boot_gen2l_pseudo_T.dat"
f=open(filename,'r')

mu,T,T_minus,T_plus=[],[],[],[]
mu_2L,T_2L,T_minus_2L,T_plus_2L=[],[],[],[]
for line in f.readlines():
    if line[0] != '#':
        x=line.split()
        mu += [float(x[0])]
        T += [float(x[1])]
        boot_samples = int(x[2])
mu_q = mu
mu = np.asarray(mu)*1000*3
T = np.asarray(T)*1000

len_mu = int(len(mu)/boot_samples)
f.close()

f=open(filename_2L,'r')
for line in f.readlines():
    if line[0] != '#':
        x=line.split()
        mu_2L += [float(x[0])]
        T_2L += [float(x[1])]
        boot_samples = int(x[2])
mu_q_2L = mu_2L
mu_2L = np.asarray(mu_2L)*1000*3
T_2L = np.asarray(T_2L)*1000

len_mu_2L = int(len(mu_2L)/boot_samples)
f.close()


# Computing mean and standard deviation for each ensemble
mu_boot = np.zeros((len_mu,boot_samples),dtype=float)
mu_q_boot = np.zeros((len_mu,boot_samples),dtype=float)
T_boot = np.zeros((len_mu,boot_samples),dtype=float)
T_mean = np.zeros(len_mu,dtype=float)
T_stdev = np.zeros(len_mu,dtype=float)
for i in range(0,len_mu):
    for b in range(0,boot_samples):
        mu_q_boot[i][b] = mu_q[i*boot_samples+b]
        mu_boot[i][b] = mu[i*boot_samples+b]
        T_boot[i][b] = T[i*boot_samples+b]
    T_mean[i] = np.mean(T_boot[i])
    T_stdev[i] = np.std(T_boot[i])
x= np.arange(mu[0],mu[len(mu)-1],10)
mu_boot = mu_boot**2

print('Gen2 - This is T mean: ',T_mean)
print('Gen2 - This is T err: ',T_stdev)

mu_boot_2L = np.zeros((len_mu_2L,boot_samples),dtype=float)
mu_q_boot_2L = np.zeros((len_mu_2L,boot_samples),dtype=float)
T_boot_2L = np.zeros((len_mu_2L,boot_samples),dtype=float)
T_mean_2L = np.zeros(len_mu_2L,dtype=float)
T_stdev_2L = np.zeros(len_mu_2L,dtype=float)
for i in range(0,len_mu_2L):
    for b in range(0,boot_samples):
        mu_q_boot_2L[i][b] = mu_q_2L[i*boot_samples+b]
        mu_boot_2L[i][b] = mu_2L[i*boot_samples+b]
        T_boot_2L[i][b] = T_2L[i*boot_samples+b]
    T_mean_2L[i] = np.mean(T_boot_2L[i])
    T_stdev_2L[i] = np.std(T_boot_2L[i])
x_2L= np.arange(mu_2L[0],mu_2L[len(mu_2L)-1],10)
mu_boot_2L = mu_boot_2L**2

print('Gen2L - This is T mean: ',T_mean_2L)
print('Gen2L - This is T err: ',T_stdev_2L)

mu_arr = np.zeros(len_mu,dtype=float)
mu_q_arr = np.zeros(len_mu,dtype=float)
for i in range(0,len_mu):
    mu_arr[i] = mu_boot[i][0]
    mu_q_arr[i] = mu_q_boot[i][0]

mu_arr_2L = np.zeros(len_mu_2L,dtype=float)
mu_q_arr_2L = np.zeros(len_mu_2L,dtype=float)
for i in range(0,len_mu_2L):
    mu_arr_2L[i] = mu_boot_2L[i][0]
    mu_q_arr_2L[i] = mu_q_boot_2L[i][0]
    

# removing the last three points from Generation 2 array since they are not included in the fit
T_boot_cut = np.zeros((len_mu-3,boot_samples),dtype=float)
for i in range(0,len_mu-3):
    for b in range(0,boot_samples):
        T_boot_cut[i][b] = T_boot[i][b]

# Computing covariance and correlation matrices - Gen2
cov_matrix = np.cov(T_boot_cut)
correlation_matrix = calculate_correlation_matrix(T_boot_cut)
plot_correlation_matrix(correlation_matrix)
plot_correlation_matrix(cov_matrix)

# Computing covariance and correlation matrices - Gen2L
cov_matrix_2L = np.cov(T_boot_2L)
correlation_matrix_2L = calculate_correlation_matrix(T_boot_2L)
plot_correlation_matrix(correlation_matrix_2L)
plot_correlation_matrix(cov_matrix_2L)

# Computing standard error on pseudocritical temperatures for different mu^2 and for each ensemble
T_err = np.sqrt(np.diag(cov_matrix))
T_err_2L = np.sqrt(np.diag(cov_matrix_2L))

T0 = T_mean[0]
T0_2L = T_mean_2L[0]

# fitting the Generation 2 data
p0 = [1.]
popt, pcov = curve_fit(fit_func, mu_arr[:(len_mu-3)], T_mean[:int(len_mu-3)], sigma=cov_matrix, absolute_sigma=True,p0=p0)

# fitting the Generation 2L data
p0_2L = [1.]
popt_2L, pcov_2L = curve_fit(fit_func_2L, mu_arr_2L[:(len_mu_2L)], T_mean_2L[:int(len_mu_2L)], sigma=cov_matrix_2L, absolute_sigma=True,p0=p0_2L)

# Finding the standard error of kappa
perr = np.sqrt(np.diag(pcov))
perr_2L = np.sqrt(np.diag(pcov_2L))

print('Gen2: ',popt[0])
print('Gen2L: ',popt_2L[0])

# Calculating correlated chi^2 
corr_chi2 = correlated_chi2(mu_arr[:(len_mu-3)], T_mean[:(len_mu-3)], cov_matrix, fit_func, popt)
corr_chi2_red = corr_chi2/(len(mu_arr[:(len_mu-3)])-1)
print('corr chi2 d.o.f = ',corr_chi2_red)
corr_chi2_2L = correlated_chi2(mu_arr_2L[:(len_mu_2L)], T_mean_2L[:(len_mu_2L)], cov_matrix_2L, fit_func_2L, popt)
corr_chi2_red_2L = corr_chi2_2L/(len(mu_arr[:(len_mu_2L)])-1)
print('2L corr chi2 d.o.r. = ',corr_chi2_red_2L)



# Plotting data with fitlines
popt_plus = popt + perr
popt_minus = popt - perr
popt_plus_2L = popt_2L + perr_2L
popt_minus_2L = popt_2L - perr_2L

plotpath = '.'

mu_q = np.asarray(mu_q_arr)
mu_q = mu_q*1000*3
mu_q_2L = np.asarray(mu_q_arr_2L)
mu_q_2L = mu_q_2L*1000*3

plt.xlabel(r'$\mu_B$ $\mathrm{[MeV]}$')
plt.ylabel(r'$T_{\rm pc}(\mu_B)$ $\mathrm{[MeV]}$')

plt.errorbar(x=-10,y=T_chiral[0],yerr=T_chiral[1],fmt='*',color='darkgray',capsize=3)
plt.errorbar(x=-10,y=T_chiral_2L[0],yerr=T_chiral_2L[1],fmt='x',color='black',capsize=3)
plt.errorbar(x=np.sqrt(mu_arr),y=T_mean,yerr=T_stdev,fmt='o',capsize=2,label='Generation 2')
plt.errorbar(x=np.sqrt(mu_arr_2L),y=T_mean_2L,yerr=T_stdev_2L,fmt='s',capsize=2,label='Generation 2L')
plt.plot(x,fit_func(x**2,*popt),color='#1f77b4')
plt.plot(x_2L,fit_func_2L(x_2L**2,*popt_2L),color='#ff7f0e')

plt.fill_between(x,fit_func(x**2,*popt_plus),fit_func(x**2,*popt_minus),alpha=0.3,color='#1f77b4')
plt.fill_between(x_2L,fit_func_2L(x_2L**2,*popt_plus_2L),fit_func_2L(x_2L**2,*popt_minus_2L),alpha=0.3,color='#ff7f0e')
plt.plot(mu_q[5:7],mu_q[5:7]/3,linestyle='--',color='black',label=r'$\mu_q/T=1$')

plt.legend(loc='best')
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/boot_fit_final_mu_method1.pdf",dpi=300)
plt.show()


# Plotting results for kappa
res = [popt[0],perr[0]]
res_2L = [popt_2L[0],perr_2L[0]]

plt.ylim(0.005,0.050,0.005)
plt.ylabel(r'$\kappa$')
line1 = plt.errorbar(x=1,y=res[0],yerr=res[1],fmt='o',capsize=2)
line2 = plt.errorbar(x=1.02,y=res_2L[0],yerr=res_2L[1],fmt='s',capsize=2)
plt.axhspan(res[0]-res[1], res[0]+res[1], alpha=0.2, color='#1f77b4')
plt.axhspan(res_2L[0]-res_2L[1], res_2L[0]+res_2L[1], alpha=0.2, color='#ff7f0e')
plt.axvline(x=1.05,linestyle='--',color = 'darkgray')
line3 = plt.errorbar(x=1,y=0.0,yerr=res[1],fmt='o',capsize=2,color='white')
line4 = plt.errorbar(x=1,y=0.0,yerr=res[1],fmt='o',capsize=2,color='white')
line5 = plt.errorbar(x=1,y=0.0,yerr=res[1],fmt='o',capsize=2,color='white')
line6 = plt.errorbar(x=1.2,y=Budapest[0],yerr=Budapest[1],fmt='p',capsize=2)

line7 = plt.errorbar(x=1.25,y=Bazavov[0],yerr=Bazavov[1],fmt='>',capsize=2)
line8 = plt.errorbar(x=1.3,y=Pisa[0],yerr=Pisa[1],fmt='<',capsize=2)
line9 = plt.errorbar(x=1.35,y=Bellwied[0],yerr=Bellwied[1],fmt='d',capsize=2)
line10 = plt.errorbar(x=1.4,y=Pisa_old[0],yerr=Pisa_old[1],fmt='^',capsize=2)
#line11 = plt.errorbar(x=1.45,y=0.020,yerr=0.004,fmt='*',capsize=2)

handles = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10]  # `None` adds a blank space
labels = ["Generation 2", "Generation 2L", "", "", "", r"Borsanyi $et$ $al.$ 2020", "HotQCD 2018", "Bonati $et$ $al.$ 2018",  "Bellwied $et$ $al.$ 2015", "Bonati $et$ $al.$ 2015"] #1508.07599
plt.xticks([])
plt.legend(handles=handles,
    labels=labels, loc='best',ncol=2)
fig=plt.gcf()
fig.savefig(plotpath+f"/plots/final_kappa.pdf",dpi=300)
plt.show()

print('Gen2 - This is kappa: ',res)
print('Gen2L - This is kappa: ',res_2L)
