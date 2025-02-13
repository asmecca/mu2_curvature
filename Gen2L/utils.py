import numpy as np
import sys

def correlated_chi_squared(y, y_fit, covariance_matrix):
    """
    Calculate the correlated chi-squared statistic.

    Parameters:
    y (np.ndarray): Observed data points.
    y_fit (np.ndarray): Fitted data points.
    covariance_matrix (np.ndarray): Covariance matrix of the observed data.

    Returns:
    float: The correlated chi-squared value.
    """
    # Calculate the residuals
    residuals = y - y_fit
    
    # Calculate the inverse of the covariance matrix
    covariance_matrix_inv = np.linalg.inv(covariance_matrix)
    
    # Calculate the correlated chi-squared statistic
    chi_squared = np.dot(residuals.T, np.dot(covariance_matrix_inv, residuals))
    
    return chi_squared

def terms_corr_boot(path,channel,Nt,boot_samples):
        term_list = []
        file_terms = path+'/boot/'+channel  #'res.g2.m0.dat'
        ft = open(file_terms,'r')
        for l in ft.readlines():
            x=l.split()
            if x != []:
                term_list += [x]
        ft.close()
        g_list=[]
        g_ii_list = []
        cross_term = []
        disc_term = []
        total_sum=[]
        complete = []
        for b in range(0,boot_samples):
            for t in range(0,Nt):
                g_list += [float(term_list[b*Nt+t][0])]                
                g_ii_list += [float(term_list[b*Nt+t][1]) + float(term_list[b*Nt+t][2]) + float(term_list[b*Nt+t][3])]
                cross_term += [float(term_list[b*Nt+t][4])]
                disc_term += [float(term_list[b*Nt+t][5])]
                total_sum += [float(term_list[b*Nt+t][6])]
                complete += [float(term_list[b*Nt+t][7])]
                
        return g_list,g_ii_list,cross_term,disc_term,total_sum,complete






def get_correlators(path_to_corrs,boot_samples,Nt_array,T_array,a_inv_gev,mu):
    corr_g_zero_A, corr_g_conn_A, corr_cross_A, corr_disc_A, corr_sum_A, corr_total_A = [],[],[],[],[],[]
    corr_g_zero, corr_g_conn, corr_cross, corr_disc, corr_sum, corr_total = [],[],[],[],[],[]
    corr_t, T_lat = [],[]    
    for i in range(0,len(Nt_array)):
        Nt = Nt_array[i]
        T = T_array[i]
        T_lat += [float(float(T)/float(a_inv_gev))]
        path = path_to_corrs+f"/{Nt}x32/analysis_mu_{mu}"        

        g_zero,g_conn,cross_term,disc_term,V_sum,total=terms_corr_boot(path,'res.vector.dat',Nt,boot_samples)
        g_zero_A,g_conn_A,cross_term_A,disc_term_A,A_sum,A_total=terms_corr_boot(path,'res.axial.dat',Nt,boot_samples)

        t=[]
        for b in range(0,boot_samples):
            for j in range(0,Nt):
                t += [j]
        corr_g_zero_A += [g_zero_A]
        corr_g_conn_A += [g_conn_A]
        corr_cross_A += [cross_term_A]
        corr_disc_A += [disc_term_A]
        corr_sum_A += [A_sum]
        corr_g_zero += [g_zero]
        corr_g_conn += [g_conn]
        corr_cross += [cross_term]
        corr_disc += [disc_term]
        corr_sum += [V_sum]
        corr_t+=[t]
        corr_total+=[total]
        corr_total_A += [A_total]
    return corr_g_zero,corr_g_zero_A,corr_total,corr_total_A

def correlated_chi2(x_data, y_data, model, theta, cov_matrix):
    """
    Calculate the correlated chi-square for a set of data points using a given covariance matrix.
    
    Parameters:
    x_data (array): Independent variable data points.
    y_data (array): Measured dependent variable data points.
    model (function): The model function, f(x, theta), to fit the data.
    theta (array): Model parameters.
    cov_matrix (array): Covariance matrix of the measurements.
    
    Returns:
    float: The correlated chi-square value.
    """
    # Invert the covariance matrix
    cov_inv = np.linalg.inv(cov_matrix)
    
    # Calculate the model values for the given parameters
    model_values = model(x_data, *theta)
    
    # Calculate the difference vector (y_data - model(x, theta))
    diff = y_data - model_values
    
    # Calculate the correlated chi-square
    chi2 = np.dot(diff, np.dot(cov_inv, diff))
    
    return chi2


# Defining a function for generating bootstrap samples from a normal distribution
def bootstrap_normal(mu, sigma, n_samples, n_bootstrap):
    """
    Generate bootstrap samples from a normal distribution with given mean (mu)
    and standard deviation (sigma), and plot the distribution of sample means.
    
    Parameters:
    mu (float): Mean of the normal distribution
    sigma (float): Standard deviation of the normal distribution
    n_samples (int): Number of synthetic data samples to generate
    n_bootstrap (int): Number of bootstrap iterations
    
    Returns:
    None: Plots the bootstrap means distribution
    """
    # Generate synthetic data from the normal distribution
    synthetic_data = np.random.normal(mu, sigma, n_samples)
    
    # Generate bootstrap samples
    bootstrap_means = []
    for _ in range(n_bootstrap):
        bootstrap_sample = np.random.choice(synthetic_data, size=n_samples, replace=True)
        bootstrap_means.append(np.mean(bootstrap_sample))
    return bootstrap_means
    # Plotting the bootstrap means distribution
    #plt.hist(bootstrap_means, bins=30, density=True, alpha=0.7, color='blue')
    #plt.axvline(mu, color='red', linestyle='--', label=f"True Mean = {mu}")
    #plt.title(f"Bootstrap Distribution of Sample Means (n = {n_bootstrap})")
    #plt.xlabel('Mean')
    #plt.ylabel('Density')
    #plt.legend()
    #plt.show()



def compute_covariance_matrix(bootstrap_samples):
    """
    Computes the covariance matrix from a list of bootstrap sample sets.
    
    Parameters:
    bootstrap_samples (list of arrays or 2D numpy array): A list where each element is an array of bootstrap samples.
                                                          If using a numpy array, it should have shape (6, n_samples).
                                                          
    Returns:
    covariance_matrix (numpy array): The 6x6 covariance matrix.
    """
    # Convert the input to a numpy array if it's not already
    bootstrap_samples = np.array(bootstrap_samples)
    
    # Calculate the covariance matrix
    covariance_matrix = np.cov(bootstrap_samples)
    
    return covariance_matrix

