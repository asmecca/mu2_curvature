import numpy as np
import sys

def correlated_chi2(x, y, cov_matrix, model_func, popt):
    """
    Calculate the chi-squared value for a correlated fit.

    Parameters:
        x (numpy.ndarray): Independent variable data.
        y (numpy.ndarray): Dependent variable data.
        cov_matrix (numpy.ndarray): Covariance matrix for the dependent variable.
        model_func (callable): The model function used for fitting.
        popt (list): Optimized parameters from the fit.

    Returns:
        float: The correlated chi-squared value.
    """
    residuals = y - model_func(x, *popt)
    cov_inv = np.linalg.inv(cov_matrix)  # Inverse of the covariance matrix
    chi2 = residuals.T @ cov_inv @ residuals
    return chi2

def naive_correlated_chi2(x_data, y_data, model, theta, cov_matrix):
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

