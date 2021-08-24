# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 14:04:36 2021

@author: LuisAlvaro
"""

# matplotlib inline
import scipy
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

def plot_contours(data, means, covs, title):
    """visualize the gaussian components over the data"""
    plt.figure()
    plt.plot(data[:, 0], data[:, 1], 'ko')

    k = means.shape[0]
    
    maxX = max(data[:,0])*1.1
    minX = min(data[:,0])*1.1
    maxY = max(data[:,1])*1.1
    minY = min(data[:,1])*1.1
    
    N_bins = 1000  # size of grid
    deltaX = (maxX-minX)/N_bins
    deltaY = (maxY-minY)/N_bins
    
    x = np.arange(minX*1.1, maxX*1.1, deltaX)
    y = np.arange(minY*1.1, maxY*1.1, deltaY)
    x_grid, y_grid = np.meshgrid(x, y)
    coordinates = np.array([x_grid.ravel(), y_grid.ravel()]).T

    col = ['green', 'red', 'indigo']
    for i in range(k):
        mean = means[:,i]
        cov = covs[i]
        z_grid = multivariate_normal(mean, cov).pdf(coordinates).reshape(x_grid.shape)
        plt.contour(x_grid, y_grid, z_grid, colors = col[i])

    plt.title(title)
    plt.tight_layout()
    plt.plot(means[0,0],means[1,0],marker='X',color='r', markersize=12)
    plt.plot(means[0,1],means[1,1],marker='o', color='y', markersize=12)

    plt.show()


def cost(data, R, Mu):
    N, D = data.shape
    K = Mu.shape[1]
    J = 0
    for k in range(K):
        J += np.dot(np.linalg.norm(data - np.array([Mu[:, k], ] * N), axis=1)**2, R[:, k])
    return J

# TODO: K-Means Assignment Step
def km_assignment_step(data, Mu):
    """ Compute K-Means assignment step
    
    Args:
        data: a NxD matrix for the data points
        Mu: a DxK matrix for the cluster means locations
    
    Returns:
        R_new: a NxK matrix of responsibilities
    """
    
    # Fill this in:
    N, D = data.shape # Number of datapoints and dimension of datapoint
    K = Mu.shape[1] # number of clusters
    r = np.zeros([N,K])

    for k in range(K):
        r[:, k] = np.linalg.norm(data - np.array([Mu[:, k], ] * N), axis=1)
    
    arg_min = r.argmin(axis=1) # argmax/argmin along dimension 1
    R_new = np.eye(K)[arg_min] # Get Cluster Assignment

    return R_new

# TODO: K-means Refitting Step
def km_refitting_step(data, R, Mu, plotstatus):
    """ Compute K-Means refitting step.
    
    Args:
        data: a NxD matrix for the data points
        R: a NxK matrix of responsibilities
        Mu: a DxK matrix for the cluster means locations
    
    Returns:
        Mu_new: a DxK matrix for the new cluster means locations
    """
    N, D = data.shape # Number of datapoints and dimension of datapoint
    K = Mu.shape[1]  # number of clusters
    Mu_new = (R.T.dot(data)/np.sum(R,axis=0)).T
    
    if (plotstatus):
        plt.scatter(data[:,0],data[:,1],c=R.argmax(axis=1))
        plt.plot(Mu_new[0,0],Mu_new[1,0],marker='X',color='r', markersize=12)
        plt.plot(Mu_new[0,1],Mu_new[1,1],marker='o', color='b', markersize=12)
        plt.show()
        print(Mu_new)
        print(np.sum(R,axis=0))
    return Mu_new

def normal_density(x, mu, Sigma):
    return np.exp(-.5 * np.dot(x - mu, np.linalg.solve(Sigma, x - mu))) \
        / np.sqrt(np.linalg.det(2 * np.pi * Sigma))
        
        
def log_likelihood(data, Mu, Sigma, Pi):
    """ Compute log likelihood on the data given the Gaussian Mixture Parameters.
    
    Args:
        data: a NxD matrix for the data points
        Mu: a DxK matrix for the means of the K Gaussian Mixtures
        Sigma: a list of size K with each element being DxD covariance matrix
        Pi: a vector of size K for the mixing coefficients
    
    Returns:
        L: a scalar denoting the log likelihood of the data given the Gaussian Mixture
    """
    # Fill this in:
    N, D = data.shape  # Number of datapoints and dimension of datapoint
    K = Mu.shape[1]    # number of mixtures
    L, T = 0., 0.
    for n in range(N):
        
        for k in range(K):
            T += Pi[k]*normal_density(data[n,:],Mu[:,k],Sigma[k]) # Compute the likelihood from the k-th Gaussian weighted by the mixing coefficients 
        L += np.log(T)
    return L

# TODO: Gaussian Mixture Maximization Step
def gm_m_step(data, Gamma, plotstatus):
    """ Gaussian Mixture Maximization Step.

    Args:
        data: a NxD matrix for the data points
        Gamma: a NxK matrix of responsibilities 
    
    Returns:
        Mu: a DxK matrix for the means of the K Gaussian Mixtures
        Sigma: a list of size K with each element being DxD covariance matrix
        Pi: a vector of size K for the mixing coefficients
    """
    # Fill this in:
    N, D = data.shape  # Number of datapoints and dimension of datapoint
    K = Gamma.shape[1]   # number of mixtures
    Nk = np.sum(Gamma, axis=0) # Sum along first axis 
    Mu = (1/Nk)*(Gamma.T.dot(data)).T
    Sigma = [np.eye(2), np.eye(2)]
    
    for k in range(K):
        weightedSum = np.zeros([D,D])
        diff = (data-np.array([Mu[:, k], ]* N)).T
        # weightedSum = np.sum(Gamma,axis=0)*(diff.dot(diff.T))
        weightedSum = (Gamma[:,k]*diff).dot(diff.T)
        Sigma[k] = weightedSum/Nk[k]
        
    Pi = Nk/N 
    if (plotstatus):
        plt.scatter(data[:,0],data[:,1],c=Gamma.argmax(axis=1))
        plt.plot(Mu[0,0],Mu[1,0],marker='X',color='r', markersize=12)
        plt.plot(Mu[0,1],Mu[1,1],marker='o', color='b', markersize=12)
        plt.show()
        print(Mu)
        print(np.sum(Gamma,axis=0))
    
    return Mu, Sigma, Pi



# TODO: Gaussian Mixture Expectation Step
def gm_e_step(data, Mu, Sigma, Pi):
    """ Gaussian Mixture Expectation Step.

    Args:
        data: a NxD matrix for the data points
        Mu: a DxK matrix for the means of the K Gaussian Mixtures
        Sigma: a list of size K with each element being DxD covariance matrix
        Pi: a vector of size K for the mixing coefficients
    
    Returns:
        Gamma: a NxK matrix of responsibilities 
    """
    # Fill this in:
    N, D = data.shape  # Number of datapoints and dimension of datapoint
    K = Mu.shape[1]   # number of mixtures
    Gamma = np.zeros([N,K]) # zeros of shape (N,K), matrix of responsibilities
    for n in range(N):
        for k in range(K):
            Dens = normal_density(data[n,:],Mu[:,k], Sigma[k])
            Gamma[n, k] = Pi[k]*Dens
        Gamma[n, :] /= np.sum(Gamma[n,:]) # Normalize by sum across second dimension (mixtures)
    return Gamma

def KMeansAlgorithm(data, labels, Mu0, PrtError):
    
    print("\n\n **** K-Means Algorithm ****\n\n")
    
    # Size
    num_samples = data.shape[0]
    
    # Plot Generated Data - Item 1(a)
    print("\n\n ---- Original Data ----\n\n")
    plt.scatter(data[:,0],data[:,1],c=labels)
    plt.show()
    
    # Call K-Mean Procedure
    # TODO: Run this cell to call the K-means algorithm
    N, D = data.shape
    K = 2
    max_iter = 100
    class_init = np.random.binomial(1., .5, size=N)
    R = np.vstack([class_init, 1 - class_init]).T

    Mu = Mu0
    
    R.T.dot(data), np.sum(R, axis=0)

    # Changed to stop after convergence
    it = 0

    tolerance = 1e-5                  # Tolerance to check convergence
    Converged = False                 # Set Control Variable
    cost_history = np.zeros(max_iter)  # Cost History of Convergence

    while (not Converged and (it < max_iter)):
        R = km_assignment_step(data, Mu)
        Mu = km_refitting_step(data, R, Mu, False) # DEBUG -> (it % 10 == 0))
    
        cost_history[it] = cost(data, R, Mu)
        if (it > 0):
            Converged = (abs(cost_history[it]-cost_history[it-1])<=tolerance)
    
        # print(it, cost_history[it])
    
        # Increment iteration
        it += 1

    print("\n>>> %s after %d iterations with Cost=%.5f" % 
          (("Convergence" if Converged else "Not converged"), 
           it, cost_history[it if Converged else (it-1)]))


    class_1 = np.where(R[:, 0])
    class_2 = np.where(R[:, 1])
    
    print("\n\nTotal No.of Iterations %d" % it)
    print("\nK-Means classified as Class 1 %d points and as Class 2 %d points\n" % (len(class_1[0]), len(class_2[0])))
    
    plt.scatter(data[:,0],data[:,1],c=R.argmax(axis=1))
    plt.plot(Mu[0,0],Mu[1,0],marker='X',color='r', markersize=12)
    plt.plot(Mu[0,1],Mu[1,1],marker='o', color='b', markersize=12)
    plt.show()
    
    # Also plot the Cost vs. iterations
    plt.figure()
    plt.title("Cost vs Iterations")
    plt.ylabel("Cost")
    plt.xlabel("Iterations")
    plt.plot(range(it), cost_history[:it])
    plt.show()

    accuracy = (1-np.count_nonzero(R.argmax(axis=1)-labels)/num_samples)*100.0
    print("\nK-Means accuracy is %5.2f%%" % accuracy)
    
    if (PrtError):
        print("\n\n------ Listing of Classification Errors ------\n")

        Err = np.nonzero(R.argmax(axis=1)-labels)[0]
        print("\n>>> Total No. of Errors: %6d\n\n" % Err.shape[0])
        dErr = {'Item No.':Err,
                'Label':labels[Err].astype(int),
                'Classif.':R.argmax(axis=1)[Err]}
        dfErr = pd.DataFrame(dErr)
    
        print(dfErr.to_string(header=True, index=False))
    
    
def EMAlgorithm(data, labels, Mu0, PrtError):
    
    print("\n\n **** E.M. Algorithm ****\n\n")

    # Size
    num_samples = data.shape[0]

    # Plot Generated Data - Item 1(a)
    print("\n\n ---- Original Data ----\n\n")
    plt.scatter(data[:,0],data[:,1],c=labels)
    plt.show()
    
    # TODO: Run this cell to call the Gaussian Mixture EM algorithm
    N, D = data.shape
    K = 2
    
    Mu = Mu0

    Sigma = [np.eye(2), np.eye(2)]
    Pi = np.ones(K) / K
    Gamma = np.zeros([N, K]) # Gamma is the matrix of responsibilities 

    max_iter  = 200
    
    # Changed to stop after convergence
    it = 0

    tolerance = 1e-5                     # Tolerance to check convergence
    Converged = False                    # Set Control Variable
    loglik_history = np.zeros(max_iter)  # Log-Likelihood History of Convergence

    while (not Converged and (it < max_iter)):
        Gamma = gm_e_step(data, Mu, Sigma, Pi)
        Mu, Sigma, Pi = gm_m_step(data, Gamma, False) # DEBUG - (it % 10 == 0))
        
        loglik_history[it] = log_likelihood(data, Mu, Sigma, Pi)
        if (it > 0):
            Converged = (abs(loglik_history[it]-loglik_history[it-1])<=tolerance)
    
        # print(it, loglik_history[it])
    
        # Increment iteration
        it += 1
    
    print("\n>>> %s after %d iterations with LogLike=%.5f" % 
        (("Convergence" if Converged else "Not converged"), 
         it, loglik_history[it if Converged else (it-1)]))
    
    class_1 = np.where(Gamma[:, 0] >= .5)
    class_2 = np.where(Gamma[:, 1] >= .5)
    
    print("\nE.M. classified as Class 1 %d points and as Class 2 %d points\n" % (len(class_1[0]), len(class_2[0])))

    plt.figure()
    plt.title("Log-Likelihood vs Iterations")
    plt.ylabel("Log-Likelihood")
    plt.xlabel("Iterations")
    plt.plot(range(it), loglik_history[:it])
    plt.show()


    accuracy = (1-np.count_nonzero(Gamma.argmax(axis=1)-labels)/num_samples)*100.0
    print("\nE.M. accuracy is %5.2f%%" % accuracy)
    
    if (PrtError):
        print("\n\n------ Listing of Classification Errors ------\n")

        Err = np.nonzero(Gamma.argmax(axis=1)-labels)[0]
        print("\n>>> Total No. of Errors: %6d\n\n" % Err.shape[0])
        dErr = {'Item No.':Err,
                'Label':labels[Err].astype(int),
                'Classif.':Gamma.argmax(axis=1)[Err]}
        dfErr = pd.DataFrame(dErr)
    
        print(dfErr.to_string(header=True, index=False))
    
    plot_contours(data, Mu, Sigma, "Contour-Plot - E.M.")


def main():
    
    # TODO: Run this cell to generate the data
    num_samples = 400
    cov = np.array([[1., .7], [.7, 1.]]) * 10
    mean_1 = [.1, .1]
    mean_2 = [6., .1]

    np.random.seed(963)  # Included to enable reproducibility

    x_class1 = np.random.multivariate_normal(mean_1, cov, num_samples // 2)
    x_class2 = np.random.multivariate_normal(mean_2, cov, num_samples // 2)
    xy_class1 = np.column_stack((x_class1, np.zeros(num_samples // 2)))
    xy_class2 = np.column_stack((x_class2, np.ones(num_samples // 2)))
    data_full = np.row_stack([xy_class1, xy_class2])
    np.random.shuffle(data_full)
    data = data_full[:, :2]
    labels = data_full[:, 2]
    
    # Setup of initial Mu0
    Mu0 = np.zeros([data.shape[1], 2])
    Mu0[:, 1] = 1.
    
    KMeansAlgorithm(data, labels, Mu0, False)
    
    EMAlgorithm(data, labels, Mu0, False)

if __name__ == "__main__":
    main()


### Config 1
### High Accuracy - populations linearly separable
### Convergence achieved in 5/11 steps for K-Means/EM respectively
#    num_samples = 600
#    cov = np.array([[1.4, .7], [.7, 1.23]]) * 10

#    mean_1 = [-5, 4]
#    mean_2 = [12.6, 1.2]
    
### Config 2 - K-Means didnt converged (bad initialization for Mu)
#    num_samples = 400
#    cov = np.array([[1.4, .7], [.7, 1.23]]) * 5

#    mean_1 = [4., -1.8]
#    mean_2 = [12.6, 1.2]
    
    ### Config 3 - Both not converged (Max-Iter achieved)
    ### Bad accuracy due to non-convergence, axis were rotated
#   num_samples = 500
#    cov = np.array([[1.4, .7], [.7, 1.23]]) * 8

#    mean_1 = [4., -1.8]
#    mean_2 = [0.6, 2.2]

### Config 4
### High Accuracy(KM) - Convergence in 4 iterations, independent of initialization means
####       due to populations linearly separable
### Low Accuracy (not convergence) - EM - Bad selection of initial Means for Mu2
###        even with population linearly separable, EM may not converge due to poor selection
###        of initial values
#    num_samples = 600
#    cov = np.array([[1.4, .7], [.7, 1.23]]) * 10

#    mean_1 = [-5, 4]
#    mean_2 = [25.6, 11.2]
 
### Config 5
###     K-Means Fast convergence (11 iterations) good selection of initial values, population
###        Not linearly separable (highly mixed)
###     E.M. Not converged due to highly mixed samples - Centroids collapsed to values near
###        each other, with consequent poor separation.
    
#    num_samples = 900
#    cov = np.array([[1.4, -0.3], [-0.3, 1.23]]) * 20

#    mean_1 = [-5.5, 4.3]
#    mean_2 = [1.6, -2.5]
    
    
    
    
    