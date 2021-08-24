# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 00:02:17 2020
Updated on Fri Dec  4 00:07:53 2020
Updated on Sun Dec  6 0:04:58:53 2020

Author - Luis Alvaro Correia

Purpose - Print Sample Paths from Gibbs Sampler from Ising Model 2D

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import time
# import pickle
import gc

matplotlib.rcParams.update({'font.size': 18})
matplotlib.rc('font',family='Arial')


'''
Plots sample path for each 
'''
def plot_sample_paths(SP, lb_beta, lb_plot):
    plt.figure(figsize=(12,10))  
    plt.xscale("log")
    plt.plot(SP[0,:], label=lb_beta); 
    plt.plot(SP[1,:]); 
    plt.plot(SP[2,:]); 
    plt.plot(SP[3,:]);
    # plt.axhline(y=trueF, linewidth=1, linestyle = 'dashed', color='black')
    plt.xlabel('log(Samples)')
    plt.legend()
    plt.savefig(lb_plot)
    plt.close()
    
def plot_Error(SP, lb_beta, lb_plot):
    plt.figure(figsize=(12,10))
    plt.xscale("log")
    plt.plot(SP[0,:], label=lb_beta); 
    plt.plot(SP[1,:]); 
    plt.plot(SP[2,:]); 
    plt.plot(SP[3,:]);
    plt.xlabel('log(Samples)')
#    plt.title(tit+' : N='+str(N)+' / K='+str(K), fontsize=18)
    plt.legend()
    plt.savefig(lb_plot)
    plt.close()

def plot_hist_paths(SP, lb_beta, lb_plot):
    # Draw Join Histogram of Sample Paths
    plt.figure(figsize=(13,10), dpi= 80)
    sns.distplot(SP[0,:], color="dodgerblue", label="SP 01", hist_kws={'alpha':.7}, kde_kws={'linewidth':3})
    sns.distplot(SP[1,:], color="orange", label="SP 02", hist_kws={'alpha':.7}, kde_kws={'linewidth':3})
    sns.distplot(SP[2,:], color="g", label="SP 03", hist_kws={'alpha':.7}, kde_kws={'linewidth':3})
    sns.distplot(SP[3,:], color="gray", label="SP 04", hist_kws={'alpha':.7}, kde_kws={'linewidth':3})
    # Decoration
    plt.title('Free Energy (Gibbs Sampling) K='+str(K)+' / '+lb_beta, fontsize=18)
    plt.legend()
#    plt.show()
    plt.savefig(lb_plot)
    plt.close()

    
''' 
####   Main Module
'''
   
def getBetaJobKFromUser():
    
    value1 = input("Please enter the value of beta: ")
    value1 = float(value1)
    YN = input(f'You entered {value1}. Is this correct[Y/N]: ')
    
    if (YN != 'Y' and YN != 'y'):
        return -1, -1, -1
    
    value2 = input("Please enter the number of job: ")
    value2 = int(value2)
    YN = input(f'You entered {value2}. Is this correct[Y/N]: ')
    if (YN != 'Y' and YN != 'y'):
        return value1, -1, -1
    
    value3 = input("Please enter the number of samples: ")
    value3 = int(value3)
    YN = input(f'You entered {value3}. Is this correct[Y/N]: ')
    if (YN != 'Y' and YN != 'y'):
        return value1, value2, -1
    
    return value1, value2, value3
    
def Print2DIsing():
    
    # Setup Global variables
    global L, N, K, beta, lBeta, NoSim, burnin, job, start_time, Version
    
    Version = '_V43_'
    
    L = 10      # Dimension of 2D Ising Model
    N = L**2   # for compatibility 
    K = 10**4  # Number of Samples per Gibbs simulation
    
    NoSim = 4  # Number of Monte Carlo Simulations per beta (= no. sample paths)
    
    lBeta = [0.01, 0.06, 0.15, 0.25, 0.75]
 
    beta, job, K = getBetaJobKFromUser()  # Get input for Beta and K from user

    if (beta == -1):
        return False
    #else:
    #    if not(beta in lBeta):
    #        print('\nInvalid Beta (%f), program terminated.' % beta)
    #        return False

    if (job == -1):
        return False

    if (K == -1):
        return False

    burnin = int(K*0.01)
    
    start_time = time.time()
    
    print("\n--->>> JOB No. %d ------------" % job)
        
    return True
        
def Plot2DIsingSamples():
    
    global SP_FGibbs_Graph
    
    #TrueZ = [3.3638430074774e07, 3.6724646193133e07, 5.9592513591093e07, 
    #         1.78755383377176e08, 4.1401944265394672e16]
    
    # Reads the Sample Paths to generate the Graphs
    SP_FGibbs_Graph = np.zeros((NoSim, K))
    
    for i in range(NoSim):
        print("\nReading Sample Path No. %d..." % (i+1))
        SP_FGibbs_Graph[i,:] = np.fromfile('SP_beta'+str(beta)+'0000_v6-CPP__ID'+str(i)+'_'+str(job)+'.dat', dtype=float)
#        with open('SP_beta'+str(beta)+'_ID'+str(i)+'_'+str(job)+'.dat','rb') as f:
#            SP_FGibbs_Graph[i,:] = pickle.load(f)
#        f.close()
    
    print("\nGenerating graph of Sample Paths...")
    
    plot_sample_paths(SP_FGibbs_Graph[:,burnin:], 'Beta='+str(beta),
                      'SP_beta'+str(beta)+Version+str(job)+'_PRT.png')
   
    #TZ = TrueZ[lBeta.index(beta)]

    #print("\nGenerating graph of Absolute Error...")
    #plot_Error(abs(SP_FGibbs_Graph[:,burnin:]-np.log(TZ))/N, 'Beta='+str(beta),
    #            'AE_beta'+str(beta)+Version+str(job)+'_PRT.png')

    #print("\nGenerating graph of Relative Error...")
    #plot_Error(abs(SP_FGibbs_Graph[:,burnin:]-np.log(TZ))/np.log(TZ), 'Beta='+str(beta),
    #           'RE_beta'+str(beta)+Version+str(job)+'_PRT.png')

    #print("\nGenerating Histograms of Sample Paths...")
    #plot_hist_paths(SP_FGibbs_Graph[:,burnin:], 'Beta='+str(beta),
    #                'HST_beta'+str(beta)+Version+str(job)+'_PRT.png')

    ExecTime = (time.time() - start_time)
    print("\n>>> Print JOB No. %d - Execution time: --- %s (in seconds) ---\n" % (job, ExecTime))
    
if Print2DIsing():
    Plot2DIsingSamples()
    del SP_FGibbs_Graph
    gc.collect()


