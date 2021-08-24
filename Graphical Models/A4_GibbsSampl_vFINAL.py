# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 00:02:17 2020
Updated on Fri Dec  4 00:07:53 2020
Updated on Sun Dec  6 0:04:58:53 2020

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import pickle
import winsound
import gc

matplotlib.rcParams.update({'font.size': 18})
matplotlib.rc('font',family='Arial')


''' My Original version - faster since it doesnt use numpy routines
'''

def AdjVert(p):
    Up = [(L-1) if p[0]==0 else p[0]-1, p[1]]
    Down = [0 if p[0]==(L-1) else p[0]+1, p[1]]
    Right = [p[0], 0 if p[1]==(L-1) else p[1]+1]
    Left = [p[0], (L-1) if p[1]==0 else p[1]-1]
    return(np.array([Up, Right, Down, Left]))

def Hamiltonian(x):
    S = 0.0
    for i in range(0,L):
        for j in range(0, L):
            A = AdjVert([i,j])
            V = [x[A[0,0],A[0,1]],x[A[1,0],A[1,1]],
                 x[A[2,0],A[2,1]],x[A[3,0],A[3,1]]]
            S = S - J*(x[i,j]*V[0]+x[i,j]*V[1]+
                       x[i,j]*V[2]+x[i,j]*V[3])
    return(S)

'''
OPTIMIZED - Returns the nearest neighbors of sigma(i,j). For nearest neighbor interactions,
these are sigma(i+1,j), sigma(i-1,j), sigma(i,j+1), and sigma(i,j-1).
'''
def get_neighbors(sigma,i,j):
    return np.array([sigma[(i+1)%L,j],
                     sigma[(i-1)%L,j],
                     sigma[i,(j+1)%L],
                     sigma[i,(j-1)%L]])    
    
'''
Calculates the value of Gamma for a given series of probability
'''
# def calcGamma(H_Gibbs, last):
#    Mag = np.sum(1/np.exp(-beta*H_Gibbs[:last]))
#    G = 1/(last*(2**N))*Mag   
#    return G

''' 
Returns p = Pr{sigma(i,j) = +1} and q = Pr{sigma(i,j) = -1}. These 
are otherwise known as the Boltzmann factors.
'''
def posterior(sigma_ij,neighbors,sigma):
#   Hp and Hm are the Hamiltonians for sigma_ij = +/- 1
    Hp = -1 * J * (neighbors[0]+neighbors[1]+neighbors[2]+neighbors[3])
    Hm = +1 * J * (neighbors[0]+neighbors[1]+neighbors[2]+neighbors[3])
#   p and q are the probability that sigma_ij is +/- 1
    p = np.exp(-beta*J*Hp)/(np.exp(-beta*J*Hp) + np.exp(-beta*J*Hm))
#    q = np.exp(-beta*J*Hm)/(np.exp(-beta*J*Hp) + np.exp(-beta*J*Hm))
    return p    

'''
Returns a zero-filled array sigma of dimension (i,j,K) which the Markov chain 
is stored in. That's not completely true. Before returning the array, the 
initial Markov state is defined in sigma[:,:,0], where each sigma(i,j,k=0) =
+/- 1 with equal probability.
'''
def initialize_sigma():
#   creates sigma(k=0), with elements +/-1 with equal probability
#   sigma_series stores sigma for every time-step
    sigma0 = np.random.choice([-1,1],size=L**2,p=[0.5,0.5]).reshape(L,L)
    sigma_series = np.zeros((L,L,K))
    H_Gibbs = np.zeros(K)
    sigma_series[:,:,0] = sigma0
    H_Gibbs[0] = Hamiltonian(sigma0)
    return sigma_series, H_Gibbs

'''
Each iteration, the Gibbs sampler selects one of the L^2 lattice elements
randomly, i.e. sigma(i,j). A new value of sigma(i,j) is then drawn from
the posterior distribution P[sigma(i,j) | all sigma(I!=i, J!=j)]. 
The posterior distribution includes only 4 sigma terms because the Ising 
model assumes nearest neighbor interactions: sigma(I=i+-1, J=j+-1). Note
that sigma being updated one (i,j) pair at a time is the characteristic
partial resampling feature of Gibbs sampling.
'''      
def Gibbs(SP):
#   sigma_series stores sigma for each step of the Markov chain
    sigma_series, H_Gibbs = initialize_sigma()
    Mag = np.zeros(K)
    Cum_Mag = 0
      # nitialize Hamiltonian 
#   selects (i,j) randomly and calls draw_sigma_ij(...) to update sigma(i,j)
    for k in range(1,K):
        i,j = np.random.randint(0,L,size=2)
        sigma_series[:,:,k] = Gibbs_update(sigma_series[:,:,k-1],i,j)
        H_Gibbs[k] = Hamiltonian(sigma_series[:,:,k])
        Mag[k] = 1/np.exp(-beta*J*H_Gibbs[k]) 
        Cum_Mag = Cum_Mag + Mag[k]
        G = Cum_Mag/(k*(2**N))
        SP[k] = np.log(1/G)/N
    return sigma_series, SP    

    
'''
Returns the new state of sigma after updating sigma(i,j) according to the
posterior distribution used in Gibbs sampling.
'''
def Gibbs_update(sigma,i,j):
#   The nearest neighbors of sigma(i,j) = sigma(i+1,j), sigma(i-1,j), 
#   sigma(i,j+1), and sigma(i,j-1).
    neighbors = get_neighbors(sigma,i,j)
#   p and q are the probabilities of being spin +1 or spin -1, respectively.
#   p + q = 1.
    p = posterior(sigma[i,j],neighbors,sigma)
#   sets sigma(i,j) according to wp and wm
    if np.random.rand() < p:
        sigma[i,j] = 1
    else:
        sigma[i,j] = -1
#   returns new state for sigma
    return sigma

def simulateFreeEnergy (SP_FGibbs):
    print( "starting Simulations")
    
    Gibbs_start_time = time.time()
    
#   simulates sigma using Gibbs sampler
    print( "Starting Gibbs sampler for Beta=%f -- %d simulations with %d samples each\n" % (beta, NoSim, K))

    sigma_Gibbs, SP_FGibbs = Gibbs(SP_FGibbs)           #   Generates K Samples using Gibbs Sampler

    ExecTime = (time.time() - Gibbs_start_time)
    print("\nSimulations (Exec. time): --- %s (in seconds) ---\n" % ExecTime)
    return SP_FGibbs

'''
Plots sample path for each 
'''
def plot_sample_paths(SP, lb_beta, lb_plot, trueF):
    plt.figure(figsize=(12,10))  
    plt.xscale("log")
    plt.plot(SP[0,:]); 
    plt.plot(SP[1,:]); 
    plt.plot(SP[2,:]); 
    plt.plot(SP[3,:]);
#    plt.axhline(y=trueF, linewidth=1, linestyle = 'dashed', color='black')
    plt.xlabel('log(Samples)')
    plt.savefig(lb_plot)
    plt.close()


''' 
####   Main Module
'''
def saveSamplePath(SP, ID, name, job):
    
    print("Saving Sample Paths...")
    
    with open(name+'_ID'+str(ID)+'_'+str(job)+'.dat','wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(SP[0,:], f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()

    print("\nConcluded!\n")
    
def getBetaKFromUser():
    
    value1 = input("Please enter the value of beta: ")
    value1 = float(value1)
    YN1 = input(f'You entered {value1}. Is this correct[Y/N]: ')
    
    if (YN1 == 'Y' or YN1 == 'y'):
       value2 = input("Please enter the number of samples: ")
       YN2 = input(f'You entered {value2}. Is this correct[Y/N]: ')
       value2 = int(value2)
       if (YN2 == 'Y' or YN2 == 'y'):
           return value1, value2
       else:
           return -1, -1
    else:
        return -1, 0

def getNFromUser():
    
    value1 = input("Please enter the size of the lattice [1-15]: ")
    value1 = int(value1)
    YN1 = input(f'You entered {value1}. Is this correct[Y/N]: ')
    
    if (YN1 == 'Y' or YN1 == 'y'):
       if not(value1 in range(15)):
           return -1
       else:
           return value1
    else:
        return -1

    
def Generate2DIsing():
    
    # Setup Global variables
    global N, L, K, beta, lBeta, J, NoSim, burnin, job, start_time, Version
    global SP_FGibbs_TMP

    
    Version = '_V43_'
    
    # L = 5      # Dimension of 2D Ising Model
    # K = 3*10**4  # Number of Samples per Gibbs simulation
    J = 0.01
    
    NoSim = 4  # Number of Monte Carlo Simulations per beta (= no. sample paths)
    
    lBeta = [0.01, 0.06, 0.15, 0.25, 0.75]
 
    beta, K = getBetaKFromUser()  # Get input for Beta and K from user
    L = getNFromUser()

    if (beta == -1):
        return False
    else:
        if not(beta in lBeta):
            print('\nInvalid Beta (%f), program terminated.' % beta)
            return False
        else:
            if (K <= 0):
                print('\nInvalid No. of Samples (%d), program terminated.' % K)
                return False
            else:
                if(L < 0):
                    print('\nInvalid size of Lattice (%d), program terminated.' % L)
                    return False
   
    N = L**2   # for compatibility 
    burnin = int(K*0.01)
    
    start_time = time.time()
    
    durationP = 800  # milliseconds
    freqP = 1000  # Hz
    durationJob = 1000  # milliseconds
    freqJob = 1500  # Hz
    
    job = int(round((start_time-int(start_time))*10**5,0))
    
    print("\n--->>> JOB No. %d ------------" % job)
        
    np.random.seed(123)   # Set seed to reproduce/test results

    SP_FGibbs_TMP = np.zeros((1, K))
    
    # Runs 04 Sample-Paths per Beta
    print('\nBLOCK - Beta = %f ---' % beta)
    for i in range(NoSim):  # plot 4 Sample paths per beta
        print('\n---- Sample Path %d - Beta = %f ---' % (i+1, beta))
        SP_FGibbs_TMP[0,:] = simulateFreeEnergy(SP_FGibbs_TMP[0,:])
        saveSamplePath(SP_FGibbs_TMP, i, 'SP_beta'+str(beta), job)
        winsound.Beep(freqP, durationP)
        gc.collect()

    ExecTime = (time.time() - start_time)
    print("\nGibbs Sampling Execution time: --- %s (in seconds) ---\n" % ExecTime)
    
    winsound.Beep(freqJob, durationJob)
    
    return True
        
def Plot2DIsingSamples():
    
    TrueZ = [3.3638430074774e07, 3.6724646193133e07, 5.9592513591093e07, 
             1.78755383377176e08, 4.1401944265394672e16]
    
    # Reads the Sample Paths to generate the Graphs
    SP_FGibbs_Graph = np.zeros((NoSim, K))
    
    for i in range(NoSim): 
        with open('SP_beta'+str(beta)+'_ID'+str(i)+'_'+str(job)+'.dat','rb') as f:
            SP_FGibbs_Graph[i,:] = pickle.load(f)
        f.close()
    
    plot_sample_paths(SP_FGibbs_Graph[:,burnin:], 'Beta='+str(beta), 'SP_beta'+str(beta)+Version+
                      str(job)+'.png',np.log(TrueZ[lBeta.index(beta)])/N)
   
    ExecTime = (time.time() - start_time)
    print("\n>>> JOB No. %d - Execution time: --- %s (in seconds) ---\n" % (job, ExecTime))
    

if Generate2DIsing():
    Plot2DIsingSamples()
    del SP_FGibbs_TMP
    gc.collect()


#SP_FGibbs_Graph = np.zeros((NoSim, 3*10**7))

#job = 59351

#for i in range(NoSim): 
#    with open('SP_beta0.25'+'_ID'+str(i)+'_'+str(job)+'.dat','rb') as f:
#        SP_FGibbs_Graph[i,:] = pickle.load(f)
#    f.close()

#plot_sample_paths(SP_FGibbs_Graph[:,burnin:], 'Beta=0.25', 'SP_beta025_V37_'+str(job)+'.png',
#                  np.log(1.78755383377176e08)/N)

