# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 00:30:31 2021
Updated on Sat Jan 30 17:26:02 2021
Updated on Mon Feb 01 10:48:04 2021
Updated on Fri Feb 05 14:43:10 2021 - Include list of weights (as of Piazza @53_1)

@author: Luis Alvaro Correia - Std No.1006508566
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# Sets variable for Plotting routines
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rc('font',family='Arial')

# Function which takes 'data' as argument and returns its randomly permuted version
# along the samples. Here we consider uniformly random permutation of training data.
def shuffle_data(data):
    sh_idx = np.random.choice(len(data['t']),size=len(data['t']), replace=False)
    data_shf = {'X':data['X'][sh_idx,:], 't':data['t'][sh_idx]}
    return(data_shf)
    
# Function which takes 'data', number of partitions as 'num-folds' and the selected
# partition fold as its arguments and returns the selected partition block 'fold' as
# 'data_fold' and the remaining data as 'data_rest' in 02 disjoint sets.
def split_data(data, num_folds, fold):
   
    # Calculates the size of each fold
    sz_fold = int(len(data['t'])/num_folds) 
    
    # Calculates the limits of validation fold
    idx_fold = list(range((fold-1)*sz_fold, fold*sz_fold-1))
    
    # Calculates the difference of the two sets
    idx_rest = list(set(range(0,len(data['t'])))-set(idx_fold))
    
    # Assign the correct portions of each folder set
    data_fold = {'X':data['X'][idx_fold,:], 't':data['t'][idx_fold]}
    data_rest = {'X':data['X'][idx_rest,:], 't':data['t'][idx_rest]}
    
    return(data_fold, data_rest)
    
# Function which takes 'data' and 'lambd' and returns the coefficients of Ridge Regression
# with penalty level 'lambd' as given in equation (2.1) in HW1.
def train_model(data, lambd):
    phi_prod = np.matmul(np.transpose(data['X']),data['X'])
    inv_phi_lmdb = np.linalg.inv(phi_prod+lambd*np.identity(len(data['X'][0,:])))
    model = np.matmul(np.matmul(inv_phi_lmdb,np.transpose(data['X'])),data['t'])
    return(model)
    
# Function which takes 'data' and 'model' and returns the predictions based in 'data' and
# 'model'.
def predict(data, model):
    predictions = np.matmul(data['X'],model)
    return(predictions)
    
# Function which takes 'data' and 'model' and returns the average squared error loss based
# in 'model'. This means that if 'data' is composed of t and Phi and hat(w) then the return
# will be the norm of ||t - Phi.dat(w)||^2/n
def loss(data, model):
    error = pow(np.linalg.norm(data['t']-np.matmul(data['X'],model)),2)/len(data['t'])
    return(error)
    
# Function that takes the training 'data', number of folds 'num_folds', and a sequence of
# lambdas as 'lambd_seq' as its arguments and return the cross validation error across all
# lambdas. Take 'lambda_seq' as evenly spaced 50 numbers over the the interval (0.02,1.5).
# This means 'cv_error' will be a vector of 50 errors corresponding to the values of 
# 'lambda_seq'.
def cross_validation(data, num_folds, lambd_seq):
    cv_error = np.zeros(len(lambd_seq))
    model = np.zeros((len(lambd_seq),len(data['X'][0,:])))
        
    data = shuffle_data(data)
    for i in range(1,len(lambd_seq)):
        lambd = lambd_seq[i]
        cv_loss_lmd = 0.0
        for fold in range(1,num_folds):
            val_cv, train_cv = split_data(data, num_folds, fold)
            model[i,:] = train_model(train_cv, lambd)
            cv_loss_lmd += loss(val_cv, model[i,:])
        cv_error[i] = cv_loss_lmd/num_folds
    return(cv_error, model)
    
# Function that takes the training and test data and return the cross validation 
# error across all lambdas for the adjusted model.
def training_test_validation(data_tr, data_tst, lambd_seq):
    cv_error_tr = np.zeros(len(lambd_seq))
    cv_error_tst = np.zeros(len(lambd_seq))
    model = np.zeros((len(lambd_seq),len(data_tr['X'][0,:])))
        
    for i in range(1,len(lambd_seq)):
        lambd = lambd_seq[i]
        model[i,:] = train_model(data_tr, lambd)
        cv_error_tr[i] = loss(data_tr, model[i,:])
        cv_error_tst[i] = loss(data_tst, model[i,:])
    
    return(cv_error_tr, cv_error_tst, model)

'''
Auxiliary function to generate sequences of lambdas for Tidge Regression 
''' 
def generate_seq_lambd(start, end, segments):
    lmbd_lst = []
    if (end < start):
        return(-1)
    lmbd = start
    increment = (end-start)/(segments+1)
    while (lmbd<end):
        lmbd_lst.append(lmbd)
        lmbd += increment
    return lmbd_lst

'''
Plots each error curve obtained after training & testing 
'''
def plot_Graph(err_train, err_tst, cv1, cv2, lambd_seq, name_plt, lb_plot):
    plt.figure(figsize=(12,10))  
    plt.title(name_plt)
    plt.plot(lambd_seq[1:], err_train[1:], label = 'Training'); 
    plt.plot(lambd_seq[1:], err_tst[1:], label = 'Test'); 
    plt.plot(lambd_seq[1:], cv1[1:], label = '5-fold')
    plt.plot(lambd_seq[1:], cv2[1:], label = '10-fold')
    plt.xlabel('Lambda')
    plt.ylabel('Error')
    plt.legend()
    plt.show()
    #plt.savefig(lb_plot)
    #plt.close()
    
'''
List Solution obtained for each procedure, for a set of lambdas
'''
def list_Model(model, lambd_seq, name_lst):
    print('\n----------- Listing of Solutions per Lambda (%s) -----------\n' % name_lst)
    for i in range(1,len(lambd_seq)):
        print('\n>>> Lambda - %2.6f\n' % lambd_seq[i])
        print(model[i,:],'\n')

# Main program - Initialization procedures, load data, runs Cross Validation
data_train = {'X':np.genfromtxt('data\\data_train_X.csv', delimiter =','),
              't':np.genfromtxt('data\\data_train_y.csv', delimiter = ',')}
data_test = {'X':np.genfromtxt('data\data_test_X.csv', delimiter =','),
             't':np.genfromtxt('data\data_test_y.csv', delimiter = ',')}

np.random.seed(123)   # Set seed to reproduce/test results

lambd_seq = generate_seq_lambd(0.02, 1.5, 50)

# Training and Test Data validation
error_tr, error_tst, model_trained = training_test_validation(data_train, data_test, lambd_seq)

# 5-Fold - Cross Validation
cross_val5, model_val5 = cross_validation(data_train, 5, lambd_seq)

# 10-Fold - Cross Validation
cross_val10, model_val10 = cross_validation(data_train, 10, lambd_seq)

# Part(b) - Print the Errors obtained
d = {'Lambda':lambd_seq,
     'Training':error_tr,
     'CV-5Fold': cross_val5,
     'CV-10Fold': cross_val10,
     'Test':error_tst
     }
df = pd.DataFrame(d)  

# Print errors, Lambdas for each process
print('\n----------- Listing of Errors per procedure -----------\n')
print(df[1:])

# Part (c) - Plot Graph and prints outputs
plot_Graph(error_tr, error_tst, cross_val5, cross_val10, 
           lambd_seq, 'Training, Test and 5/10-Fold CV', 'TTCV_5_10Fold_v4')

print('\n----------- Processing Summary --------------')


print('\nThe value of lambda proposed by Cross-Validation Procedure is as follows:')
print('\n>>> Minimum Error Lambda for CV-5Fold=%2.6f' % lambd_seq[df[1:].idxmin()['CV-5Fold']])
print('\n>>> Minimum Error Lambda for CV-10Fold=%2.6f' % lambd_seq[df[1:].idxmin()['CV-10Fold']])

print('\n\nThe value of Lambda with best performance after regular training:')
print('\n>>> Minimum Error Lambda after Training=%2.6f' % lambd_seq[df[1:].idxmin()['Test']])

print('\n\nThe errors of best performance K-Fold models over Test Data are:')
print('\n>>> Test Data with Optimum Model 5-Fold=%2.6f' % loss(data_test, model_val5[df[1:].idxmin()['CV-5Fold']]))
print('\n>>> Test Data with Optimum Model 10-Fold=%2.6f' % loss(data_test, model_val10[df[1:].idxmin()['CV-10Fold']]))

print('\n\n------ END PROCESSING ------\n')

print('\n\n------ APPENDIX - List of Solutions per Procedure ------\n')

list_Model(model_trained, lambd_seq, 'TRAINING')