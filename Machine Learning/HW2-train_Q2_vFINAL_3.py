# -*- coding: utf-8 -*-
"""
@author: Luis Alvaro Correia - Std No.1006508566
"""

from __future__ import absolute_import
from __future__ import print_function
from future.standard_library import install_aliases
install_aliases()
import numpy as np
import pandas as pd
from scipy.special import logsumexp
import os
import gzip
import struct
import array
import matplotlib.pyplot as plt
import matplotlib.image
from urllib.request import urlretrieve

def download(url, filename):
    if not os.path.exists('data'):
        os.makedirs('data')
    out_file = os.path.join('data', filename)
    if not os.path.isfile(out_file):
        urlretrieve(url, out_file)


def mnist():
    base_url = 'http://yann.lecun.com/exdb/mnist/'

    def parse_labels(filename):
        with gzip.open(filename, 'rb') as fh:
            magic, num_data = struct.unpack(">II", fh.read(8))
            return np.array(array.array("B", fh.read()), dtype=np.uint8)

    def parse_images(filename):
        with gzip.open(filename, 'rb') as fh:
            magic, num_data, rows, cols = struct.unpack(">IIII", fh.read(16))
            return np.array(array.array("B", fh.read()), dtype=np.uint8).reshape(num_data, rows, cols)

    for filename in ['train-images-idx3-ubyte.gz',
                     'train-labels-idx1-ubyte.gz',
                     't10k-images-idx3-ubyte.gz',
                     't10k-labels-idx1-ubyte.gz']:
        download(base_url + filename, filename)

    train_images = parse_images('data/train-images-idx3-ubyte.gz')
    train_labels = parse_labels('data/train-labels-idx1-ubyte.gz')
    test_images = parse_images('data/t10k-images-idx3-ubyte.gz')
    test_labels = parse_labels('data/t10k-labels-idx1-ubyte.gz')

    return train_images, train_labels, test_images[:10000], test_labels[:10000]


def load_mnist(N_data=None):
    partial_flatten = lambda x: np.reshape(x, (x.shape[0], np.prod(x.shape[1:])))
    one_hot = lambda x, k: np.array(x[:, None] == np.arange(k)[None, :], dtype=int)
    train_images, train_labels, test_images, test_labels = mnist()
    train_images = (partial_flatten(train_images) / 255.0 > .5).astype(float)
    test_images = (partial_flatten(test_images) / 255.0 > .5).astype(float)
    K_data = 10
    train_labels = one_hot(train_labels, K_data)
    test_labels = one_hot(test_labels, K_data)
    if N_data is not None:
        train_images = train_images[:N_data, :]
        train_labels = train_labels[:N_data, :]

    return train_images, train_labels, test_images, test_labels


def plot_images(images, ax, ims_per_row=5, padding=5, digit_dimensions=(28, 28),
                cmap=matplotlib.cm.binary, vmin=None, vmax=None):
    """Images should be a (N_images x pixels) matrix."""
    N_images = images.shape[0]
    N_rows = np.int32(np.ceil(float(N_images) / ims_per_row))
    pad_value = np.min(images.ravel())
    concat_images = np.full(((digit_dimensions[0] + padding) * N_rows + padding,
                             (digit_dimensions[1] + padding) * ims_per_row + padding), pad_value)
    for i in range(N_images):
        cur_image = np.reshape(images[i, :], digit_dimensions)
        row_ix = i // ims_per_row
        col_ix = i % ims_per_row
        row_start = padding + (padding + digit_dimensions[0]) * row_ix
        col_start = padding + (padding + digit_dimensions[1]) * col_ix
        concat_images[row_start: row_start + digit_dimensions[0],
                      col_start: col_start + digit_dimensions[1]] = cur_image
        cax = ax.matshow(concat_images, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.xticks(np.array([]))
        plt.yticks(np.array([]))
    return cax


def save_images(images, filename, **kwargs):
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)
    plot_images(images, ax, **kwargs)
    fig.patch.set_visible(False)
    ax.patch.set_visible(False)
    plt.savefig(filename)
    plt.close()  # Included
    
def train_log_regression(images, labels, learning_rate, max_iter):
    """ Used in Q1
        Inputs: train_images, train_labels, learning rate,
        and max num of iterations in gradient descent
        Returns the trained weights (w/o intercept)"""
    N_data, D_data = images.shape
    K_data = labels.shape[1]
    weights = np.zeros((D_data, K_data))
    
    # YOU NEED TO WRITE THIS PART
    
    ### - Start Coding 01 - ###
    cross_ent_history = np.zeros(max_iter)
    accuracy_history = np.zeros(max_iter)
    grad_history = np.zeros(max_iter)
    it = 0
    converged = False
    
    # v.FINAL2 : Here, differently in the previous version i defined the
    # tolerance to check if accuracy stabilizes from previous step. 
    tolerance = 1e-5
    
    while (not converged and (it < max_iter)):
        # Calculate the log-Softmax for the calculated weights
        lsmax = log_softmax(images, weights)
        
        # Calculate the prediction
        psmax = predict(lsmax)
        
        # Calculate the Cross Entropy over the prediction
        cent = cross_ent(psmax, labels)

        # Calculate the Gradient
        grad = (images.T).dot(np.exp(lsmax)-labels)

        # Update Weights with the current gradient descent value
        weights -= grad*learning_rate
        
        # Stores Gradient, Cross-Entropy and Accuracy History
        grad_history[it] = np.linalg.norm(grad)
        cross_ent_history[it] = cent  
        accuracy_history[it] = accuracy(lsmax, labels)
        
        # Check for convergence under a tolerance defined
        # converged = (grad_history[it]<=tolerance)  - OLD Convergence Criteria
        
        # v-FINAL2: Initially I was checking the size of gradient 
        # to check convergence, but I noticed it takes lots of interations
        # to became small, so I changed to check the convergence of 
        if (it > 0):
            converged = (abs(accuracy_history[it]-accuracy_history[it-1])<=tolerance)
        
        # Increment iteration
        it += 1

    # Save Cross-Entrpy for later reference
    plt.figure(figsize=(10,7))  
    plt.title('Cross entropy History Learning rate='+
              str(learning_rate)+' / Max-Iter='+str(max_iter), fontsize=12)
    plt.plot(cross_ent_history); 
    plt.xlabel('Iteration')
    plt.ylabel('Cross-Entropy')
    # plt.show()
    plt.savefig('Cross-Entropy-GD.png')
    plt.close()    

    # Save Training Accuracy for later reference
    plt.figure(figsize=(10,7))  
    plt.title('Training Accuracy History | Learning rate='+
              str(learning_rate)+' / Max-Iter='+str(max_iter), fontsize=12)
    plt.plot(accuracy_history[:it]); 
    plt.xlabel('Iteration')
    plt.ylabel('Accuracy')
    # plt.show()
    plt.savefig('Trn_Accuracy-GD.png')
    plt.close()    
        
    #### - End Coding 01 - ###
    
    w0 = None # No intercept for log-reg
    return weights, w0


def train_gda(images, labels):
    """ Used in Q2
        Inputs: train_images, train_labels
        Returns the trained weights, the intercept, and D x K class means, 
        and D x D common covariance matrix."""
    N_data, D_data = images.shape
    K_data = labels.shape[1]
    
    # YOU NEED TO WRITE THIS PART
    ### - Start Coding 02 - ###
    
    # Calculating the number of items for each class
    N_k = np.sum(labels, axis=0)
    
    # Estimator for the Priors of each class
    Pi_k = N_k/N_data
    
    # Estimator for Mu_k
    Mu = ((labels.T).dot(images)).T[:,:K_data] / N_k[:K_data]
    
    # Initialize Sigma with non-zeros to avoid singulatity
    Sigma_hat = np.identity(D_data)*1e-6
    
    # Estimate Sigma_k and Sigma_hat

    for k in range(K_data):
        # Select Images from Class-k
        Idx_k = np.where(labels[:,k])[0]
        
        # Estimator for Sigma_k
        Sigma_hat_k = (1/N_k[k])*(images[Idx_k,:]-Mu[:,k]).T.dot((images[Idx_k,:]-Mu[:,k]))
        
        Sigma_hat += (N_k[k]/N_data)*Sigma_hat_k
    
    Inv_Sigma_hat = np.linalg.inv(Sigma_hat)
    
    # Initializing Intercept wk0 and weights for all Classes
    w0 = np.zeros(K_data)
    weights = np.zeros((D_data, K_data))
    for k in range(K_data):
        # Calculating Intercept wk0 for each Classe
        w0[k] = -0.5*((Mu[:,k].T).dot(Inv_Sigma_hat)).dot(Mu[:,k])+np.log(Pi_k[k])
        
        # Calculating the weights for each class
        weights[:,k] = Inv_Sigma_hat.dot(Mu[:,k])
            
    #### - End Coding 02 - ###
    
    return weights, w0, Mu, Sigma_hat

def train_gda2(images, labels):
    """ Used in Q2
        Inputs: train_images, train_labels
        Returns the trained weights, the intercept, and D x K class means, 
        and D x D common covariance matrix."""
    N_data, D_data = images.shape
    K_data = labels.shape[1]
    
    n_k = np.sum(labels, axis = 0)
    Mu = images.T @ (labels / n_k)
    delta = (images.T - Mu @ (labels).T)/np.sqrt(N_data)
    Sigma = delta @ delta.T + 1/N_data * np.eye(D_data)
    weights = np.linalg.inv(Sigma) @ Mu
    w0 = np.diag(-1/2*Mu.T @ np.linalg.inv(Sigma) @ Mu) + np.log(n_k/N_data)
    
    return weights, w0, Mu, Sigma


def log_softmax(images, weights, w0=None):
    """ Used in Q1 and Q2
        Inputs: images, and weights
        Returns the log_softmax values."""
    if w0 is None: w0 = np.zeros(weights.shape[1])

    # YOU NEED TO WRITE THIS PART

    ### - Start Coding - ###
    XW = images.dot(weights) + w0
    return XW - np.array([logsumexp(XW,axis=1),]).T
    ### - End Coding - ###


def cross_ent(log_Y, train_labels):
    """ Used in Q1
        Inputs: log of softmax values and training labels
        Returns the cross entropy."""

    # YOU NEED TO WRITE THIS PART
    
    ### - Start Coding - ###
    return(- np.sum(train_labels * (log_Y)))
    ### - End Coding - ###

def predict(log_softmax):
    """ Used in Q1 and Q2
        Inputs: matrix of log softmax values
        Returns the predictions"""

    # YOU NEED TO WRITE THIS PART
    
    ### - Start Coding - ###
    values = log_softmax.argmax(axis=1)
    return(np.eye(log_softmax.shape[1])[values])
    ### - End Coding - ###
    
def accuracy(log_softmax, labels):
    """ Used in Q1 and Q2
        Inputs: matrix of log softmax values and 1-of-K labels
        Returns the accuracy based on predictions from log likelihood values"""
    
    # YOU NEED TO WRITE THIS PART
    ### - Start Coding - ###
    notmtc = np.count_nonzero(labels.argmax(axis=1)-log_softmax.argmax(axis=1))
    return(1-notmtc/labels.shape[0])
    ### - End Coding - ###

def main():
    N_data = 60000 # Num of samples to be used in training
    # Set this to a small number while experimenting.
    # For log reg, finally use the entire training dataset for training (N_data=None).
    # For gda, use as many training samples as your computer can handle.
    
    print('\nLoading MNIST...\n')

    train_images, train_labels, test_images, test_labels = load_mnist(N_data)

    # Q1: train logistic regression
    #learning_rate, max_iter = .00001, 100
    #Proc = 'Multi-Class Logistic Regression'
    #ProcID = 'MC-LogReg'

    # Q1(c) - Training Multiclass Logistic Regression
    # print('\nProcessing %s...\n' % Proc)
    #weights, w0 = train_log_regression(train_images, train_labels, learning_rate, max_iter)
    
    # Q1(d) - Saving the Weights as 10 images
    #print('\nSaving Images...\n')
    #save_images(weights.T, 'weights.png')

    # Q2(b) - Train gaussian discriminant
    Proc = 'Gaussian Discriminant Analysis'
    ProcID = 'GDA'
    print('Processing %s...\n' % Proc)
    # weights, w0, Mu, Sigma = train_gda(train_images, train_labels)
    weights, w0, Mu, Sigma = train_gda2(train_images, train_labels)
    save_images(Mu.T, 'means.png')

    # Q2(e) - Using the Generative model generate 10 samples of digit '0' and '3'
    
    # Generate images of No. 0
    np.random.seed(123)
    new_digit = 0
    print('\nGenerating & Saving Images...\n')
    new_images = np.random.multivariate_normal(Mu[:, new_digit], Sigma, 10)
    save_images((new_images > .5).astype(float), 'new_images_0.png')

    # Generate images of No. 3
    new_digit = 3
    new_images = np.random.multivariate_normal(Mu[:, new_digit], Sigma, 10)
    save_images((new_images > .5).astype(float), 'new_images_3.png')

    # Q1-Q2(c)  - Generate Predictions no Training and Test dats-sets
    #           - Report Errors and Accuracy for both processes
    
    log_softmax_train = log_softmax(train_images, weights, w0)
    log_softmax_test = log_softmax(test_images, weights, w0)
    
    train_accuracy = accuracy(log_softmax_train, train_labels)*100.0
    test_accuracy = accuracy(log_softmax_test[:10000], test_labels[:10000])*100.0
    
    print('\nProcessing Report...\n')

    f=open('Summary_'+ProcID+'.prn','w')
    f.write("\n------ Processing Summary (%s) -----------\n" % Proc)
    f.write("\nTraining accuracy is %5.2f%%" % train_accuracy)
    f.write("\nTest accuracy is %5.2f%%\n" % test_accuracy)

    f.write("\n------ Listing of Errors from TRAINING procedure (%d Samples) ------\n" % N_data)

    ErrTrain = np.nonzero(train_labels.argmax(axis=1)-log_softmax_train.argmax(axis=1))[0]
    f.write("\n>>> Total No. of Errors: %6d\n" % ErrTrain.shape[0])
    f.write("\n>>> Note: Only first 120 will be printed due to limitation on Crowdmark\n")
    f.write("          to manage high no. of pages in uploaded PDF format.\n\n")
    dtrain = {'Img No.':ErrTrain,
              'Label':train_labels[ErrTrain].argmax(axis=1),
              'Classif.':log_softmax_train[ErrTrain].argmax(axis=1)}
    dftrain = pd.DataFrame(dtrain)
    f.write(dftrain[:120].to_string(header=True, index=True))
    #f.write(dftrain.to_string(header=True, index=True))

    f.write("\n\n------ Listing of Errors from TESTING procedure (%d Samples)------\n" % 10000)
    ErrTest = np.nonzero(test_labels.argmax(axis=1)-log_softmax_test[:10000].argmax(axis=1))[0]
    f.write("\n>>> Total No. of Errors: %6d\n" % ErrTest.shape[0])
    f.write("\n>>> Note: Only first 120 will be printed due to limitation on Crowdmark\n")
    f.write("          to manage high no. of pages in uploaded PDF format.\n\n")
    dtest = {'Img No.':ErrTest,
             'Label':test_labels[ErrTest].argmax(axis=1),
             'Classif.':log_softmax_test[ErrTest].argmax(axis=1)}
    dftest = pd.DataFrame(dtest)
    f.write(dftest[:120].to_string(header=True, index=True))
    #f.write(dftest.to_string(header=True, index=True))

    f.close()
    
    print('\nProcessing concluded!')
    
if __name__ == '__main__':
    main()
