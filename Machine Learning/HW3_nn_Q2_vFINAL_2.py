""" Instruction:

In this section, you are asked to train a NN with different hyperparameters.
To start with training, you need to fill in the incomplete code. There are 3
places that you need to complete:
a) Backward pass equations for an affine layer (linear transformation + bias).
b) Backward pass equations for ReLU activation function.
c) Weight update equations.

After correctly fill in the code, modify the hyperparameters in "main()".
You can then run this file with the command: "python nn.py" in your terminal.
The program will automatically check your gradient implementation before start.
The program will print out the training progress, and it will display the
training curve by the end. You can optionally save the model by uncommenting
the lines in "main()".
"""

"""
@author: Luis Alvaro Correia - Std No.1006508566
"""

from utils import load_data, load_model, load_stats, save, display_plot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image
import pandas as pd
import seaborn as sns

def init_nn(num_inputs, num_hiddens, num_outputs):
    """ Initializes neural network's parameters.
    :param num_inputs: Number of input units
    :param num_hiddens: List of two elements, hidden size for each layers.
    :param num_outputs: Number of output units
    :return: A dictionary of randomly initialized neural network weights.
    """
    
    W1 = 0.1 * np.random.randn(num_inputs, num_hiddens[0])
    W2 = 0.1 * np.random.randn(num_hiddens[0], num_hiddens[1])
    W3 = 0.01 * np.random.randn(num_hiddens[1], num_outputs)
    b1 = np.zeros((num_hiddens[0]))
    b2 = np.zeros((num_hiddens[1]))
    b3 = np.zeros((num_outputs))
    model = {
        "W1": W1,
        "W2": W2,
        "W3": W3,
        "b1": b1,
        "b2": b2,
        "b3": b3
    }
    return model


def affine(x, w, b):
    """ Computes the affine transformation.
    :param x: Inputs (or hidden layers)
    :param w: Weight
    :param b: Bias
    :return: Outputs
    """
    y = x.dot(w) + b
    return y


def affine_backward(grad_y, x, w):
    """ Computes gradients of affine transformation.
    Hint: you may need the matrix transpose np.dot(A, B).T = np.dot(B, A) and (A.T).T = A
    :param grad_y: Gradient from upper layer
    :param x: Inputs from the hidden layer
    :param w: Weights
    :return: A tuple of (grad_h, grad_w, grad_b)
        WHERE
        grad_x: Gradients wrt. the inputs/hidden layer.
        grad_w: Gradients wrt. the weights.
        grad_b: Gradients wrt. the biases.
    """
    #####################################################################
    # TODO:                                                             #
    # Complete the function to compute the gradients of affine          #
    # transformation.                                                   #
    #####################################################################

    # Calculating the gradients
    grad_x = grad_y.dot(w.T).reshape(x.shape)
    grad_w = x.reshape(x.shape[0], w.shape[0]).T.dot(grad_y)
    grad_b = np.sum(grad_y, axis=0) 
    #####################################################################
    #                       END OF YOUR CODE                            #
    #####################################################################
    return grad_x, grad_w, grad_b


def relu(x):
    """ Computes the ReLU activation function.
    :param z: Inputs
    :return: Activation of x
    """
    return np.maximum(x, 0.0)


def relu_backward(grad_y, x):
    """ Computes gradients of the ReLU activation function wrt. the unactivated inputs.
    :param grad_y: Gradient of the activation.
    :param x: Inputs
    :return: Gradient wrt. x
    """
    #####################################################################
    # TODO:                                                             #
    # Complete the function to compute the gradients of relu.           #
    #####################################################################

    # Calculating the gradient
    grad_x = grad_y * (x > 0)
    #####################################################################
    #                       END OF YOUR CODE                            #
    #####################################################################
    return grad_x


def softmax(x):
    """ Computes the softmax activation function.
    :param x: Inputs
    :return: Activation of x
    """
    return np.exp(x) / np.exp(x).sum(axis=1, keepdims=True)


def nn_forward(model, x):
    """ Runs a forward pass.
    :param model: Dictionary of all the weights.
    :param x: Input to the network.
    :return: Dictionary of all intermediate variables.
    """
    z1 = affine(x, model["W1"], model["b1"])
    h1 = relu(z1)
    z2 = affine(h1, model["W2"], model["b2"])
    h2 = relu(z2)
    y = affine(h2, model["W3"], model["b3"])
    var = {
        "x": x,
        "z1": z1,
        "h1": h1,
        "z2": z2,
        "h2": h2,
        "y": y
    }
    return var


def nn_backward(model, err, var):
    """ Runs the backward pass.
    :param model: Dictionary of all the weights.
    :param err: Gradients to the output of the network.
    :param var: Intermediate variables from the forward pass.
    :return: None
    """
    dE_dh2, dE_dW3, dE_db3 = affine_backward(err, var["h2"], model["W3"])
    dE_dz2 = relu_backward(dE_dh2, var["z2"])
    dE_dh1, dE_dW2, dE_db2 = affine_backward(dE_dz2, var["h1"], model["W2"])
    dE_dz1 = relu_backward(dE_dh1, var["z1"])
    _, dE_dW1, dE_db1 = affine_backward(dE_dz1, var["x"], model["W1"])
    model["dE_dW1"] = dE_dW1
    model["dE_dW2"] = dE_dW2
    model["dE_dW3"] = dE_dW3
    model["dE_db1"] = dE_db1
    model["dE_db2"] = dE_db2
    model["dE_db3"] = dE_db3
    return


def nn_update(model, eta):
    """ Update NN weights.
    :param model: Dictionary of all the weights.
    :param eta: Learning rate
    :return: None
    """
    #####################################################################
    # TODO:                                                             #
    # Complete the function to update the neural network's parameters.  #
    # Your code should look as follows                                  #
    # model["W1"] = ...                                                 #
    # model["W2"] = ...                                                 #
    # ...                                                               #
    #####################################################################
    
    # Updating the NN
    model["W1"] -= model["dE_dW1"]*eta
    model["W2"] -= model["dE_dW2"]*eta
    model["W3"] -= model["dE_dW3"]*eta
    model["b1"] -= model["dE_db1"]*eta
    model["b2"] -= model["dE_db2"]*eta
    model["b3"] -= model["dE_db3"]*eta
    
    #####################################################################
    #                       END OF YOUR CODE                            #
    #####################################################################
    return


def train(model, forward, backward, update, eta, num_epochs, batch_size):
    """ Trains a simple MLP.
    :param model: Dictionary of model weights.
    :param forward: Forward prop function.
    :param backward: Backward prop function.
    :param update: Update weights function.
    :param eta: Learning rate.
    :param num_epochs: Number of epochs to run training for.
    :param batch_size: Mini-batch size, -1 for full batch.
    :return: A tuple (train_ce, valid_ce, train_acc, valid_acc)
        WHERE
        train_ce: Training cross entropy.
        valid_ce: Validation cross entropy.
        train_acc: Training accuracy.
        valid_acc: Validation accuracy.
    """
    inputs_train, inputs_valid, inputs_test, target_train, target_valid, \
        target_test = load_data("data\\toronto_face.npz")
    rnd_idx = np.arange(inputs_train.shape[0])
    
    train_ce_list = []
    valid_ce_list = []
    train_acc_list = []
    valid_acc_list = []
    num_train_cases = inputs_train.shape[0]
    if batch_size == -1:
        batch_size = num_train_cases
    num_steps = int(np.ceil(num_train_cases / batch_size))
    for epoch in range(num_epochs):
        np.random.shuffle(rnd_idx)
        inputs_train = inputs_train[rnd_idx]
        target_train = target_train[rnd_idx]
        for step in range(num_steps):
            # Forward pass.
            start = step * batch_size
            end = min(num_train_cases, (step + 1) * batch_size)
            x = inputs_train[start: end]
            t = target_train[start: end]

            var = forward(model, x)
            prediction = softmax(var["y"])

            train_ce = -np.sum(t * np.log(prediction)) / float(x.shape[0])
            train_acc = (np.argmax(prediction, axis=1) ==
                         np.argmax(t, axis=1)).astype("float").mean()
            print(("Epoch {:3d} Step {:2d} Train CE {:.5f} "
                   "Train Acc {:.5f}").format(
                epoch, step, train_ce, train_acc))

            # Compute error.
            error = (prediction - t) / float(x.shape[0])

            # Backward prop.
            backward(model, error, var)

            # Update weights.
            update(model, eta)

        valid_ce, valid_acc = evaluate(
            inputs_valid, target_valid, model, forward, batch_size=batch_size)
        print(("Epoch {:3d} "
               "Validation CE {:.5f} "
               "Validation Acc {:.5f}\n").format(
            epoch, valid_ce, valid_acc))
        train_ce_list.append((epoch, train_ce))
        train_acc_list.append((epoch, train_acc))
        valid_ce_list.append((epoch, valid_ce))
        valid_acc_list.append((epoch, valid_acc))
    display_plot(train_ce_list, valid_ce_list, "Cross Entropy", number=0)
    display_plot(train_acc_list, valid_acc_list, "Accuracy", number=1)

    train_ce, train_acc = evaluate(
        inputs_train, target_train, model, forward, batch_size=batch_size)
    valid_ce, valid_acc = evaluate(
        inputs_valid, target_valid, model, forward, batch_size=batch_size)
    test_ce, test_acc = evaluate(
        inputs_test, target_test, model, forward, batch_size=batch_size)
    print("CE: Train %.5f Validation %.5f Test %.5f" %
          (train_ce, valid_ce, test_ce))
    print("Acc: Train {:.5f} Validation {:.5f} Test {:.5f}".format(
        train_acc, valid_acc, test_acc))

    stats = {
        "train_ce": train_ce_list,
        "valid_ce": valid_ce_list,
        "test_ce": test_ce,           # Included for reporting
        "train_acc": train_acc_list,
        "valid_acc": valid_acc_list,
        "test_acc": test_acc          # Included for reporting
    }

    return model, stats


def evaluate(inputs, target, model, forward, batch_size=-1):
    """ Evaluates the model on inputs and target.
    :param inputs: Inputs to the network
    :param target: Target of the inputs
    :param model: Dictionary of network weights
    :param forward: Function for forward pass
    :param batch_size: Batch size
    :return: A tuple (ce, acc)
        WHERE
        ce: cross entropy
        acc: accuracy
    """
    num_cases = inputs.shape[0]
    if batch_size == -1:
        batch_size = num_cases
    num_steps = int(np.ceil(num_cases / batch_size))
    ce = 0.0
    acc = 0.0
    for step in range(num_steps):
        start = step * batch_size
        end = min(num_cases, (step + 1) * batch_size)
        x = inputs[start: end]
        t = target[start: end]
        prediction = softmax(forward(model, x)["y"])
        ce += -np.sum(t * np.log(prediction))
        acc += (np.argmax(prediction, axis=1) == np.argmax(
            t, axis=1)).astype("float").sum()
    ce /= num_cases
    acc /= num_cases
    return ce, acc


def check_grad(model, forward, backward, name, x):
    """ Check the gradients.
    """
    # np.random.seed(0)
    var = forward(model, x)
    loss = lambda y: 0.5 * (y ** 2).sum()
    grad_y = var["y"]
    backward(model, grad_y, var)
    grad_w = model["dE_d" + name].ravel()
    w_ = model[name].ravel()
    eps = 1e-7
    grad_w_2 = np.zeros(w_.shape)
    check_elem = np.arange(w_.size)
    np.random.shuffle(check_elem)
    # Randomly check 20 elements.
    check_elem = check_elem[:20]
    for ii in check_elem:
        w_[ii] += eps
        err_plus = loss(forward(model, x)["y"])
        w_[ii] -= 2 * eps
        err_minus = loss(forward(model, x)["y"])
        w_[ii] += eps
        grad_w_2[ii] = (err_plus - err_minus) / 2. / eps
    np.testing.assert_almost_equal(grad_w[check_elem], grad_w_2[check_elem],
                                   decimal=3)

def prtStats(stats, hyperP, ProcID):
    """ Process Summary Statistics
    """
    print('\nProcessing Report...\n')
    
    L = hyperP['num_epochs']

    f=open('Summary_NN_'+ProcID+'.prn','w')
    f.write("\n------ Processing Summary (%s) -----------\n" % ProcID)
    f.write('\nHyperparameters: ' % hyperP)
    f.write("\n>>> num_hiddens: [%s]" % ', '.join(map(str, hyperP['num_hiddens'])))
    f.write('\n>>> num_epochs: %d' % hyperP['num_epochs'])
    f.write('\n>>> eta: %.4f' % hyperP['eta'])
    f.write('\n>>> batch-size: %d' % hyperP['batch_size'])
    
    
    f.write("\n\nTraining accuracy is %.5f" % stats[1]["train_acc"][L-1][1])
    f.write("\nValidation accuracy is %.5f" % stats[1]["valid_acc"][L-1][1])
    f.write("\nTest accuracy is %.5f\n" % stats[1]["test_acc"])

    f.write("\n\nTraining CE is %.5f" % stats[1]["train_ce"][L-1][1])
    f.write("\nValidation CE is %.5f" % stats[1]["valid_ce"][L-1][1])
    f.write("\nTest CE is %.5f\n" % stats[1]["test_ce"])
    f.write("\n--------------------------------------------------\n")
    
    # Listing of TRAINING Errors for iteration 1-50 and 951-1000
    NTrain = len(stats[1]["train_acc"])
    f.write("\n------ Listing of Errors from TRAINING procedure (%d Iterations) ------\n" % NTrain)
    f.write("\n>>> Note: Only first 50 and last 50 errors will be printed due to\n")
    f.write("          limitation on Crowdmark to manage high no. of pages in \n")
    f.write("          uploaded PDF format.\n\n")
    
    dtrain = {'First50':np.ones(50)-list(zip(*stats[1]["train_acc"][0:50]))[1],
              'Last50':np.ones(50)-list(zip(*stats[1]["train_acc"][NTrain-50:]))[1]}
    dftrain = pd.DataFrame(dtrain)
    f.write(dftrain.to_string(header=True))

    NValid = len(stats[1]["valid_acc"])
    f.write("\n\n------ Listing of Errors from TESTING procedure (%d Iterations)------\n" % NValid)
    f.write("\n>>> Note: Only first 50 and last 50 errors will be printed due to\n")
    f.write("          limitation on Crowdmark to manage high no. of pages in \n")
    f.write("          uploaded PDF format.\n\n")
    
    dvalid = {'First50':np.ones(50)-list(zip(*stats[1]["valid_acc"][0:50]))[1],
              'Last50':np.ones(50)-list(zip(*stats[1]["valid_acc"][NValid-50:]))[1]}
    dfvalid = pd.DataFrame(dvalid)
    f.write(dfvalid.to_string(header=True))

    f.close()
    
    save_plot(stats[1]["train_ce"], stats[1]["valid_ce"], 
              "Cross Entropy", "CE-"+ProcID+".png", number=2)
    save_plot(stats[1]["train_acc"], stats[1]["valid_acc"], 
              "Accuracy", "ACC-"+ProcID+".png", number=3)
    
def prtImages (model, stats, hyperP, PID):
    """ Prints images within a criteria
    """
    threshold = .5  # Identify predictions with confidence below 50%
    inputs_train, inputs_valid, inputs_test, target_train, target_valid, \
        target_test = load_data("data\\toronto_face.npz")
    print('Images loaded :)')
    
    # re-Calculates Prediction of Test Images
    prediction = softmax(nn_forward(model, inputs_test)["y"])
    L0 = np.argmax(prediction, axis=1)    # List of predictions
    res0 = prediction[np.arange(prediction.shape[0]), L0]  # List of Confidence NN
    
    # Calculates the Test Errors - generates DF with data
    ErrTest = np.nonzero(target_test.argmax(axis=1)-prediction.argmax(axis=1))[0]
    dErrTest = {'Img No.':ErrTest,
              'Label':target_test[ErrTest].argmax(axis=1)+1,
              'ConfNN':res0[ErrTest],
              'Classif.':prediction[ErrTest].argmax(axis=1)+1}
    dfErrTest = pd.DataFrame(dErrTest)
    np.savetxt(r'dfErrTestNN'+PID+'.prn', dfErrTest.values, fmt='%03d %d %.5f %d', delimiter='\t')
    
    # Plot ScatterPlot of Classification vs. Label - Error Test
    plt.figure(figsize=(10,7))
    sns.relplot(x='Label', y='Classif.', size='ConfNN',
            sizes=(40, 400), alpha=.5, palette="muted",
            height=6, data=dfErrTest)    
    plt.savefig('ErrorTestGraph'+PID+'.png')
    plt.close()    

    # Calculates the Test Samples with lowest NN Confidence
    L1 = np.where(res0<threshold)[0]   # List of lowest confidence of NN
    dLowConf = {'Img No.':L1,
              'Label':target_test[L1].argmax(axis=1)+1,
              'ConfNN':prediction[np.arange(prediction.shape[0]), L0][L1],
              'Classif.':prediction[L1].argmax(axis=1)+1,
              'Correct':False}
    dLowConf['Correct'] = (dLowConf['Label']==dLowConf['Classif.'])
    dfLowConf = pd.DataFrame(dLowConf)
    np.savetxt(r'dfLowConfNN'+PID+'.prn', dfLowConf.values, fmt='%03d %d %.5f %d %s', delimiter='\t')

    # Plot ScatterPlot of Classification vs. Label - Lowest Confidence
    plt.figure(figsize=(10,7))
    sns.relplot(x='Label', y='Classif.', size='ConfNN',
            sizes=(40, 400), alpha=.5, palette="muted",
            height=6, data=dfLowConf)    
    plt.savefig('LowConfGraph'+PID+'.png')
    plt.close()    
    
    LS = dfLowConf.sort_values(by=['ConfNN'])['Img No.']
    # Print Latex Table for Overleaf
    print((dfLowConf.sort_values(by=['ConfNN'])).to_latex(index=False))
    
    # Plot Images with lowest confidence
    save_images(abs(inputs_test[LS,:]-1), 'imgFaceLowConf'+PID+'.png')
        

def plot_images(images, ax, ims_per_row=5, padding=5, digit_dimensions=(48, 48),
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
    
def save_plot(train, valid, y_label, pltname, number=0):
    """ Save Plot.
    :param train: Training statistics
    :param valid: Validation statistics
    :param y_label: Y-axis label of the plot
    :param number: The number of the plot
    :return: None
    """
    plt.figure(number, figsize=(10,7))
    plt.clf()
    train = np.array(train)
    valid = np.array(valid)
    plt.plot(train[:, 0], train[:, 1], "b", label="Train")
    plt.plot(valid[:, 0], valid[:, 1], "g", label="Validation")
    plt.xlabel("Epoch")
    plt.ylabel(y_label)
    plt.legend()
    plt.savefig(pltname)
    plt.close()    

    

def main():
    """ Trains a neural network.
    :return: None
    """
    PID = "RunQ2-2"
    model_file_name = "nn_model-"+PID+".npz"
    stats_file_name = "nn_stats-"+PID+".npz"
    loadedM = False

    # Hyper-parameters. Modify them if needed.
    num_hiddens = [16, 32]
    eta = 0.01
    num_epochs = 1000 # Number of iterations
    batch_size = 100 
    
    # Hyper-parameters
    hyperP = {
        "num_hiddens": num_hiddens,
        "eta": eta,
        "num_epochs": num_epochs,           
        "batch_size": batch_size         
    }

    # Input-output dimensions.
    num_inputs = 2304
    num_outputs = 7

    # Initialize model.
    np.random.seed(123)  # Included to enable reproducibility
    model = init_nn(num_inputs, num_hiddens, num_outputs)

    # New functions to implement 'load' Model and Stats due to problems with numpy 1.16.4
    if (loadedM):
        model = load_model(model_file_name)
        stats = load_stats(stats_file_name)
    
    # Check gradient implementation.
    print("Checking gradients...")
    np.random.seed(894)   # Included to enable reproducibility
    x = np.random.rand(10, 48 * 48) * 0.1
    check_grad(model, nn_forward, nn_backward, "W3", x)
    check_grad(model, nn_forward, nn_backward, "b3", x)
    check_grad(model, nn_forward, nn_backward, "W2", x)
    check_grad(model, nn_forward, nn_backward, "b2", x)
    check_grad(model, nn_forward, nn_backward, "W1", x)
    check_grad(model, nn_forward, nn_backward, "b1", x)

    # Train model.
    np.random.seed(6941)   # Included to enable reproducibility
    if (not loadedM):
        stats = train(model, nn_forward, nn_backward, nn_update, eta,
                  num_epochs, batch_size)
        save(model_file_name, model)
        save(stats_file_name, stats)
    
    #Print Statistics of the model
    prtStats(stats, hyperP, PID)

    #2.5 - Print Images of the model
    if (loadedM):
        prtImages(model, stats, hyperP, PID)

if __name__ == "__main__":
    main()
