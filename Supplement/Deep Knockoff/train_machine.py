# 3/1/2024
# train the deep knockoff machines
import numpy as np
import cvxpy
import cvxopt
import sys
sys.path.append('/python_package/DeepKnockoffs/')
from DeepKnockoffs import KnockoffMachine
from DeepKnockoffs import GaussianKnockoffs

# input parameters for which data to train on
dir_data = "/result/deep/training_data/"
dir_machine = "/result/deep/machine/"
ind = int(sys.argv[1])

# iterate through each of the knockoff machine
names = ["control_neg", "case_neg", "control_pos", "case_pos"]

# compute the correlation
TrainingData = np.loadtxt(dir_data + names[ind] + ".txt", delimiter = " ")

# Compute the empirical covariance matrix of the training data
SigmaHat = np.cov(TrainingData, rowvar=False)
# Initialize generator of second-order knockoffs
second_order = GaussianKnockoffs(SigmaHat, mu=np.mean(TrainingData,0), method="sdp")
# Measure pairwise second-order knockoff correlations 
corr = (np.diag(SigmaHat) - np.diag(second_order.Ds)) / np.diag(SigmaHat)

p = TrainingData.shape[1]
n = TrainingData.shape[0]

# set up parameters for the machine
pars = dict()
# Number of epochs
pars['epochs'] = 50
# Number of iterations over the full data per epoch
pars['epoch_length'] = 20
# Data type, either "continuous" or "binary"
pars['family'] = "continuous"
# Dimensions of the data
pars['p'] = p
# Size of the test set
pars['test_size']  = int(0.1*n)
# Batch size
pars['batch_size'] = int(0.2*n)
# Learning rate
pars['lr'] = 0.01
# When to decrease learning rate (unused when equal to number of epochs)
pars['lr_milestones'] = [pars['epochs']]
# Width of the network (number of layers is fixed to 6)
pars['dim_h'] = int(10*p)
# Penalty for the MMD distance
pars['GAMMA'] = 1
# Penalty encouraging second-order knockoffs
pars['LAMBDA'] = 1
# Decorrelation penalty hyperparameter
pars['DELTA'] = 1
# Kernel widths for the MMD measure (uniform weights)
pars['alphas'] = [1.,2.,4.,8.,16.,32.,64.,128.]
# Target pairwise correlations between variables and knockoffs
pars['target_corr'] = corr

# train the machine on the data
checkpoint_name = dir_machine + names[ind] 
# Where to print progress information
logs_name = dir_machine  + names[ind]  + "_progress.txt"
print(logs_name)
# store the machine 
machine = KnockoffMachine(pars, checkpoint_name=checkpoint_name, logs_name=logs_name)
print("Fitting the knockoff machine..." + names[ind])

machine.train(TrainingData)





