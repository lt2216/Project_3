
# Command-line arguments: NREPL (1)

from sys import argv
import os
import numpy as np
from sklearn.metrics import mean_squared_error
import math
import numpy as np
from matplotlib import pyplot as plt
plt.switch_backend('agg')

from sklearn.metrics import accuracy_score


def compute_accuracy_score(path):                                           # Compute accuracy statistics and store with informative variable names.

    y_pred = np.load(path + "/y_pred.npy")
    y_pred = MLE(y_pred)
    y_true = np.load(path + "/y_truth.npy")

    score = accuracy_score(y_true, y_pred)

    return score

def MLE(y):
    y_mle = []
    for n in y:
        mle=0
        for lab,i in enumerate(n):
            mle += i*lab
        y_mle.append(int(round(mle)))
    return y_mle

def compute_SSE(path):                                           # Compute accuracy statistics and store with informative variable names.

    y_pred = np.load(path + "/y_pred.npy")
    y_pred = MLE(y_pred)
    y_truth = np.load(path + "/y_truth.npy")

    MSE = mean_squared_error(y_truth, y_pred)
    SSE = MSE*len(y_truth)
    #RMSE =  math.sqrt(MSE)
    #Pearson_coef = np.corrcoef(y_truth, y_pred)[0,1]

    return SSE #, Pearson_coef


##################### Loop through directory structure #####################


option = raw_input('Absolute values or comparisons? [1: Absolute; 2: Comparisons] : ')


#rootdir = "/Users/Luis/Desktop/MSc-Bioinformatics-and-Theoretical-Systems-Biology/PROJECTS/THESIS/Work/Data/test_preprocessing/NREPL=10000"  #argv[1]

rootdir = "."

epoch1 = []
epoch2 = []
epoch3 = []

##### Make a list of all the sorting options per epoch #####

for root, dir, file in os.walk(rootdir,topdown=False):          # get name of all the files in the file tree.
    for f in file:
        if "SORTING_41" in root and "y_pred" in f:
            #path = os.path.join(root, f)                        # create path for each y_pred or y_truth file.
            if "1epoch" in root:                                # store paths in arrays corresponding to their epoch (there will be a different plot for each epoch).
                epoch1.append(root)
            elif "2epoch" in root:
                epoch2.append(root)
            elif "3epoch" in root:
                epoch3.append(root)

labels_list = ["epoch1", "epoch2", "epoch3"]

acc_score_epoch1 = np.zeros(4)
acc_score_epoch2 = np.zeros(4)
acc_score_epoch3 = np.zeros(4)

SSE_Cat_loss_epoch1 = np.zeros(4)
SSE_Cat_loss_epoch2 = np.zeros(4)
SSE_Cat_loss_epoch3 = np.zeros(4)

for epoch in labels_list:                                  # Choose an epoch number
    for path in eval(epoch):                     # For each path in the array of the chosen epoch:
        if "no_ordering" in path:
            SSE = compute_SSE(path)
            eval("SSE_Cat_loss_" + epoch)[0] = SSE
            acc_score = compute_accuracy_score(path)
            eval("acc_score_" + epoch)[0] = acc_score

        elif "order_rows_and_cols" in path:
            SSE = compute_SSE(path)
            eval("SSE_Cat_loss_" + epoch)[3] = SSE
            acc_score = compute_accuracy_score(path)
            eval("acc_score_" + epoch)[3] = acc_score

        elif "order_rows" in path:
            SSE = compute_SSE(path)
            eval("SSE_Cat_loss_" + epoch)[1] = SSE
            acc_score = compute_accuracy_score(path)
            eval("acc_score_" + epoch)[1] = acc_score

        elif "order_cols" in path:
            SSE = compute_SSE(path)
            eval("SSE_Cat_loss_" + epoch)[2] = SSE
            acc_score = compute_accuracy_score(path)
            eval("acc_score_" + epoch)[2] = acc_score


##################### Produce plot #####################


arr1 = np.array((SSE_Cat_loss_epoch1, SSE_Cat_loss_epoch2, SSE_Cat_loss_epoch3), dtype=float)
arr2 = np.array((acc_score_epoch1, acc_score_epoch2, acc_score_epoch3), dtype=float)

if option == 1:
    arr2 *= 100

##### Calculate means #####

SSE_mean = np.mean(arr1, axis=0)
acc_score_mean = np.mean(arr2, axis=0)

##### Calculate standard deviations #####

sd1 = np.std(arr1, axis=0)
sd2 = np.std(arr2, axis=0)

if option == 2: #scale means and sds

    SSE_mean_minima = np.min(SSE_mean)
    acc_mean_minima = np.min(acc_score_mean)

    SSE_mean /= SSE_mean_minima
    acc_score_mean /= acc_mean_minima


    ##### Calculate and scale standard deviations #####

    sd1 /= np.square(SSE_mean_minima)
    sd2 /= np.square(acc_mean_minima)


##### Plot SSEs #####

fig1, ax1 = plt.subplots(figsize=(12, 8))

n_groups = 4
index = np.arange(n_groups)
bar_width = 0.23
opacity = 0.8

plt.grid(axis='y', linestyle='--')

feature1 = plt.bar(index, SSE_mean, bar_width, linewidth=0.8, alpha=opacity, color='tan', edgecolor = "k", ecolor='k', capsize=5, yerr=sd1)

plt.tight_layout()
file_name = "SSE_comparison_sorting.png"
plt.savefig(rootdir + "/" + file_name)


##### Plot accuracies #####

fig2, ax2 = plt.subplots(figsize=(12, 8))

n_groups = 4
index = np.arange(n_groups)
bar_width = 0.23
opacity = 0.8

plt.grid(axis='y', linestyle='--')

feature2 = plt.bar(index, acc_score_mean, bar_width, linewidth=0.8, alpha=opacity, color='tan', edgecolor = "k", ecolor='k', capsize=5, yerr=sd2)

plt.tight_layout()
file_name = "Accuracy_comparison_sorting.png"
plt.savefig(rootdir + "/" + file_name)
