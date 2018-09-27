
import gzip
import numpy as np
import keras
from matplotlib import pyplot as plt
import os
import scipy.stats as ss
import pandas as pd
from sklearn.metrics import mean_squared_error, accuracy_score

def MLE(y):
    y_mle = []
    for n in y:
        mle=0
        for lab,i in enumerate(n):
            #mle += i*lab
            mle += i*(float(class_label_dict[lab])+5)
        y_mle.append(round(mle,2))
    return y_mle


def convert_label_to_pmf(y, nb_classes, class_label_dict):

	y_max = y												#Conserve the scalar vector
	y= keras.utils.np_utils.to_categorical(y,nb_classes)	#Lazy way of generating an array with the desired dimensions

	x = np.array(class_label_dict.keys()).astype(np.float)	#Vector of labels
	xU, xL = x + 0.5, x - 0.5

	for n,lab in enumerate(y_max):
		prob = ss.norm.cdf(xU, scale=0.6,loc=lab) - ss.norm.cdf(xL, scale=0.6,loc=lab)
		prob = prob / prob.sum() #normalize the probabilities so their sum is 1
		y[n] = prob

	return y


options = raw_input('Choose sorting option [0:none; 1:rows; 2:columns; 3:both] : ')
option2 = raw_input('Absolute values or comparison? [0: absolute; 1: comparisons] : ')

rootdir = "."


options_dir = {}
options_dir['0'] = "no_ordering"
options_dir['1'] = "order_rows"
options_dir['2'] = "order_cols"
options_dir['3'] = "order_rows_and_cols"


for root, dir, file in os.walk(rootdir,topdown=False):                                   # get name of all the files in the file tree.
    for f in file:
        if options_dir[options] in root and "train_im" in f:
            if "1epoch" in root:
                epoch1_path = root
            elif "2epoch" in root:
                epoch2_path = root
            elif "3epoch" in root:
                epoch3_path = root

class_label_dict = dict(zip(range(len(np.arange(0,410,10))+1), np.unique(map(int,np.arange(0,410,10)))))


SSE_X1 = np.zeros(3)
acc_X1 = np.zeros(3)

SSE_X2 = np.zeros(3)
acc_X2 = np.zeros(3)

SSE_X3 = np.zeros(3)
acc_X3 = np.zeros(3)


for testing in range(1,4):

    #### Load testing data
    testing_dir = eval("epoch" + str(testing) + "_path")
    f = gzip.GzipFile(testing_dir + '/test_im.npy.gz', 'r')
    X_test = np.load(f)
    X_test = X_test.astype(float)
    f.close()

    X_test = X_test.reshape(X_test.shape[0], 128, 128, 1)
    y_test = np.load(testing_dir + '/test_par.npy')
    y_test = convert_label_to_pmf(y_test, 41, class_label_dict)
    y_truth = np.argmax(y_test, axis=1)

    for training in range(1,4):

        ##### Load model
        training_dir = eval("epoch" + str(training) + "_path")
        CNN = keras.models.load_model(training_dir + '/06CNN_SORTING_41_' + options_dir[options] + '/CNN_KL.h5')

        ##### Make the prediction
        Y_pred = CNN.predict(X_test,batch_size=None, verbose=1)

        ##### Compute performance metrics
        mle_pred = MLE(Y_pred)
        score = accuracy_score(y_truth, mle_pred)
        MSE = mean_squared_error(y_truth, mle_pred)
        SSE = MSE*len(y_truth)

        eval("SSE_X" + str(testing))[training-1] = SSE
        eval("acc_X" + str(testing))[training-1] = score



if option2 == 1:

    # Normalize SSEs against the correct training-testing model combination:

    SSE_X1 /= (SSE_X1[0])
    SSE_X2 /= (SSE_X2[1])
    SSE_X3 /= (SSE_X3[2])

    acc_X1 /= (acc_X1[0])
    acc_X2 /= (acc_X2[1])
    acc_X3 /= (acc_X3[2])


##################### Plot #####################

fig1, ax1 = plt.subplots(figsize=(10, 6))

n_groups = 3
index = np.arange(n_groups)
bar_width = 0.23
opacity = 0.8
plt.grid(axis='y', linestyle='--')

feature1 = plt.bar(index, SSE_X1, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='tan', label='Testing model 1')
feature2 = plt.bar(index + bar_width, SSE_X2, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='wheat', label='Testing model 2')
feature3 = plt.bar(index + bar_width*2, SSE_X3, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='brown', label='Testing model 3')


plt.xlabel('Training model', fontsize=20)
if option2 == 1:
    plt.ylabel('Within-group normalised SSE', fontsize=20)
else:
    plt.ylabel('SSE', fontsize=20)
plt.title('Robutness analysis - relative SSE', fontsize=23)
plt.xticks(index + bar_width, ('1', '2', '3'), fontsize=20)

plt.legend(prop={'size':18}, fontsize=20)

plt.tight_layout()
file_name = "SSE_Comparison_Robustness.png"
plt.savefig(rootdir + "/" + file_name)


##########

fig2, ax2 = plt.subplots(figsize=(10, 6))

n_groups = 4
index = np.arange(n_groups)
bar_width = 0.23
opacity = 0.8
plt.grid(axis='y', linestyle='--')

feature1 = plt.bar(index, acc_X1, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='tan', label='Testing model 1')
feature2 = plt.bar(index + bar_width, acc_X2, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='wheat', label='Testing model 2')
feature3 = plt.bar(index + bar_width*2, acc_X3, bar_width, linewidth=1, edgecolor = "k", alpha=opacity, color='brown', label='Testing model 3')


plt.xlabel('Training model', fontsize=20)
if option2 == 1:
    plt.ylabel('Within-group normalised accuracy', fontsize=20)
else:
    plt.ylabel('Accuracy', fontsize=20)
plt.title('Robustness analysis - relative Accuracy', fontsize=23)
plt.xticks(index + bar_width, ('1', '2', '3'), fontsize=20)

plt.tight_layout()
file_name = "Acc_Comparison_Robustness.png"
plt.savefig(rootdir + "/" + file_name)
