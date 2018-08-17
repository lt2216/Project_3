#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Carries out convolutional neural network analysis
"""

from sys import argv
import shutil

import gzip
import glob
import os
import shutil
import itertools
import PIL
import math
import numpy as np
import ntpath
import random
import csv
import json
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams #Note: Thhe version number has to be 2.2.2
import pydot
import graphviz
import tensorflow
import keras#Note: Thhe version number has to be 2.1.3
from keras import backend as K
import sklearn
from sklearn.model_selection import train_test_split
import scipy.stats
import scipy.stats as ss
import pandas as pd
import subprocess
import tqdm
import time
import subprocess
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sn
import pandas  as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma
from keras.models import Model
import operator
from mpl_toolkits.axes_grid1 import ImageGrid

################################################
matplotlib.rcParams.update({'figure.autolayout': True})
################################################

################################################
keras.backend.set_image_dim_ordering('tf')
#K.set_image_dim_ordering('tf')
keras.backend.set_image_data_format('channels_last')
#K.set_image_data_format('channels_last') #e' quello giusto per Theano
print("Image_data_format is " + keras.backend.image_data_format())
################################################
#from random import *

#for i in range(1,12):
#	print("##########", argv[i], "##########")


										##################################################################################################
										########################################## FUNCTIONS #############################################
										##################################################################################################


##########################################
#########FUNCTIONS FOR USER INPUT#########
##########################################

def listdir_nohidden(path):
	"""
	Gets list of files in a directory whilst ignoring hidden files

	Keyword Arguments:
		path (string) -- path to directory in which to get list of files

	Returns:
		 (string) -- list of directory contents with full paths

	"""
	return glob.glob(os.path.join(path, '*'))

def getfile(filetype):
	"""
	Checks if the input file/folder exists: if not asks for new file path until correct one found

	Keyword Arguments:
		filetype (string) -- type of file looking for, eg 'reference genome'
		if not folder.startswith('resize'):

	Returns:
		user_input (string) -- path to specifed files

	"""
	loop = True
	order = 'Path to the %s: ' %(filetype)
	error = 'Path entered does not exist. Enter the path again"'
	while loop:
		user_input = raw_input(order)
		user_input = user_input.strip()
		result = os.path.exists(user_input)
		if result == True:
			loop = False
			file_loc = user_input

		else:
			print error
	return file_loc

def inputNumber(parameter_name):
	"""
		Function to get threshold value and check if it is within user specified limits

		Keyword Arguments:
		parameter_name (string) -- name of parameter being specified

		Returns:
		value (int) -- user specified paramter value

		"""
	order = 'Write the %s paramter value as a positive integer: ' %(parameter_name)
	while True:
		try:
			userInput = int(raw_input(order))
			if userInput < 0:
				print("Not a positive integer! Try again.")
				continue
		except ValueError:
			print("Write postive integer in numerical form! Try again.")
			continue
		else:
			return userInput
			break

def options_menu(title, options):

	"""
		Creates menus of user specified options

		Keyword Arguments:
		title (string) -- Name of menu
		options (string) -- Specified user options

		Returns:
		options (string) -- Specified user options
		"""
	#print width_screen * "-"
	#print(title.center(width_screen))
	#	print '{:^{width_screen}}'.format(title,width_screen)
	#print width_screen * "-"
	for x in range(len(options)):
		print str(x+1) + ". {}".format(options[x])
	#print width_screen * "-"
	return(options)


##########################################
###########IMAGE PREPERATION##############
##########################################

def resize_image(parent_directory, img_columns,img_rows):
	"""
	Function to resize images

	Keyword Arguments:
		parent_directory (string) -- directory contaiing images to be resized
	"""
	#List of img files in the directory for the specific simulation file
	listing = os.listdir(parent_directory)
	#Number of img files in the directory for the specific simulation file
	num_samples = len(listing)
	#Create Directory to store the re-sized images
	dir_nam =   parent_directory + '/resize_col_'  + str(img_columns) + '_row_' + str(img_rows)
	if not os.path.exists(dir_nam):
		os.makedirs(dir_nam)
	#Open image and resize --> (img_rows,img_columns)
	for file in listing:
		if file.endswith('.bmp'):
			im = PIL.Image.open(parent_directory + '/' + file)
			im_resized = im.resize((img_columns,img_rows))
			im_resized.save(dir_nam + '/' + file)
	return;


def classification_report_csv(report):
	report_data = []
	lines = report.split('\n')
	for line in lines[2:-3]:
		row = {}
		row_data = line.split('      ')
		row['class'] = row_data[1]
		row['precision'] = float(row_data[2])
		row['recall'] = float(row_data[3])
		row['f1_score'] = float(row_data[4])
		row['support'] = float(row_data[5])
		report_data.append(row)
	dataframe = pd.DataFrame.from_dict(report_data)
	dataframe.to_csv(CNN_dir +'/'+'classification_report.csv', index = False)


import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt

def convert_label_to_pmf(y, nb_classes, class_label_dict):

	y_max = y												#Conserve the scalar vector
	y= keras.utils.np_utils.to_categorical(y,nb_classes)	#Lazy way of generating an array with the desired dimensions

	x = np.array(class_label_dict.keys()).astype(np.float)	#Vector of labels
	xU, xL = x + 0.5, x - 0.5
	for n,lab in enumerate(y_max):
		prob = ss.norm.cdf(xU, loc=lab, scale = 1.5) - ss.norm.cdf(xL, loc=lab, scale = 1.5)
		prob = prob / prob.sum() #normalize the probabilities so their sum is 1
		y[n] = prob

	return y


									##################################################################################################
									############################################# MAIN ###############################################
									##################################################################################################


##########################################
################GET FILE PATHS############
##########################################

#Get the location of the simulations
simulation_files = argv[1] #'./NREPL=4/Simulations_1epoch'
#Get the location of the results directory
results_directory =  argv[2] #'./NREPL=4/Images_1epoch/M_m/no_ordering'
#Using the path to the results directory we get the path to the simulation images
simulation_images = results_directory + '/simulation_data'

##########################################
###########DATA_SET INFORMATION###########
##########################################
'''
NOTE: We assume that the values of NITER, NCHROMS and NREF are constant between the simulation files
'''
#Open the index file
with open(simulation_files + '/data.json', 'r') as fp:
	data_index = json.load(fp)
#Determine which simulation files are currently in ananlysis: Deleted simulations are still kept in index file [marked as 'inactive']
active_dict = filter(lambda x: x.get('active') == 'active', data_index)
#Get set of the unqiue selection_coefficients from those active simulation files
selection_coefficients = set([d['SAA'] for d in active_dict if 'SAA' in d])
NITER = int(active_dict[0]['NITER'])
NCHROMS = int(active_dict[0]['NCHROM'])
NREF = int(active_dict[0]['NREF'])

##########################################
###########DETERMINE IMAGE SIZE###########
##########################################

#Setting the number of channels
channel_loop = True
while channel_loop:
	channel_option_1, channel_option_2 = options_menu("Channel Number", ['Black & White Images [1 Channel]','RGB Images [3 Channels]'])
	channel_options = int(argv[3]) #1
	if channel_options ==1:
		print "Choice: {} has been selected".format(channel_option_1)
		channels = 1
		channel_loop = False
	else:
		print("Wrong option selection. Enter your choice again")

#Setting the width and height of the images
Resize_loop_01 =True
while Resize_loop_01:
	option_1, option_2,option_3, option_4, option_5,option_6 = options_menu("Resize Options", ['Raw Image Dimensions','Square 32 x 32','Square 64 x 64','Square 96 x 96','Square 128 x 128','Square 256 x 256'])
	resize_options = int(argv[4]) #5
	if resize_options==1:
		print "Choice: {} has been selected".format(option_1)
		Resize_loop_01=False
		Resize_loop_02 =True
		#Open the file containing the dimension data
		with open(simulation_images + '/img_dimension.json', 'r') as fp:
			dimensions = json.load(fp)
		row_size = NCHROMS
		while Resize_loop_02:
			raw_option_1, raw_option_2, raw_option_3 = options_menu("Raw Dimension Options", ['Mean','Max','Min'])
			raw_resize_options = int(argv[5]) #1
			if raw_resize_options == 1:
				print "Choice: {} has been selected".format(raw_option_1)
				col_size= int(dimensions.get('mean'))
				Resize_loop_02 =False
			elif raw_resize_options == 2:
				print "Choice: {} has been selected".format(raw_option_2)
				col_size = int(dimensions.get('max'))
				Resize_loop_02 =False
			elif raw_resize_options == 3:
				print "Choice: {} has been selected".format(raw_option_3)
				col_size = int(dimensions.get('min'))
				Resize_loop_02 =False
			else:
				print("Wrong option selection. Enter your choice again")
	elif resize_options==2:
		print "Choice: {} has been selected".format(option_2)
		col_size = 32
		row_size = 32
		Resize_loop_01=False
	elif resize_options==3:
		print "Choice: {} has been selected".format(option_3)
		col_size = 64
		row_size = 64
		Resize_loop_01=False
	elif resize_options==4:
		print "Choice: {} has been selected".format(option_4)
		col_size = 96
		row_size = 96
		Resize_loop_01=False
	elif resize_options==5:
		print "Choice: {} has been selected".format(option_5)
		col_size = 128
		row_size = 128
		Resize_loop_01=False
	elif resize_options==6:
		print "Choice: {} has been selected".format(option_6)
		col_size = 256
		row_size = 256
		Resize_loop_01=False
	else:
		print("Wrong option selection. Enter your choice again")

##########################################
########DETERMINE TRAINING METHOD#########
##########################################

training_loop = True
while training_loop:
	option_1, option_2,option_3 = options_menu("Training Options", ['True Data Labels','Perturbation','Kullback–Leibler_Divergence'])
	training_options = int(argv[12]) #3
	if training_options == 1:
		options_training_name = option_1
		training_loop = False
	elif training_options == 2:
		options_training_name = option_2
		training_loop = False
	elif training_options == 3:
		options_training_name = option_3
		training_loop = False
	else:
		print("Incorrect option selected. Enter choice again")


########################################################################
#########################GET NAME TO SAVE MODEL#########################
########################################################################

order = 'Write the name in which to save the model'
userInput = "CNN_KL"		#str(raw_input(order))
model_file_name = os.path.splitext(userInput)[0]


#####################################################
#####CREATING DIRECTORIES FOR CNN RESULTS############
#####################################################

#Create directory for the CNN
if training_options == 3:
	CNN_dir = results_directory + '/CNN_' + argv[11] + '_PMF'
	#CNN_dir = results_directory + '/CNN_' + 'M_m' + '_PMF'

elif training_options == 1:
	CNN_dir = results_directory + '/CNN_' + argv[11] + '_Categorical'
	#CNN_dir = results_directory + '/CNN_' + 'M_m' + '_Categorical'


if not os.path.exists(CNN_dir):
		os.makedirs(CNN_dir)

class_label_dict = dict(zip(range(len(np.arange(0,410,10))+1), np.unique(map(int,np.arange(0,410,10)))))


############## Convert to 41 classes to 5 classes ############## For this, create new dic, convert previous label arrays.

#class_label_dict = dict(zip(range(len(np.arange(0,400,80))+1), np.unique(map(int,np.arange(0,400,80)))))  #For 5 classes.

label40_label5_dict = {}
counter = 0
old_labels = np.arange(41)
new_labels = np.zeros(len(old_labels)).astype(int)

new_labels[0:8] = 0
new_labels[8:16] = 1
new_labels[16:24] = 2
new_labels[24:32] = 3
new_labels[32:] = 4

for i,l in enumerate(old_labels):
	label40_label5_dict[l] = new_labels[i]


f = gzip.GzipFile(results_directory + '/train_im.npy.gz', 'r')
X_train = np.load(f)
X_train = X_train.astype(float)
f.close()
y_train = np.load(results_directory + '/train_par.npy')
#y_train = [label40_label5_dict[old_label] for old_label in y_train]

f = gzip.GzipFile(results_directory + '/val_im.npy.gz', 'r')
X_val = np.load(f)
X_val = X_val.astype(float)
f.close()
y_val = np.load(results_directory + '/val_par.npy')
#y_val = [label40_label5_dict[old_label] for old_label in y_val]

f = gzip.GzipFile(results_directory + '/test_im.npy.gz', 'r')
X_test = np.load(f)
X_test = X_test.astype(float)
f.close()
y_test = np.load(results_directory + '/test_par.npy')
#y_test = [label40_label5_dict[old_label] for old_label in y_test]


#X_train.shape[0]: The number of samples used for training
X_train = X_train.reshape(X_train.shape[0], row_size, col_size, channels)
X_val = X_val.reshape(X_val.shape[0], row_size, col_size, channels)
X_test = X_test.reshape(X_test.shape[0], row_size, col_size, channels)

nb_classes = len(class_label_dict.keys())
#nb_classes = len(np.unique(label40_label5_dict.values()))

if training_options == 1:

	y_train= keras.utils.np_utils.to_categorical(y_train,nb_classes)
	y_val= keras.utils.np_utils.to_categorical(y_val,nb_classes)
	y_test= keras.utils.np_utils.to_categorical(y_test,nb_classes)

elif training_options == 3:

							## Convert to pmf ###
	y_train = convert_label_to_pmf(y_train, nb_classes, class_label_dict)
	y_val = convert_label_to_pmf(y_val, nb_classes, class_label_dict)
	y_test = convert_label_to_pmf(y_test, nb_classes, class_label_dict)


############################################################################
############################ CNN ARCHITECTURE ##############################
############################################################################

						############################
						######## PARAMETERS ########
						############################

######batch_size --> number of samples used per gradient update (default = 32)
#######epochs --> number of epochs to train the model; an epoch is an iteration over the entire x and y data provided
batch_size_parameter = int(argv[6]) #4
#batch_size_parameter = 32
epochs = int(argv[7]) #2
#epochs = 20
filters = int(argv[8]) #32
#filters = 32
kernel_size = int(argv[9]) #32
#kernel_size = 3
pooling_size = int(argv[10]) #2
#pooling_size = 2

			#######################################################
			######## INITIATE SESSION FOR MULTICORE USAGE #########
			#######################################################

from keras.backend import tensorflow_backend as K

#sess=tensorflow.Session(config=tensorflow.ConfigProto(intra_op_parallelism_threads=12))
sess=tensorflow.Session(config=tensorflow.ConfigProto(intra_op_parallelism_threads=12))
K.set_session(sess)

						###################################
						######## CNN ARCHTIECTURE #########
						###################################

from keras.layers import LeakyReLU

#This model uses a linear sequential format
CNN = keras.models.Sequential()

#ROUND ONE OF PATTERN: CONV, MAX POOL, DROPOUT
#Declare the input layer
CNN.add(keras.layers.convolutional.Convolution2D(filters,(kernel_size, kernel_size), strides=(1,1),padding='same', data_format='channels_last', input_shape=(row_size,col_size,channels)))
CNN.add(LeakyReLU(alpha=0.1))
#Add pooling layer: Max Pooling method (Parameter reduction)
CNN.add(keras.layers.convolutional.MaxPooling2D(pool_size=(pooling_size,pooling_size)))
#Add dropout layer: regularises model to prevent over-fitting
CNN.add(keras.layers.core.Dropout(rate=0.7)) #0.7

#ROUND TWO OF PATTERN: CONV, MAX POOL, DROPOUT
#Add convolutional layer
CNN.add(keras.layers.convolutional.Convolution2D(filters,(kernel_size, kernel_size), strides=(1,1),padding='same', data_format='channels_last'))
CNN.add(LeakyReLU(alpha=0.1))
#Add pooling layer: Max Pooling method (Parameter reduction)
CNN.add(keras.layers.convolutional.MaxPooling2D(pool_size=(pooling_size,pooling_size)))
#Add dropout layer: regularises model to prevent over-fitting
CNN.add(keras.layers.core.Dropout(rate=0.7)) #0.7

#ROUND THREE OF PATTERN: CONV, MAX POOL, DROPOUT
#Add convolutional layer
CNN.add(keras.layers.convolutional.Convolution2D(filters,(kernel_size, kernel_size), strides=(1,1),padding='same', data_format='channels_last'))
CNN.add(LeakyReLU(alpha=0.1))
#Add pooling layer: Max Pooling method (Parameter reduction)
CNN.add(keras.layers.convolutional.MaxPooling2D(pool_size=(pooling_size,pooling_size)))
#Add dropout layer: regularises model to prevent over-fitting
CNN.add(keras.layers.core.Dropout(rate=0.7)) #0.7

#ADDING A FULLY CONNECTED LATER
#Weights of previous layers are flattened(made 1D) before passing to the fully connected dense layer
CNN.add(keras.layers.core.Flatten())
CNN.add(keras.layers.core.Dense(128))
CNN.add(LeakyReLU(alpha=0.1))
#Add dropout layer: regularises model to prevent over-fitting
#CNN.add(keras.layers.core.Dropout(rate=0.5))
#Add output layer
CNN.add(keras.layers.core.Dense(nb_classes, activation='softmax'))
#Summary of the CNN architecture
CNN.summary()

#IMAGE OF THE MODEL
keras.utils.plot_model(CNN, to_file=CNN_dir +'/'+'model.png')

						###########################
						######## COMPILE ##########
						###########################

#Compiling the neural network using KL divergence as a loss function.
#Adams optimizer is still used
CNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

						###########################
						######## TRAINING #########
						###########################

#Fitting the model as before
hist = CNN.fit(X_train, y_train, batch_size=batch_size_parameter, epochs=epochs, verbose=1, validation_data=(X_val, y_val))

					#######################################
					############# SAVE MODEL ##############
					#######################################

CNN.save(CNN_dir +'/'+ model_file_name +'.h5')

#Save the class_label association
with open(CNN_dir +'/'+ model_file_name + "_" + options_training_name + '_class_label_association.json', 'w') as fp:
	json.dump(class_label_dict, fp)

						############################
						######## PLOT TRAINING #####
						############################

#Set the proportion
golden_size = lambda width: (width, 2. * width / (1 + np.sqrt(5)))
fig, ax = plt.subplots(figsize=golden_size(6))
#Create dataframe of history object
hist_df = pd.DataFrame(hist.history)
#Plot the dataframe
hist_df.plot(ax=ax)
#Set the x and y labels
ax.set_ylabel('Accuracy/Loss')
ax.set_xlabel('# epochs')
#ax.set_ylim(.99*hist_df[1:].values.min(), 1.1*hist_df[1:].values.max())
plt.savefig(CNN_dir +'/' + model_file_name + "_" +  options_training_name + "_"  + "loss_acc.pdf")


							######################
							###### TESTING #######
							######################

#Evaluate the model using unseen testing data
score = CNN.evaluate(X_test ,y_test, batch_size = None, verbose=1)
print('Test loss:', score[0])
print('Test accuracy:', score[1])

				#########################################
				############### CNN PLOTS ###############
				#########################################



			######## CNN LOSSES AND ACCURACIES:#############


train_loss=hist.history['loss']
val_loss=hist.history['val_loss']
train_acc=hist.history['acc']
val_acc=hist.history['val_acc']
xc=range(epochs)
x_axis =np.zeros(len(xc))
for x,i in enumerate (xc):
	x_axis[i]=x+1

rcParams['axes.titlepad'] = 20
plt.figure(1,figsize=(7,5),facecolor='white')
plt.plot(x_axis,train_loss)
plt.plot(x_axis,val_loss)
plt.xlabel('Epoch', fontsize=12)
plt.ylabel('Loss', fontsize=12)
plt.title('Training loss and validation loss',fontsize=12)
plt.grid(True)
plt.legend(['Training loss','Validation loss'],fontsize=12)
plt.style.use(['classic'])
plt.savefig(CNN_dir +'/'+ "Loss.eps")


plt.figure(2,figsize=(7,5),facecolor='white')
plt.plot(x_axis,train_acc)
plt.plot(x_axis,val_acc)
plt.xlabel('Epoch',fontsize=12)
plt.ylabel('Accuracy',fontsize=12)
plt.title('Training accuracy and validation accuracy',fontsize=12)
plt.grid(True)
plt.legend(['Training accuracy','Validation accuracy'],fontsize=12,loc=4)
plt.style.use(['classic'])
plt.savefig(CNN_dir +'/'+ "Accuracy.eps")


#Get predictions using the testing data: this data we know the true label
Y_pred = CNN.predict(X_test,batch_size=None, verbose=1)
#Convert to class labels
y_pred = np.argmax(Y_pred, axis=1)
#The true labels are contained within the y_test_data array: we need to convert to class label
y_truth = np.argmax(y_test, axis=1)

np.save(CNN_dir + "/y_truth",y_truth)
np.save(CNN_dir + "/y_pred",y_pred)


			############## CONFUSION MATRIX ##############


#Get the confusion matrix: Note: the number of labels is dependent on the number of classes in the testing data
cm = confusion_matrix(y_truth,y_pred,labels = class_label_dict.keys())
cm = np.nan_to_num(cm.astype('float') / cm.sum(axis=1)[:,np.newaxis])
#Create a dataframe of the confusion matrix
df_cm = pd.DataFrame(cm,index= np.arange(len(class_label_dict.values())))

#PLOT THE DATAFRAME:

#Plot the confusion matrix
plt.figure(figsize = (10,7))
#Centralise the tick marks
tick_marks = np.arange(len(class_label_dict.values()))
cen_tick = tick_marks + 0.5
#Plot the dataframe
sn.set(font_scale=1.4)
sn.heatmap(df_cm, annot=False,annot_kws={"size": 12})
#Add the labels and the postions of the tick marks
plt.xticks(cen_tick, class_label_dict.values(), rotation=45, fontsize=8)
plt.yticks(cen_tick, class_label_dict.values(),rotation=45, fontsize=8)
plt.ylabel('True label')
plt.xlabel('Predicted label')
#Save the plot
plt.savefig(CNN_dir +'/'+ model_file_name + "_" + options_training_name + "Confusion Matrix.pdf",bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()


				############# CLASSIFICATION REPORT ##############


selection_coefficients_str = map(str, class_label_dict.values())
#This could be formatted to be in a given order
cr = classification_report(y_truth, y_pred, labels= class_label_dict.keys(), target_names = selection_coefficients_str )
classification_report_csv(cr)
