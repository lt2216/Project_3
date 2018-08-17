
# -*- coding: utf-8 -*-


from sys import argv # For passing arguments to the script.

import sklearn
from Bio import SeqIO
from PIL import Image
import numpy as np #Note: There version number of numpy has to be above 1.13.3, else the error ' unique() got an unexpected keyword argument 'axis'' will occur
import copy
import sys # Package contains functions that allow user input from the command line to be used in programme
import os
import os.path
import json
import platform
import ntpath
import gzip
import csv
import allel
#Packages required for progress bar
from tqdm import tqdm, trange
from time import sleep
#Getting the terminal screen size: used when formatting user input
import subprocess



										##################################################################################################
										########################################## FUNCTIONS #############################################
										##################################################################################################


##########################################
#########FUNCTIONS FOR USER INPUT#########
##########################################

def get_gene():
	"""
	Function to get the name of the gene being investigated:
	NOTE: Probably doesn't have to be a function

	Returns:
		gene_name (string) -- name of the gene under investigation

	"""
	order = 'Name of the gene under investigation: '
	gene_name  = raw_input(order)
	gene_name = gene_name.strip()
	return gene_name;

def getfile(filetype):
	"""
	Checks if the input file/folder exists: if not asks for new file path until correct one found

	Keyword Arguments:
		filetype (string) -- type of file looking for, eg 'reference genome'

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


def create_results_folder(results_folder_path, record_file_path,gene = 'na'):
	"""
 	Checks if specifieid folder exists and if not creates it.

	Keyword Arguments:
		results_folder_path (string) -- path to the specific folder
 		record_file_path (string) -- name of file being processed
 		gene (string) -- name of the gene under investigation(DEFAULT: na)

 	Returns:
 		save_results (string) -- path to specific results folder
	"""
	filename, file_extension = os.path.splitext(record_file_path)# extract the filename and the file extension
	file_base = ntpath.basename(filename)
	fasta_path = results_folder_path + '/' + 'real_data'
	sim_path = results_folder_path + '/' + 'simulation_data'
	gene_path = fasta_path + '/' + gene
	sim_file_path = sim_path + '/' + file_base

	if file_extension == '.fasta':
		save_results = gene_path
		if not os.path.exists(fasta_path):
			os.makedirs(fasta_path)
		if not os.path.exists(gene_path):
			os.makedirs(gene_path)
	else:
		if file_extension == '.txt':
			save_results = sim_file_path
			if not os.path.exists(sim_path):
				os.makedirs(sim_path)
			if not os.path.exists(sim_file_path):
				os.makedirs(sim_file_path)
		if file_extension == '.gz':
			save_results = os.path.splitext(sim_file_path)[0]
			if not os.path.exists(sim_path):
				os.makedirs(sim_path)
			if not os.path.exists(os.path.splitext(sim_file_path)[0]):
				os.makedirs(os.path.splitext(sim_file_path)[0])
	return save_results;

def get_threshold(threshold_type,lower_limit,upper_limit):
	"""
	Gets threshold value and checks if it is within user specified limits

	Keyword Arguments:

		lower_limit (float) -- lower limit of the theshold
		threshold_type (string) -- the type of threshold
		upper_limit (float) -- upper limit of the theshold

	Returns:
		user_input (float) -- user specified threshold value

	"""
	order = 'Type the %s threshold as a float: ' %(threshold_type)
	error_limits = 'The %s threshold entered is not within the limits of %d and %d.' %(threshold_type,lower_limit,upper_limit)
	user_input = float(input(order))
	if user_input > upper_limit or user_input < lower_limit:
		print error_limits
		get_threshold(threshold_type, lower_limit, upper_limit)  #Repeat the function recursively
	else:
		return user_input


#################################################
#########FUNCTIONS FOR  FASTA PROCESSING#########
#################################################

def frequencies (alleles,n_gen,n_positions,matrix):
	"""
	Calculates Allelic frequencies

	Keyword Arguments:
		alleles (list) -- list of different alleles
		n_gen (int) -- total number of genomes
		n_positions (int) -- total number of nucleotide positions in individual sequence
		matrix (matrix) -- matrix of allelic matrix positions

	Returns:
		freq (array) -- allele frequencies
		n_alleles (array) -- number of alleles
	"""
	freq = np.zeros((4,n_positions))
	for i,letter in enumerate(alleles):
		freq[i,:] = matrix[:,:,i].sum(axis=0)/n_gen
		Vfreq = copy.copy(freq)
		Vfreq[Vfreq>0]=1
		n_alleles = Vfreq.sum(axis=0)
	return [freq,n_alleles]

def delete_loci(n_positions,matrix,freq,n_alleles,ref_genome,positions):
	"""
	Deletes Allelic Positions: USED FOR READ DATA IN FASTA FILES

	Keyword Arguments:
		matrix (matrix) -- matrix of allelic matrix positions
		n_alleles (int) -- number of alleles
		n_gen (int) -- total number of genomes
		n_positions (int) -- total number of nucleotide positions in individual sequence
		positions () --
		ref_genome (list) -- referencce genome

	Returns:
		freq (array) -- allele frequencies
		mask () --
		matrix (matrix) -- matrix of ?
		n_alleles (int) -- number of alleles
		n_positions (int) -- total number of nucleotide positions in individual sequence
		ref_genome (list) -- reference genome
	"""
	mask = np.ones(n_positions, dtype=bool)
	mask[positions] = False
	matrix = matrix [:,mask,:]
	n_positions = matrix.shape[1]
	freq = freq [:,mask]
	n_alleles = n_alleles [mask]
	ref_genome = ref_genome[mask]
	return [mask,n_positions,matrix,freq,n_alleles,ref_genome]

def process_fasta(ref,records,alleles,p_order,populations,delete_d,delete_mono,delete_bi,delete_tri,delete_quadri):
	"""
	Processes the Fasta files: Calculates the allele frequencies etc.

	Keyword Arguments:
		alleles (matrix) -- matrix of allelic matrix positions
		delete_d (int) -- number of alleles
		delete_mono () --
		delete_bi () --
		delete_tri () --
		delete_quadri () --
		populations (list) -- population included in analysis
		p_order (list) -- ordered populations
		records (string) -- path to the recorded genome sequences
		ref (string) -- path to the reference genome sequence

	Returns:
		freq (array) -- Frequnecy of allleles at each loci
		matrix (matrix) -- Matrix of ?
		Mm (array) -- Major and Minor alleles
		n (int) -- Number of populations
		n_alleles (array) -- Number of alleles at each genome position
		n_positions (int) -- Number of loci in each genomes
		populations (list) -- Populations included in analysis
		ref_genome (array) -- Sequence of the reference genome
		tot_n_gen (list) -- Total number of genomes for each population
	"""
	#EXTRACT ALL THE POPULATIONS FROM THE .FASTA FILE
	if populations == []:
		pops = []
		for i in range(len(records)):
			pops.append(records[i].name)
			populations = np.unique(pops)
	#n --> Number of populations in the file/Number of populations specificed by the user
	n = len(populations)
	#ORDER THE POPULATIONS
	populations_ordered = [0]*n
	index = 0
	for population in (p_order):
		if (population in populations):
			populations_ordered[index] = population
			index += 1
	#populations --> List containing the ordered populations
	populations = populations_ordered
	#READ THE GENOMES OF THE POPULATIONS
	genomes = []
	#tot_n_gen --> List containing the number of indviduals in each population
	#(ordine dato da populations)
	tot_n_gen = [0]*n
	for p,population in enumerate(populations):
		#p_n_gen --> Counter for individuals of every populations
		p_n_gen = 0
		for i in range(len(records)):
			if records[i].name == population:
				p_n_gen += 1
				genome = []
				for j in range(len(records[i].seq)):
					genome.append(records[i].seq[j])
				genomes.append(genome)
		tot_n_gen[p] = p_n_gen
	#total_genomes --> List of genomes of all population, one after the other.
	total_genomes = np.array(genomes)
	#READ THE REFERENCE GENOME
	ref_genome = []
	for i in range (len(ref.seq)):
		ref_genome.append (ref.seq[i])
	ref_genome = np.array(ref_genome)
	#Create a matrix of allelic matrix positions:
	#the first layer in the third dimension refers to the A and has 1 where there is an A
	#n_gen --> total number of individuals (all populations together)
	n_gen = total_genomes.shape[0]
	n_positions = total_genomes.shape[1]
	matrix = np.zeros((n_gen,n_positions,4))
	for i,letter in enumerate(alleles):
		xy = np.where(total_genomes == letter)
		xyz = list(xy)
		z = np.ones((len(xy[0])),dtype=int)*i
		xyz.append(z)
		matrix[xyz] = 1
	#Calculation of allele frequencies
	freq, n_alleles = frequencies (alleles,n_gen,n_positions,matrix)
	#DELETE mono / bi / tri / allelic positions
	arguments = [delete_mono,delete_bi,delete_tri,delete_quadri]
	if delete_mono==delete_bi==delete_tri==delete_quadri==False:
		pass
	else:
		for i,arg in enumerate(arguments):
			if arg == True:
				positions = np.where(n_alleles==i+1)
				mask,n_positions,matrix,freq,n_alleles,ref_genome = delete_loci(n_positions,matrix,freq,n_alleles,ref_genome,positions)
	#CALCULATION MINOR AND MAIOR
	Mm = np.zeros((2,n_positions), dtype=np.str)
	freq1 = copy.copy(freq)
	for k in range(n_positions):
		max1 = freq[:,k].argmax(axis=0)
		freq1[max1,k] = 0
		s = freq1[:,k].sum(axis=0)
		if s == 0:
			max2 = max1
		else:
			max2 = freq1[:,k].argmax(axis=0)
		for i,letter in enumerate(alleles):
			if max1 == i:
				Mm [0,k] = letter
			if max2 == i:
				Mm [1,k] = letter
	return [populations,matrix,tot_n_gen,ref_genome,freq,Mm,n_positions,n,n_alleles]

#################################################
#########FUNCTIONS FOR  IMAGE PROCESSING#########
#################################################
def options_menu(title, options):

	"""
	Creates menus of user specified options

	Keyword Arguments:
		title (string) -- Name of menu
		options (string) -- Specified user options

	Returns:
		options (string) -- Specified user options
	"""

	for x in range(len(options)):
		print str(x+1) + ". {}".format(options[x])

	return(options)


def order_data(im_matrix):
	"""
	Sorts matrix containg sequence data

	Keyword Arguments:
		im_matrix (Array [2D]) -- Array containing sequence data

	Returns:
		b (Array [2D]) -- Sorted array containing sequence data
	"""
	#u: Sorted Unique arrays
	#index: Index of 'im_matrix' that corresponds to each unique array
	#count: The number of instances each unique array appears in the 'im_matrix'
	u,index,count = np.unique(im_matrix,return_index=True,return_counts=True,axis=0)
	#b: Intitialised matrix the size of the original im_matrix[where new sorted data will be stored]
	b = np.zeros((np.size(im_matrix,0),np.size(im_matrix,1)),dtype=int)
	#c: Frequency table of unique arrays and the number of times they appear in the original 'im_matrix'
	c = np.stack((index,count), axis= -1)
	# The next line sorts the frequency table based mergesort algorithm
	c = c[c[:,1].argsort(kind='mergesort')]
	pointer = np.size(im_matrix,0)-1
	for j in range(np.size(c,0)):
		for conta in range(c[j,1]):
			b [pointer,:] = im_matrix[c[j,0]]
			pointer -= 1
	return b

def delete_simulation(n_positions,matrix,freq,positions):
	"""
	Deletes Allelic Positions: USE FOR SIMULATED DATA

	Keyword Arguments:
		freq (Array) -- frequency of alleles at each loci of the gene under investigation
		matrix (matrix) -- matrix of allelic matrix positions
		n_positions (int) -- total number of nucleotide positions in individual sequence
		pos () --
		positions () --

	Returns:
		freq (Array) -- frequency of alleles at each loci of the gene under investigation
		matrix (matrix) -- matrix of allelic matrix positions
		n_positions (int) -- total number of nucleotide positions in individual sequence
		pos () --
		positions () --
	"""
	mask = np.ones(n_positions, dtype=bool)
	mask[positions[0]] = False
	matrix = matrix[:,mask]
	n_positions = matrix.shape[1]
	freq = freq[mask]
	return [matrix,n_positions,freq]

def full_colour (n_gen,n_positions,matrix,population, gene_images_path):
	"""
	DESCRIBE THE FUNCTION HERE

	Keyword Arguments:
		n_gen () --
		n_positions (int) -- number of loci in genome sequence
		population (string) -- population
		gene_images_path (string) -- Path to directory where images are to be stored

	"""
	if not os.path.exists(gene_images_path + '/full_color'):
		os.makedirs(gene_images_path + '/full_color')
	global c
	colours = [[100,0,0,0],[0,100,0,0],[0,0,100,0],[0,0,0,100]]
	colour_genomes = np.zeros((n_gen,n_positions,4))
	for i,colour in enumerate(colours):
		index = np.where(matrix[:,:,i]==1)
		colour_genomes[index] = colour
	colour_genomes_uint8 = np.uint8(colour_genomes)
	colour_genomes_im = Image.fromarray(colour_genomes_uint8, mode = 'CMYK').convert('RGB')
	string = gene_images_path  + "/full_color/" + population + ".bmp"
	colour_genomes_im.save (string)
	return;

def black_white (n_gen,n_positions,alleles,matrix,ref_genome,population,p_threshold,apply_threshold,gene_images_path):
	"""
	DESCRIBE THE FUNCTION HERE

	Keyword Arguments:
		alleles (List) -- Alleles to include in investigated
		apply_threshold (Boolean) -- Apply threshold (True/False)
		matrix (Matrix) --
		n_gen () --
		n_positions (int) -- Number of loci in the gene sequence
		population (string) --
		p_threshold (float) -- Threshold for ?
		ref_genome (array) -- Nucleotide sequence of the reference genome
		gene_images_path (string) -- Path to directory where images are to be stored

	"""
	if not os.path.exists(gene_images_path + '/black_white'):
		os.makedirs(gene_images_path + '/black_white')

	bw_genomes = np.zeros((n_gen,n_positions))
	freq_a = np.zeros((1,n_positions))
	for k in range(n_positions):
		for i,letter in enumerate(alleles):
			if ref_genome[k] == letter:
				#bw_genomes ha 1 dove l'allele è quello ancestrale
				bw_genomes[:,k] = matrix[:,k,i]
				freq_a[0,k] = bw_genomes[:,k].sum(axis=0)/n_gen
	if apply_threshold == True:
		#elimino le posizioni in cui la frequenza dell'allele ancestrale è minore di p_threshold
		positions = np.where(freq_a <= p_threshold)
		mask = np.ones(n_positions, dtype=bool)
		mask[positions[1]] = False
		bw_genomes = bw_genomes [:,mask]
		#n_p_positions = bw_genomes.shape[1]
	#ORDINO LE RIGHE
	#bw_genomes = order_data(bw_genomes)
	bw_genomes_uint8 = np.uint8(bw_genomes)
	bw_genomes_im = Image.fromarray (bw_genomes_uint8*255, mode = 'L')
	string = gene_images_path + '/black_white/' + population + ".bmp"
	bw_genomes_im.save (string)
	return;

def black_white_Mm (n_gen,n_positions,Mm,matrix,alleles,population,gene_images_path):
	"""
	Generation of black and white images

	Keyword Arguments:
		alleles (List) -- Alleles to include in investigated
		matrix (Matrix) --
		Mm (Array) -- Major and Minor alleles at each loci
		n_gen () --
		n_positions (int) -- Number of loci in the gene sequence
		population (string) --
		gene_images_path (string) -- Path to directory where images are to be stored

	"""
	if not os.path.exists(gene_images_path + '/black_white_Mm'):
		os.makedirs(gene_images_path + '/black_white_Mm')

	bw_Mm_genomes = np.zeros((n_gen,n_positions))
	for k in range(n_positions):
		for i,letter in enumerate(alleles):
			if Mm[0,k]==letter:
				#bw_Mm_genomes ha 1 dove l'allele e' quello maior
				bw_Mm_genomes[:,k]=matrix[:,k,i]
	#ORDINO LE RIGHE
	#bw_Mm_genomes = order_data(bw_Mm_genomes)
	bw_Mm_genomes_uint8 = np.uint8(bw_Mm_genomes)
	bw_Mm_genomes_im = Image.fromarray (bw_Mm_genomes_uint8*255, mode = 'L')
	string = gene_images_path + '/black_white_Mm/' + population + ".bmp"
	bw_Mm_genomes_im.save (string)
	return;


def colour_freq (freq,n_positions):
	"""
	DESCRIBE THE FUNCTION HERE

	Keyword Arguments:
		freq () --
		n_positions (int) -- Number of loci in the gene sequence

	Returns:
		colour_freq_uint8 () --

	"""
	freq = np.reshape(freq.T, (1,n_positions,4))
	colour_freq_uint8 = np.uint8(freq)
	return colour_freq_uint8


def Mm_freq (Mm,n_positions,freq):
	"""
	Calculates the frequencies of major and minor alleles

	Keyword Arguments:
		freq () --
		Mm () --
		n_positions (int) -- Number of loci in the gene sequence

	Returns:
		m_uint8 () --

	"""
	m = np.zeros((1,n_positions))
	for k in range(n_positions):
		for i,letter in enumerate(alleles):
			if Mm[0,k]==letter:
				#m e' tanto piu'bianco quanto piu' e' alta la frequenza del maior
				#che e' come dire e' tanto piu' nero quanto piu' e' alta la freq. del minor
				m[0,k] = np.rint(freq[i,k]*2.55)
	m_uint8 = np.uint8(m)
	return m_uint8


def Ad_freq (n_gen,n_positions,alleles,matrix,ref_genome,population):
	"""
	DESCRIBE THE FUNCTION HERE

	Keyword Arguments:
		alleles () --
		matrix () --
		n_gen () --
		n_positions (int) -- Number of loci in the gene sequence
		population (string) -- Population where sequenced individual came from
		ref_genome (array) -- Sequence of the reference genome

	Returns:
		freq_d_uint8 () --

	"""
	bw_genomes = np.zeros((n_gen,n_positions))
	freq_a = np.zeros((1,n_positions))
	for k in range(n_positions):
		for i,letter in enumerate(alleles):
			if ref_genome[k] == letter:
				#bw_genomes ha 1 dove l'allele è quello ancestrale
				bw_genomes[:,k] = matrix[:,k,i]
				freq_a[0,k] = bw_genomes[:,k].sum(axis=0)/n_gen
	#freq_d e' tanto piu' nero quanto piu' e' alta la freq. dell'allele derivato
	freq_d = np.rint(freq_a*100)
	freq_d = np.rint(freq_d*2.55)
	freq_d_uint8 = np.uint8(freq_d)
	return freq_d_uint8


def generate_images(populations,matrix,tot_n_gen,ref_genome,freq,Mm,n_positions,n,n_alleles,p_threshold,apply_threshold, gene_images_path):

	"""
	DESCRIBE THE FUNCTION HERE : CREATES IMAGES FROM REAL DATA IN FASTA FILES

	Keyword Arguments:
		apply_threshold (Boolean) -- Apply threshold (True/False)
		freq (array) -- Frequnecy of allleles at each loci
		gene_images_path (string) -- path to directory where images are saved
		matrix (matrix) -- Matrix of ?
		Mm (array) -- Major and Minor alleles
		n (int) -- Number of populations
		n_alleles (array) -- Number of alleles at each genome position
		n_positions (int) -- Number of loci in each genomes
		p_threshold (float) --  ?
		populations (list) -- Populations included in analysis
		ref_genome (array) -- Sequence of the reference genome
		tot_n_gen (list) -- Total number of genomes for each population
	"""

	for p,population in enumerate(populations):
		from_gen = sum(tot_n_gen [:p])
		to_gen = from_gen + tot_n_gen[p]
		n_gen = (to_gen - from_gen)
		p_matrix = matrix[from_gen:to_gen,:,:]
		p_freq,p_n_alleles = frequencies (alleles,n_gen,n_positions,p_matrix)
		p_freq = np.rint(p_freq*100)
		if population == 'STU':
			colour_freq_uint8 = colour_freq (p_freq,n_positions)
			poli_freq_uint8 = np.concatenate((poli_freq_uint8, colour_freq_uint8),axis=0)
			poli_freq_im = Image.fromarray(poli_freq_uint8, mode = 'CMYK').convert('RGB')
			string = gene_images_path  + "/poli_freq.bmp"
			poli_freq_im.save (string)

			m_uint8 = Mm_freq (Mm,n_positions,p_freq)
			poli_m_uint8 = np.concatenate((poli_m_uint8, m_uint8*0),axis=0)
			poli_m_im = Image.fromarray(poli_m_uint8, mode = 'L')
			string = gene_images_path  + "/poli_m.bmp"
			poli_m_im.save (string)

			freq_d_uint8 = Ad_freq (n_gen,n_positions,alleles,p_matrix,ref_genome,population)
			poli_d_uint8 = np.concatenate((poli_d_uint8, freq_d_uint8),axis=0)
			poli_d_im = Image.fromarray(poli_d_uint8, mode = 'L')
			string = gene_images_path  + "/poli_d.bmp"
			poli_d_im.save (string)
		if population == 'TSI':
			colour_freq_uint8 = colour_freq (p_freq,n_positions)
			poli_freq_uint8 = np.concatenate((poli_freq_uint8, colour_freq_uint8),axis=0)
			poli_freq_im = Image.fromarray(poli_freq_uint8, mode = 'CMYK').convert('RGB')
			string = gene_images_path  + "/poli_freq.bmp"
			poli_freq_im.save (string)

			m_uint8 = Mm_freq (Mm,n_positions,p_freq)
			poli_m_uint8 = np.concatenate((poli_m_uint8, m_uint8*0),axis=0)
			poli_m_im = Image.fromarray(poli_m_uint8, mode = 'L')
			string = gene_images_path  + "/poli_m.bmp"
			poli_m_im.save (string)

			freq_d_uint8 = Ad_freq (n_gen,n_positions,alleles,p_matrix,ref_genome,population)
			poli_d_uint8 = np.concatenate((poli_d_uint8, freq_d_uint8),axis=0)
			poli_d_im = Image.fromarray(poli_d_uint8, mode = 'L')
			string = gene_images_path  + "/poli_d.bmp"
			poli_d_im.save (string)
		if p==0:
			poli_freq_uint8 = colour_freq (p_freq,n_positions)
			poli_m_uint8 = Mm_freq (Mm,n_positions,p_freq)
			poli_d_uint8 = Ad_freq (n_gen,n_positions,alleles,p_matrix,ref_genome,population)
		elif p == n-1:
			colour_freq_uint8 = colour_freq (p_freq,n_positions)
			poli_freq_uint8 = np.concatenate((poli_freq_uint8, colour_freq_uint8),axis=0)
			poli_freq_im = Image.fromarray(poli_freq_uint8, mode = 'CMYK').convert('RGB')
			string = gene_images_path  + "/poli_freq.bmp"
			poli_freq_im.save (string)

			m_uint8 = Mm_freq (Mm,n_positions,p_freq)
			poli_m_uint8 = np.concatenate((poli_m_uint8, m_uint8),axis=0)
			poli_m_im = Image.fromarray(poli_m_uint8, mode = 'L')
			string = gene_images_path  + "/poli_m.bmp"
			poli_m_im.save (string)

			freq_d_uint8 = Ad_freq (n_gen,n_positions,alleles,p_matrix,ref_genome,population)
			poli_d_uint8 = np.concatenate((poli_d_uint8, freq_d_uint8),axis=0)
			poli_d_im = Image.fromarray(poli_d_uint8, mode = 'L')
			string = gene_images_path  + "/poli_d.bmp"
			poli_d_im.save (string)
		else:
			colour_freq_uint8 = colour_freq (p_freq,n_positions)
			poli_freq_uint8 = np.concatenate((poli_freq_uint8, colour_freq_uint8),axis=0)
			m_uint8 = Mm_freq (Mm,n_positions,p_freq)
			poli_m_uint8 = np.concatenate((poli_m_uint8, m_uint8),axis=0)
			freq_d_uint8 = Ad_freq (n_gen,n_positions,alleles,p_matrix,ref_genome,population)
			poli_d_uint8 = np.concatenate((poli_d_uint8, freq_d_uint8),axis=0)

		full_colour (n_gen,n_positions,p_matrix,population, gene_images_path)
		black_white_Mm(n_gen,n_positions,Mm,p_matrix,alleles,population,gene_images_path)
		black_white (n_gen,n_positions,alleles,p_matrix,ref_genome,population,p_threshold,apply_threshold,gene_images_path)
	return;

def image_simulation(path1,path2,S, N, file_name, NCHROMS, threshold, apply_threshold,sort,maj_min,img_columns,img_rows,coef_label_dic):
	"""
	Generates images from iterations of simulation files
		- Deals with both txt files and gzip txt files
		- Calculates summary statistics for each iteration

	Keyword Arguments:
		apply_threshold (Boolean) -- Whether or not to apply p-threshold
		col_order (Boolean) -- Whether or not to order the columns
		file_name (string) -- The name of the simulation file being processed( either txt or txt.gz)
		NCHROMS (int) --
		N (int) -- N parameter of simulation
		n_alleles (array) -- Number of alleles at each genome position
		path1 (string) -- Path to directory where the simulation files exist
		path2 (string) -- Path to directory where produced image should be stored
		S (float) -- Selection Co-efficient of simulation
		threshold (float) -- Threshold value ?
		row_order (Boolean) -- Whether or not to order the columns
		maj_min (Boolean) -- Whether or not to colour my major/minor alleles

	Returns:
		simulation_error (list) -- List of erronous simulation files and the iteration with error in
		statistics_list (list) -- List of dictionaries containing summary statitics of simulations

	"""
	global once
	global nsl
	simulation_error = []
	statistics_list = []
	dim = []

	##################################################
	#############OPENING THE SIMULATION FILES#########
	##################################################

	#Suffix of g_zip files (Compressed)
	gzip_suffix = ".gz"
	#Suffix of txt files (Uncompressed)
	txt_suffix = ".txt"
	#we import and open the file
	if file_name.endswith(gzip_suffix):
		with gzip.open(path1 + file_name, 'rb') as f:
			file = f.read()
		if type(file) == str:
			#gzip files might need to be processed to be in correct format
			file = file.splitlines(True)
	elif file_name.endswith(txt_suffix):
		file = open(path1 + file_name).readlines()

	##################################################
	##########INDEXING THE FILES BY INTERATION########
	##################################################

	#we look for the caracter // inside the file
	find = []
	for i, string in enumerate(file):
		if string == '//\n':
			find.append(i+3)

	##################################################
	###GENERATE ONE IMAGE PER SIMULATION ITERATION####
	##################################################

	for ITER, pointer in enumerate(find):
		try:
			###########################
			####CREATE CHROM MATRIX####
			###########################

			n_columns = len(list(file[pointer]))-1
			croms = np.zeros((NCHROMS,n_columns),dtype=int)
			for j in range(NCHROMS):
				f = list(file[pointer + j])
				del f[-1]
				position_it = file[pointer - 1].split()
				del position_it[0]
				position_it = np.array(position_it, dtype='float')
				position_it = position_it*N
				F = np.array(f,dtype=int)
				if j == 0:
					crom_array = F
				else:
					crom_array = np.vstack((crom_array,F))
				croms[j,:]=F
			n_pos = np.size(croms,1)

			###########################
			#####APPLY THRESHOLD#######
			###########################

			if apply_threshold == True:
				#Count the number of derived alleles at each position
				count = croms.sum(axis=0,dtype=float)
				#Calculate the frrequency of the drived allele for each position
				freq = count/float(NCHROMS)
				for i in range(n_pos):
					if freq[i] > 0.5:
						freq[i] = 1-freq[i]
				#freq is now a vector that contains the minor allele frequency for each position
				#we delete the positions in which the minor allele frequency is <= threshold
				positions = np.where(freq<=threshold)
				croms,n_pos,freq = delete_simulation(n_pos,croms,freq,positions)

			###########################
			###COLOUR BY MAJOR/MINOR###
			###########################

			if maj_min == True:
				#Calculate the Major and the minor allele for each position of the matrix/array
				#Traspose the matrix/array
				transponse_array_croms = np.transpose(croms)
				#Record the Major and Minor allele for each allelic position
				maj_allele = []
				minor_allele = []
				for i in range(len(transponse_array_croms)):
					freq_data = np.unique(transponse_array_croms[i], return_counts = True)
					index_max =  np.argmax(freq_data[1])
					if index_max == 0:
						maj_allele.append(0)
						minor_allele.append(1)
					if index_max == 1:
						maj_allele.append(1)
						minor_allele.append(0)

				#Black and white image:
				#Simulation File: 0 = ancestrial, 1 = Derived (White encoded by 1, Black encoded by 0)
				#If the major allele is 0, we want to change 0 with 1 and vice verasa (1 = Major, 0 = Minor)
				#If the major allele is 1, no changes need to be made as 1 would by default be coded to be white

				matrix_maj_min_col = np.ones((n_pos,NCHROMS),dtype=int)
				for row in range(len(transponse_array_croms)):
					if maj_allele[row] == 1:
						matrix_maj_min_col[row,:] = transponse_array_croms[row]
					if maj_allele[row] == 0:
						matrix_maj_min_col[row,:] = matrix_maj_min_col[row,:] - transponse_array_croms[row]
				#Transpose the matrix so that the rows are the NCHROM and the columns are n_pos
				croms = np.transpose(matrix_maj_min_col)

			if maj_min == False:
				#Black and white image:
				#Simulation File: 0 = ancestrial, 1 = Derived (White encoded by 1, Black encoded by 0)
				#We want the opposite(ancestrial = white & derived = black) : hence we need to change 0 with 1 and vice versa before producing the image
				all1 = np.ones((NCHROMS,n_pos))
				croms = all1 - croms

			###########################
			####ORDER ROWS/COLUMNS#####
			###########################

			if sort == 2:
			#Sort the matrix by row (chromosome)
				croms = order_data(croms)

			if sort == 3:
			#Sort the matrix by column (genetic posistion)
				croms_transpose = croms.transpose()
				croms_transpose = order_data(croms_transpose)
				croms = croms_transpose.transpose()

			if sort == 4:
				#First: sort the matrix by row (chromosome)
				croms = order_data(croms)
				#Second: sort the matrix by column (genetic posistion)
				croms_transpose = croms.transpose()
				croms_transpose = order_data(croms_transpose)
				croms = croms_transpose.transpose()

			######################
			###IMAGE GENERATION###
			######################

			#Create image from the simulations
			bw_croms_uint8 = np.uint8(croms)
			bw_croms_im = Image.fromarray(bw_croms_uint8*255, mode = 'L')

			#####################
			## IMAGE RESIZINNG ##
			#####################

			bw_croms_im_resized = bw_croms_im.resize((img_columns,img_rows))
			dim.append(bw_croms_im_resized.size[0])

			##############################
			## CONVERSION TO FLAT ARRAY ##
			##############################

			arr = np.asarray(bw_croms_im_resized).flatten()
			arr = arr.astype('float32')

			im_matrix.append(arr)
			param_index.append(coef_label_dic[str(int(S))])
			iter_index.append([S,ITER])

		except:
			simulation_error.append(pointer)
			continue

	return(simulation_error,dim)

									##################################################################################################
									##################################################################################################
									############################################# MAIN ###############################################
									##################################################################################################
									##################################################################################################

#Get the type of input file
data_menu =True
while data_menu:
	#option_1, option_2, option_3 = options_menu("Raw Data Type", ['Simulation Files','Fasta Files', 'Simulation & Fasta Files'])
	data_type = int(argv[1]) #1
	if data_type==1:
		data_menu=False
	else:
		print("Wrong option selection. Enter your choice again")


if data_type == 1:
	#If simulation: Ask for folder or file name
	data_input = argv[2] #"NREPL=4/Simulations_1epoch"

	#These paramters are independent of whether the user is interested in generating images of a file or a folder
	simulations_results_outout = argv[3] #"NREPL=4/Images_1epoch/M_m/no_ordering"
	#Determine type of sorting to be carried out
	loop_2=True
	while loop_2:
		#option_1, option_2, option_3, option_4 = options_menu('Sorting Options', ['No Sorting','Sort Rows', 'Sort Columns', 'Sort Rows & Columns'])   ## Displays menu
		sort = int(argv[4]) #1
		if sort==1:
			loop_2=False
		elif sort==2:
			loop_2=False
		elif sort==3:
			loop_2=False
		elif sort==4:
			loop_2=False
		else:
			print("Wrong option selection. Enter your choice again")
	#Determine if threshold is to be applied and if so get the threshold
	loop_3=True
	while loop_3:
		#option_1, option_2 = options_menu("Apply Threshold", ['Yes','No'])
		app_thes = int(argv[5]) #2
		if app_thes==1:
			apply_threshold = True
			threshold  = get_threshold('threshold',0 ,1)
			print "Choice: A threshold of {} will be applied ".format(threshold)
			loop_3=False
		elif app_thes==2:
			apply_threshold = False
			threshold = 0.05
			loop_3=False
		else:
			print("Wrong option selection. Enter your choice again")
	#Major/Minor or Ancestral/Derived
	loop_4=True
	while loop_4:
		#option_1, option_2 = options_menu("Feature Selection", ['Major/Minor','Ancestral/Derived'])
		feature = int(argv[6]) #1
		if feature==1:
			maj_min = True
			loop_4=False
		elif feature==2:
			maj_min = False
			loop_4=False
		else:
			print("Wrong option selection. Enter your choice again")

	if os.path.isdir(data_input):
		dimensions = [] #What does this do?
		file_error = []
		col_error_parent = []

		with open(data_input + '/data.json', 'r') as fp:
			data = json.load(fp)
		#Count the number of files to be processed
		filecounter = 0
		for filename in os.listdir(data_input):
			filecounter += 1


		im_matrix = []
		param_index = []
		iter_index = []

		coef_label_dic = {}
		for n,l in enumerate(np.arange(0,410,10)):
			coef_label_dic[str(l)] = n

		for filename in os.listdir(data_input):	# For each simulation file:
			try:
				if filename.endswith('.txt') or filename.endswith('.gz'):
					if not filename.startswith('.'):
						string = '/%s' %(filename)
						specific_dict = filter(lambda x: x.get('name') == filename, data)
						S = float(specific_dict[0]['SAA'])/2
						NCHROMS = int(specific_dict[0]['NCHROM'])
						NREF = int(specific_dict[0]['NREF'])
						sim,dimen = image_simulation(data_input, simulations_results_outout,S, NREF, string, NCHROMS, threshold, apply_threshold, sort ,maj_min,128,128,coef_label_dic)

						if sim:
							d = {"simulation_file":filename, "column/s": sim}
							col_error_parent.append(d)
						dimensions.extend(dimen)

			except:
				file_error.append(string)
				continue


		#Save the errors
		if col_error_parent:
			if not os.path.exists(simulations_results_outout + '/errors'):
				os.makedirs(simulations_results_outout + '/errors')
			with open(simulations_results_outout + '/errors/simulation_error.json', 'w') as fp:
				json.dump(col_error_parent, fp, sort_keys=True, indent=4)



		####### Shuffle dataset #######

		im_matrix_shuffled, param_index_shuffled, iter_index = sklearn.utils.shuffle(im_matrix,param_index,iter_index,random_state=2)


		####### Split shuffled dataset #######

		data_split = [50,25,25]

		train_size = len(param_index_shuffled)/2
		val_size = len(param_index_shuffled)/4
		test_size = val_size

		train_im = im_matrix_shuffled[0:train_size]
		for i in train_im: i /= 255
		train_par = param_index_shuffled[0:train_size]

		val_im = im_matrix_shuffled[train_size:(train_size+val_size)]
		for i in val_im: i /= 255
		val_par = param_index_shuffled[train_size:(train_size+val_size)]

		test_im = im_matrix_shuffled[(train_size+val_size):]
		for i in test_im: i /= 255
		test_par = param_index_shuffled[(train_size+val_size):]


		train_im = np.asarray(train_im)
		train_im = train_im.astype(bool)
		f = gzip.GzipFile(simulations_results_outout + '/train_im.npy.gz', 'w')
		np.save(f, train_im, allow_pickle=False)
		f.close()
		np.save(simulations_results_outout + '/train_par', train_par, allow_pickle=False)



		val_im = np.asarray(val_im)
		val_im = val_im.astype(bool)
		f = gzip.GzipFile(simulations_results_outout + '/val_im.npy.gz', 'w')
		np.save(f, val_im, allow_pickle=False)
		f.close()
		np.save(simulations_results_outout + '/val_par', val_par, allow_pickle=False)



		test_im = np.asarray(test_im)
		test_im = test_im.astype(bool)
		f = gzip.GzipFile(simulations_results_outout + '/test_im.npy.gz', 'w')
		np.save(f, test_im, allow_pickle=False)
		f.close()
		np.save(simulations_results_outout + '/test_par', test_par, allow_pickle=False)


		iter_index = np.asarray(iter_index)
		np.save(simulations_results_outout + '/iter_index', iter_index, allow_pickle=False)
