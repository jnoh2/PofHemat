# File name: 8_mRNA_adult.py
# Last Edit: 10/09/2019
#
# Coded by: Joseph Noh
# Coded for: ZaroNoh et al. 2019
#
# PLEASE REFER TO THE READ_ME TO SEE HOW EACH .PY FILE RELATES TO
# THE PAPER
# 
# https://github.com/jnoh4/PofHemat
#
# Description: 
# Using the compiled data, compare to mRNA expressions and
# assume a K = 2 GMM to infer differentially translated genes
# 
# Also compare the distribution to a K = 1 distribution
# and report on the errors

import sys
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import linregress
from sklearn import mixture

# Constants used for directory names
#
# raw_data - file with raw data
# gen_data - file with generated data
# imp_data - file with data imported from uniprot.org
# analysis - file with analysis
# read_cell_dir - additional extension for file to read data from
# write_cell_dir - additional extension for file to write data to
# data_dir - directory of the data
raw_data_dir = './raw_data/'
gen_data_dir = './gen_data/'
imp_data_dir = './imp_data/'
analysis_dir = './analysis/'
read_cell_dir = 'adult_aged/'
write_cell_dir = 'adult/mRNA/'
data_dir = '3_gene_intensity_navg.txt'
mRNA_dir = './raw_data/adult/mRNA.csv'

# Change this for analysis
cell_order = ['HSC', 'MPPa', 'MPPb', 'MPPc']
rainbow_color = {'HSC' : 'lightcoral', 'MPPa' : 'orange', 'MPPb' : 'forestgreen', 
'MPPc' : 'mediumturquoise'}

# Constant for ease
PI = np.pi
def flush():
	sys.stdout.flush()
STAT_COUNT = 1000 # Number of times to run a statistical model when finding the average

# Reads the raw data from the given directory and returns the list of cells, genes, and the data
# Removes cells that are not used; removes genes that are not used
def read_data_ordered(directory):
	cells = [] # List of cells used
	genes = [] # List of genes
	data = [] # Intensity values
	with open(directory, 'rb') as r:
		reader = csv.reader(r, delimiter = ',')
	# Record all cell names in the file
		for row in reader:
			for i in range(len(row) - 2):
				cell = row[i + 2]
				if cell.find('Aged') != -1:
					cell = 'a' + cell[4:]
				cells.append(cell)
			break
	# Record all gene names in the file
		for row in reader:
			genes.append(row[1])
			data.append(row[2:])
	cells = np.array(cells)
	data = np.array(data, dtype = float)

	# Only include cells that are used
	data_copy = []
	for i in range(len(cell_order)):
		cell = cell_order[i]
		index = (cells == cell)
		data_copy.append(data[:, index])
	data = np.transpose(np.array(np.squeeze(data_copy), dtype = float))

	# Remove genes that are not expressed in any cell types
	data_TF = (data > 0).astype(int)
	data_TF = np.sum(data_TF, axis = 1)
	data_TF = (data_TF > 0).astype(bool)
	data = data[data_TF, :]
	genes = np.array(genes, dtype = str)[data_TF]

	return genes, data

# Retrieve mappings generated in 2_gene_prot_mappings.py
# prot_gene_map - UniprotID : gene name mapping
# gene_prot_map - gene name : UniprotID mapping
# prot_aliases - List of sets of UniprotID (k) to sets of gene name (k + 1) mappings where 0 <= k <= n / 2
def retrieve_mappings():
	prot_gene_map = {}; gene_prot_map = {}; prot_aliases = []

	# UniprotID : gene name mapping bidirectionally
	with open(gen_data_dir + '2_prot_gene_single.txt', 'rb') as f:
		reader = csv.reader(f, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
		for row in reader:
			prot_gene_map[row[0]] = row[1]
			gene_prot_map[row[1]] = row[0]

	# {UniprotIDs} : {gene names}
	with open(gen_data_dir + '2_prot_gene_group.txt', 'rb') as f:
		reader = csv.reader(f, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
		for row in reader:
	# Parse out row into {UniprotIDs} and {gene names}
			length = len(row)
			prots = set(('dummy',)) # 'dummy' used to initialize set of set
			genes = set(('dummy',))
			counter = 0 # Tracks column along row of a {UniprotIDs} : {gene names}
			while row[counter] != 'PROT:GENE': # Indicates end of list of UniprotIDs
				prots.add(row[counter])
				counter += 1
			counter += 1
			while counter < length: # Now at genes
				genes.add(row[counter])
				counter += 1
			prots.remove('dummy')
			genes.remove('dummy')
	# Add to prot_aliases
			prot_aliases.append(prots)
			prot_aliases.append(genes)

	return prot_gene_map, gene_prot_map, prot_aliases

# Reads the raw mRNA data from the given directory and returns an ordered, TPM mRNA
# dataset, removing cells that are not used; removes genes that are not used. Records 
# genes that are not used.
def read_mRNA_ordered(directory, genes, prot_gene_map, gene_prot_map, prot_aliases):
	mRNA_only_genes = [] # Genes expressed in mRNA, but not in proteins
	mRNA_unmapped = [] # Genes expressed in mRNA, but unable to be mapped to anything
	mRNA_genes = [] # List of genes expressed in mRNA AND in proteins
	mRNA_data = [] # Data from the mRNA file
	mRNA_cells = [] # Order of cells as given in the mRNA file; used only for temporary reasons

	with open(directory, 'rb') as r:
		reader = csv.reader(r, delimiter = ',')
		for row in reader:
			for i in range(len(row) - 1):
				mRNA_cells.append(row[i + 1])
			break
		for row in reader:
			gene = row[0] # mRNA gene name
			if gene in gene_prot_map: # Already mapped properly
				if gene not in genes: # Yes gene, no protein
					mRNA_only_genes.append(gene)
				else:
					mRNA_genes.append(gene)
					ROI = row[1:]
					for i in range(len(ROI)):
						if ROI[i] == '':
							ROI[i] = '0'
					mRNA_data.append(ROI)
			else: # Has not yet been mapped properly
				# Search in prot_aliases
				for i in range(len(prot_aliases) / 2):
					index = (2 * i) + 1
					if gene in prot_aliases[index]:
						for ele in prot_aliases[index]:
							if ele in gene_prot_map:
								gene = ele
								break
						break
				if gene not in gene_prot_map: # If gene is never found
					mRNA_unmapped.append(gene)
				else:
					if gene not in genes: # Yes gene, no protein
						mRNA_only_genes.append(gene)
					else:
						mRNA_genes.append(gene)
						ROI = row[1:]
						for i in range(len(ROI)):
							if ROI[i] == '':
								ROI[i] = '0'
						mRNA_data.append(ROI)

	# Order based on cells
	mRNA_data = np.array(np.squeeze(mRNA_data), dtype = float)
	mRNA_genes = np.array(mRNA_genes, dtype = str)
	new_data = []
	for i in range(len(cell_order)):
		cell = cell_order[i]
		index = mRNA_cells.index(cell)
		new_data.append(mRNA_data[:, index])
	mRNA_data = np.transpose(np.array(np.squeeze(new_data), dtype = float))

	# Record mRNA expressed, protein unexpressed
	w1 = open(analysis_dir + write_cell_dir + '8_mRNA_no_prot.csv', 'w+')
	for i in range(len(mRNA_only_genes)):
		w1.write(mRNA_only_genes[i] + '\n')
	w1.close()
	# Record mRNA expressed, no mappable protein
	w2 = open(analysis_dir + write_cell_dir + '8_mRNA_unmapped.csv', 'w+')
	for i in range(len(mRNA_unmapped)):
		w2.write(mRNA_unmapped[i] + '\n')
	w2.close()

	return mRNA_data, mRNA_genes

# Order the mRNA_data such that it is compatible with genes. Add 0 for any time the gene is not expressed
def expand_mRNA_genes(genes, mRNA_genes, mRNA_data):
	new_mRNA_data = []
	for i in range(len(genes)):
		gene = genes[i]
		if gene in mRNA_genes:
			index = mRNA_genes.tolist().index(gene)
			new_mRNA_data.append(mRNA_data[index,:])
		else:
			new_mRNA_data.append(np.zeros(len(mRNA_data[0,:]), dtype = float))

	return np.array(np.squeeze(new_mRNA_data), dtype = float)

# Find the log 2 of fold change differences between mRNA and protein
# of cell_order[i] against combinations of all other cells 
def comb_fold(i, data, mRNA_data, prot_mRNA_TF):
	fold_diff_cont = [] # Tracks fold changes in an ordered, continuous list
	fold_diff_arr = [] # Tracks fold changes in an array
	highs_names = [] # Order of cells compared to ith cell
	highs_genes = [] # Ordered indices of gene names wrt 'genes' that are recorded in fold_diff_cont
	# Compare against all other cell types except for itself
	for j in range(len(cell_order)):
		if j != i:
	# Find mutual expressions and calculate the difference between log2 fold changes for protein and mRNA
			both_express_ind = (prot_mRNA_TF[:, i] * prot_mRNA_TF[:, j]).astype(bool)
			prot_fold_change = np.log2(data[both_express_ind, i] / data[both_express_ind, j])
			mRNA_fold_change = np.log2(mRNA_data[both_express_ind, i] / mRNA_data[both_express_ind, j])
			fold_change_diff = prot_fold_change - mRNA_fold_change
	# Record the results
			fold_diff_cont = fold_diff_cont + fold_change_diff.tolist()
			fold_diff_arr.append(fold_change_diff)
			highs_names.append(cell_order[j])
			highs_genes.append(np.sort(np.where(both_express_ind))) # Sorting assumes that high_in_former is still in order
	fold_diff_dist = np.array(fold_diff_cont, dtype = float)

	return fold_diff_dist, fold_diff_arr, highs_names, highs_genes

# Given data, use the distribution to fit to a GMM num_times times and
# return the average parameters
def model_GMM(fold_diff_dist, num_times, num_distributions):
	mean1_avg = 0 # averaged mu of distribution 1
	mean2_avg = 0 # averaged mu of distribution 2
	var1_avg = 0 # averaged variance of distribution 1
	var2_avg = 0 # averaged variance of distribution 2
	weight1_avg = 0 # averaged weight of distribution 1
	weight2_avg = 0 # average weight of distribution 2
	
	for j in range(num_times):
		if j % 100 == 0:
			print('GMM Run: ' + str(j) + ' out of ' + str(num_times)); flush()
		g = mixture.GaussianMixture(n_components = num_distributions)
		g.fit(fold_diff_dist.reshape((-1, 1)))
		covariances = np.squeeze(g.covariances_)
		means = np.squeeze(g.means_)
		weights = np.squeeze(g.weights_)
	# Combine parameters of distributions into tuples and sort by covariance
	# Assumes lower variance is distribution 1
		stats = []
		for k in range(len(means)):
			stats.append((means[k], covariances[k], weights[k]))
		stats.sort(key = lambda x: x[1]) # Ascending
		mean1_avg += stats[0][0] / float(num_times)
		mean2_avg += stats[1][0] / float(num_times)
		var1_avg += stats[0][1] / float(num_times)
		var2_avg += stats[1][1] / float(num_times)
		weight1_avg += stats[0][2] / float(num_times)
		weight2_avg += stats[1][2] / float(num_times)
	# Correct for not summing to 1.0
	weight2_avg = 1.0 - weight1_avg

	return mean1_avg, mean2_avg, var1_avg, var2_avg, weight1_avg, weight2_avg

# Determines the overexpression & underexpression cut-offs based on the GMM parameters
def cutoff(mean1_avg, mean2_avg, var1_avg, var2_avg):
	# Knowing that distribution 2 is the biologically interesting distribution, determine possible overexpression & underexpression cut-offs
	OE1 = mean1_avg + 3 * np.sqrt(var1_avg)
	UE1 = mean1_avg - 3 * np.sqrt(var1_avg)
	OE2 = mean2_avg + 2 * np.sqrt(var2_avg)
	UE2 = mean2_avg - 2 * np.sqrt(var2_avg)
	# Select the more stringent cut-offs
	if OE1 > OE2:
		OE_select = OE1
	else: 
		OE_select = OE2
	if UE1 < UE2:
		UE_select = UE1
	else:
		UE_select = UE2

	return OE_select, UE_select

# Make a graph and return the figure and axes. Contains predetermined variables,
# unless otherwise stated.
def make_graph(title, x_axis, y_axis, width, height, 
	scale_width, scale_height, font_title, font_axes, font_ticks, 
	pad_title, pad_x, pad_y, shift_x, shift_y):

	fig, ax = plt.subplots(figsize = (width ,height))

	ax.set_title(title, size = font_title, pad = pad_title)
	ax.set_xlabel(x_axis, size = font_axes, labelpad = pad_x)
	ax.set_ylabel(y_axis, size = font_axes, labelpad = pad_y)

	plt.setp(ax.get_xticklabels(), fontsize = font_ticks)
	plt.setp(ax.get_yticklabels(), fontsize = font_ticks)

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	box = ax.get_position()
	ax.set_position([box.x0 + shift_x * box.width, box.y0 + shift_y * box.width, box.width * scale_width, box.height * scale_height])

	return fig, ax

# Plot the distribution given
def plot_distribution(i, fold_diff_dist, mean1_avg, mean2_avg, single_mean, var1_avg, var2_avg, single_std, weight1_avg, weight2_avg):
	# Vectorizable functions to use for plotting distributions
	def GMM_dist(x):
		return (weight1_avg * (np.exp(-((x - mean1_avg)**2) / (2 * var1_avg)) / (np.sqrt(2 * PI * var1_avg))) 
			+ weight2_avg * (np.exp(-((x - mean2_avg)**2) / (2 * var2_avg)) / (np.sqrt(2 * PI * var2_avg))))
	def single_dist(x):
		return (np.exp(-((x - single_mean)**2) / (2 * single_std * single_std)) / (np.sqrt(2 * PI * single_std * single_std))) 
	def GMM_dist1(x):
		return weight1_avg * (np.exp(-((x - mean1_avg)**2) / (2 * var1_avg)) / (np.sqrt(2 * PI * var1_avg))) 
	def GMM_dist2(x):
		return weight2_avg * (np.exp(-((x - mean2_avg)**2) / (2 * var2_avg)) / (np.sqrt(2 * PI * var2_avg)))

	# Define parameters of legend
	x_shift = 0.68 # Increase -> 
	y_shift = 0.53 # Increase -> 
	font_size = 11.75
	
	# Plot the histogram against GMM & save the figure
	fig, ax = make_graph(title = cell_order[i] + ' Protein/mRNA', 
		x_axis = 'Log2 Fold Change (Prot - mRNA)', y_axis = 'Normalized Count',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, font_ticks = 12.0, 
		pad_title = 5.0, pad_x = 6.0, pad_y = -2.0, shift_x = 0.08, shift_y = 0.10) # pad_value = 
	n, bins, patches = ax.hist(fold_diff_dist, density = 1, bins = 500, alpha = 1, color = 'silver')
	ax.plot(bins, np.vectorize(GMM_dist)(bins), 'k', label = 'K = 2', color = rainbow_color[cell_order[i]])
	fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)
	plt.savefig(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_K2')
	plt.close(fig)

	# Plot the histogram against single Gaussian model & save the figure
	fig, ax = make_graph(title = cell_order[i] + ' Protein/mRNA', 
		x_axis = 'Log2 Fold Change (Prot - mRNA)', y_axis = 'Normalized Count',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, font_ticks = 12.0, 
		pad_title = 5.0, pad_x = 6.0, pad_y = -2.0, shift_x = 0.08, shift_y = 0.10) # pad_value = 
	n, bins, patches = ax.hist(fold_diff_dist, density = 1, bins = 500, alpha = 1, color = 'silver')
	ax.plot(bins, np.vectorize(single_dist)(bins), 'k', label = 'K = 1', color = 'teal')
	fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)
	plt.savefig(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_K1')
	plt.close(fig)
	
	# Plot the histogram against split Gaussians & save the figure
	fig, ax = make_graph(title = cell_order[i] + ' Protein/mRNA', 
		x_axis = 'Log2 Fold Change (Prot - mRNA)', y_axis = 'Normalized Count',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, font_ticks = 12.0, 
		pad_title = 5.0, pad_x = 6.0, pad_y = -2.0, shift_x = 0.08, shift_y = 0.10) # pad_value = 
	n, bins, patches = ax.hist(fold_diff_dist, density = 1, bins = 500, alpha = 1, color = 'silver')
	ax.plot(bins, np.vectorize(GMM_dist1)(bins), 'k', label = 'K = 2/1', color = 'green')
	ax.plot(bins, np.vectorize(GMM_dist2)(bins), 'k', label = 'K = 2/2', color = 'blue')
	fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)
	plt.savefig(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_K2_Split')
	plt.close(fig)

# Given the parameters of the distributions, generate data points
def generate_distributions(mean1_avg, var1_avg, n1, mean2_avg, var2_avg, n2, single_mean, single_std, length):
	dist1 = np.random.normal(mean1_avg, np.sqrt(var1_avg), n1) # Distribution 1
	dist2 = np.random.normal(mean2_avg, np.sqrt(var2_avg), n2) # Distribution 2
	dist_comb = np.concatenate((dist1, dist2)) # Distribution 1 + 2
	dist_single = np.random.normal(single_mean, single_std, length) # Single Gaussian Distribution
	dist_single_2 = np.random.normal(single_mean, single_std, length)

	# For debugging purposes
	if length != len(dist_comb):
		print('Debugger: length data != length of GMM generated distribution. Needs debugging'); flush()
	else:
		print('Debugger: length data == length of GMM generated distribution. No problems'); flush()

	return dist_single, dist_single_2, dist_comb

# Given the OE and UE cutoff values, record the proteins that are relatively over- or under- expressed
def record_OE_UE(i, OE, UE, fold_diff_arr, highs_names, highs_genes, genes):
	w1 = open(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_overexpressions.txt', 'w+')
	w2 = open(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_underexpressions.txt', 'w+')
	w1.write('Overexpress is above ' + str(OE) + ' and underexpress is below ' + str(UE) + '\n')
	w2.write('Overexpress is above ' + str(OE) + ' and underexpress is below ' + str(UE) + '\n')

	# For cell type cell_order[i], record against cell type of high_names[k]
	for k in range(len(fold_diff_arr)):
		cell = highs_names[k] # Cell type to compare against
		high = fold_diff_arr[k] # Fold changes in proteins that are detected in both cell types
		gene_indices = highs_genes[k] # Indices of genes that are expressed in both cell types, order-matched to the data
		high_TF = high > OE # Fold changes making the overexpression cutoff
		low_TF = high < UE # Fold changes making the underexpression cutoff

		w1.write('Compared against ' + highs_names[k] + '(%i): ' %(np.sum(high_TF.astype(int))) + '\n')
		w1.write('. '.join(genes[np.squeeze(gene_indices)[np.squeeze(high_TF)]]) + '\n\n')

		w2.write('Compared against ' + highs_names[k] + '(%i): ' %(np.sum(low_TF.astype(int))) + '\n')
		w2.write('. '.join(genes[np.squeeze(gene_indices)[np.squeeze(low_TF)]]) + '\n\n')

	w1.close()
	w2.close()

# Given the possible distributions and the data, make the QQ plots
def plot_QQ(i, dist_single, dist_single_2, dist_comb, fold_diff_dist):
	fig, ax = make_graph(title = cell_order[i] + ' QQ', 
		x_axis = 'Log2 Fold Diff', y_axis = 'Log2 Fold Diff',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, font_ticks = 12.0, 
		pad_title = 5.0, pad_x = 4.0, pad_y = 1.0, shift_x = 0.15, shift_y = 0.10) 
	
	# Plot the null distribution
	ax.plot(np.sort(dist_single), np.sort(dist_single_2), 'k', label = 'null', color = 'black', alpha = 0.9)
	# Plot the single Gaussian against the actual distribution
	ax.plot(np.sort(fold_diff_dist), np.sort(dist_single), 'k', label = 'K = 1', color = 'teal', alpha = 0.9)
	# Plot the GMM model against the actual distribution
	ax.plot(np.sort(fold_diff_dist), np.sort(dist_comb), 'k', label = 'K = 2', color = rainbow_color[cell_order[i]], alpha = 0.9)
	
	# Define parameters of legend
	x_shift = 0.68 # Increase -> 
	y_shift = 0.33 # Increase -> 
	font_size = 11.75

	fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)
	plt.savefig(analysis_dir + write_cell_dir + '8_' + cell_order[i] + '_QQPlot')
	plt.close(fig)

# Given the parameters of all the distributions, find the average errors of linearly regressing QQ plots
def get_linreg_errors(mean1_avg, var1_avg, n1, mean2_avg, var2_avg, n2, single_mean, single_std, fold_diff_dist):
	single_errors = []
	GMM_errors = []
	null_errors = []

	# Run regression 1,000 times to determine the average error
	for l in range(STAT_COUNT):
		if l % 500 == 0:
			print('Linreg Run: ' + str(l) + ' out of ' + str(STAT_COUNT)); flush()
		dist1 = np.random.normal(mean1_avg, np.sqrt(var1_avg), n1)
		dist2 = np.random.normal(mean2_avg, np.sqrt(var2_avg), n2)
		dist_comb = np.concatenate((dist1, dist2))
		dist_single = np.random.normal(single_mean, single_std, n1 + n2)
		dist_single_2 = np.random.normal(single_mean, single_std, n1 + n2)
		s1, i1, r1, p1, e1 = linregress(np.sort(fold_diff_dist), np.sort(dist_single))
		s2, i2, r2, p2, e2 = linregress(np.sort(fold_diff_dist), np.sort(dist_comb))
		s3, i3, r3, p3, e3 = linregress(np.sort(dist_single), np.sort(dist_single_2))
		single_errors.append(e1)
		GMM_errors.append(e2)
		null_errors.append(e3)

	return single_errors, GMM_errors, null_errors

# Given boxplot_gatherings of error values of QQ plots, make boxplot & save the file
def plot_errors(boxplot_gathering):
	fig, ax = make_graph(title = 'QQ Plot Squared Error', 
		x_axis = 'QQ Plot', y_axis = 'Squared Error',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, font_ticks = 12.0, 
		pad_title = 5.0, pad_x = 4.0, pad_y = 3.0, shift_x = 0.19, shift_y = 0.10) 

	# Make boxplots & make their color appropriate
	for i in range(len(cell_order)):
		DOI = list(boxplot_gathering[i])
		bp = ax.boxplot(DOI, sym = '', labels = ['K = 1', 'K = 2', 'null'], patch_artist = True, showbox = False)
		plt.setp(bp['medians'][0], color = rainbow_color[cell_order[i]])
		plt.setp(bp['medians'][1], color = rainbow_color[cell_order[i]])
		plt.setp(bp['medians'][2], color = rainbow_color[cell_order[i]])

	# Make the legend for the boxplot
	leg = []
	for i in range(len(cell_order)):
		leg.append(mpatches.Patch(color = rainbow_color[cell_order[i]], label = cell_order[i]))

	# Define parameters of legend
	x_shift = 0.76 # Increase -> 
	y_shift = 0.57 # Increase -> 
	font_size = 11.75
	color_length = 3.0
	fig.legend(handles = leg, loc = 'center left', bbox_to_anchor = (x_shift, y_shift), 
		frameon = False, handlelength = 0.3, fontsize = font_size)

	plt.savefig(analysis_dir + write_cell_dir + '8_sqrErrors.png')
	plt.close(fig)

# Main function
def main():
	# Starting message
	print('Starting [8_mRNA_adult]\n'); flush()

	# Read the data generated in 3_data_compile.py
	genes, data = read_data_ordered(gen_data_dir + read_cell_dir + data_dir)
	print('Regular data reading complete...\n'); flush()
	# Retrieve mappings generated in 2_gene_prot_mappings.py
	prot_gene_map, gene_prot_map, prot_aliases = retrieve_mappings()
	# Read the raw TPM mRNA data & return as TPM
	mRNA_data, mRNA_genes = read_mRNA_ordered(mRNA_dir, genes, prot_gene_map, gene_prot_map, prot_aliases)
	# Order the mRNA_data such that it is compatible with genes. Add 0 for any time the gene is not expressed
	mRNA_data = expand_mRNA_genes(genes, mRNA_genes, mRNA_data)
	print('mRNA data reading complete...\n'); flush()

	# Determine locations where both prot and mRNA are found
	data_TF = (data > 0).astype(int)
	mRNA_data_TF = (mRNA_data > 0).astype(int)
	prot_mRNA_TF = data_TF * mRNA_data_TF

	# Open files to record R_squared values of QQ plot and
	# statistics involved with the GMM
	w1 = open(analysis_dir + write_cell_dir + '8_QQ_R_squared.csv', 'w+')
	w2 = open(analysis_dir + write_cell_dir + '8_GMM_parameters.csv', 'w+')
	w1.write(',K = 1,K = 2,null\n')
	w2.write(',mu1,var1,weight1,mu2,var2,weight2,UE,OE\n')

	# List of tuples containing R-squared values for each cell-type comparison
	# (K = 1, K = 2, null distribution)
	boxplot_gathering = []

	# Go through all cell types present
	for i in range(len(cell_order)):
		print(cell_order[i], 'beginning...'); flush()
	# Find the log 2 of fold change differences between mRNA and protein of cell_order[i] against combinations of all other cells 
		fold_diff_dist, fold_diff_arr, highs_names, highs_genes = comb_fold(i, data, mRNA_data, prot_mRNA_TF)
		
	# Record statistics
		length = len(fold_diff_dist); single_mean = np.mean(fold_diff_dist); single_std = np.std(fold_diff_dist)

	# Find parameters of GMM-derived distributions
		print(cell_order[i], 'GMM...'); flush()
		mean1_avg, mean2_avg, var1_avg, var2_avg, weight1_avg, weight2_avg = model_GMM(fold_diff_dist, STAT_COUNT, 2)
		n1 = int(math.ceil(length * weight1_avg)) # Number of data points from distribution 1
		n2 = int(length - n1) # Number of data points from distribution 2
	
	# Make sure var2_avg > var1_avg
		if var2_avg > var1_avg:
			print('Debugger: var2 > var1; no issues'); flush()
		else: 
			print('Debugger: var1 > var2; code assumption breaks down -- get checked!'); flush()
	# Determine cutoff values
		OE, UE = cutoff(mean1_avg, mean2_avg, var1_avg, var2_avg)

	# Write the parameters and cutoffs into file
		w2.write(cell_order[i] + ',' + str(mean1_avg) + ',' + str(var1_avg) + ',' + str(weight1_avg)
			+ ',' + str(mean2_avg) + ',' + str(var2_avg) + ',' + str(weight2_avg) + ',' + str(UE) + ',' + str(OE)+ '\n')

	# Plot the possible distributions against the data histogram
		plot_distribution(i, fold_diff_dist, mean1_avg, mean2_avg, single_mean, var1_avg, var2_avg, single_std, weight1_avg, weight2_avg)
	
	# Generate distributions for the GMM & the single Gaussian model
		dist_single, dist_single_2, dist_comb = generate_distributions(mean1_avg, var1_avg, n1, mean2_avg, var2_avg, n2, single_mean, single_std, length)
		
	# Record the proteins that make the cutoff values of overexpression and underexpression
		print(cell_order[i], 'OE/UE...'); flush()
		record_OE_UE(i, OE, UE, fold_diff_arr, highs_names, highs_genes, genes)
	
	# Make a QQ plot of the null distribution (Gaussian vs Gaussian), K = 1 (Gaussian vs data) and K = 2 (GMM vs data)
		plot_QQ(i, dist_single, dist_single_2, dist_comb, fold_diff_dist)
		
	# Write the means of the errors of linear regression for QQ plots
		print(cell_order[i], 'linreg errors...'); flush()
		single_errors, GMM_errors, null_errors = get_linreg_errors(mean1_avg, var1_avg, n1, mean2_avg, var2_avg, n2, single_mean, single_std, fold_diff_dist)
		w1.write(cell_order[i] + ',' + str(np.mean(single_errors)) + ',' + str(np.mean(GMM_errors)) + ',' + str(np.mean(null_errors)) + '\n')
	# Record the error data for this specific data type
		boxplot_gathering.append((single_errors, GMM_errors, null_errors))

		print(cell_order[i], 'complete!'); print('\n'); flush()
		
	# Plot the errors of the QQ plots for all cell types
	plot_errors(boxplot_gathering)
	
	# Closing message
	print('\n[8_mRNA_adult] complete')

main()