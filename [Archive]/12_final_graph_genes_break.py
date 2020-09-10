# File name: 13_final_graph_genes_break.py
# Last Edit: 10/23/2019
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
# Using the GMM parameters derived from the program 6_GMM_adult.py,
# make a graph analogous to the one produced in 12_final_graph_genes_mRNA_break.py,
# but purely with protein comparisons. It is also similar to the graph
# produced in 9_final_graph.py, but indicating where some of the genes 
# of interest lie. The genes of interest are currently identified only 
# by comparing against MPPa, but this can be generalized with minor modifications to the code.

import sys
import numpy as np
import csv
import matplotlib.pyplot as plt

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
write_cell_dir = 'adult/distribution/'
data_dir = '3_gene_intensity_navg.txt'
data_dir_2 = 'adult/distribution/6_GMM_parameters.csv'

# Change this for analysis
cell_order = ['HSC', 'MPPa', 'MPPb', 'MPPc', 'CLP', 'CMP', 'MEP', 'GMP']
rainbow_color = {'HSC' : 'lightcoral', 'MPPa' : 'orange', 'MPPb' : 'forestgreen', 
'MPPc' : 'mediumturquoise', 'CLP' : 'dodgerblue', 'CMP' : 'darkblue', 'MEP' : 'violet', 'GMP' : 'crimson'}
MPPa_genes_reg = ['Hprt1']
MPPa_genes = ['Adnp', 'Ciapin1', 'Ssbp2']
MPPa_tick_pads = [0.0, 8.0, 0.0, 0.0] # Difference of 7 is sufficient

# Constant for ease
PI = np.pi
def flush():
	sys.stdout.flush()
STAT_COUNT = 1000 # Number of times to run a statistical model when finding the average

# Reads previously generated GMM parameters and returns the stats
# for the indicated cell type
def read_stats(cell, directory):
	with open(directory, 'rb') as r:
		reader = csv.reader(r, delimiter = ',')
		for row in reader:
			I_mu1 = row.index('mu1'); I_mu2 = row.index('mu2'); I_var1 = row.index('var1'); I_var2 = row.index('var2')
			I_weight1 = row.index('weight1'); I_weight2 = row.index('weight2'); I_UE = row.index('UE'); I_OE = row.index('OE')
			break;
		for row in reader:
			if row[0] == cell:
				return (float(row[I_mu1]), float(row[I_var1]), float(row[I_weight1]), 
					float(row[I_mu2]), float(row[I_var2]), float(row[I_weight2]), 
					float(row[I_UE]), float(row[I_OE]))

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

# Find the log 2 of fold change differences between mRNA and protein
# of cell_order[i] against combinations of all other cells 
def comb_fold(i, data):
	fold_changes_cont = [] # Tracks fold changes in an ordered, continuous list
	fold_changes_arr = [] # Tracks fold changes in an array
	highs_names = [] # Order of cells compared to ith cell
	highs_genes = [] # Ordered indices of gene names wrt 'genes' that are recorded in fold_changes_cont
	# Compare against all other cell types except for itself
	for j in range(len(cell_order)):
		if j != i:
	# Find mutual expressions and calculate the difference between log2 fold changes for protein and mRNA
			cell1_express_ind = (data[:, i] > 0).astype(int)
			cell2_express_ind = (data[:, j] > 0).astype(int)
			both_express_ind = (cell1_express_ind * cell2_express_ind).astype(bool)
			high_in_former = np.log2(data[both_express_ind, i] / data[both_express_ind, j])
	# Record the results
			fold_changes_cont = fold_changes_cont + high_in_former.tolist()
			fold_changes_arr.append(high_in_former)
			highs_names.append(cell_order[j])
			highs_genes.append(np.sort(np.where(both_express_ind))) # Sorting assumes that high_in_former is still in order
	fold_changes_dist = np.array(fold_changes_cont, dtype = float)

	return fold_changes_dist, fold_changes_arr, highs_names, highs_genes

# Directly modify the values and labels vector to include the gene label
# and the log protein expression value of the gene of interest
def get_gene_values(genes, highs_genes, highs_names, fold_changes_arr, OE, UE, values, labels):
	# Assumes that the gene is neither overexpressed nor underexpressed
	for i in range(len(MPPa_genes_reg)):
		gene = MPPa_genes_reg[i]
		my_list = genes[highs_genes[highs_names.index('MPPa')]].tolist()[0]
		if gene in my_list and gene:
			index = my_list.index(gene)
			values.append(fold_changes_arr[i][index])
			labels.append(gene)

	# Assumes that the gene should be overexpressed or underexpressed, and error-checks otherwise
	for i in range(len(fold_changes_arr)):
		if highs_names[i] == 'MPPa':
			high = fold_changes_arr[i] # Fold changes in proteins that are detected in both cell types
			gene_indices = highs_genes[i] # Indices of genes that are expressed in both cell types, order-matched to the data
			gene_names = genes[np.squeeze(gene_indices)].tolist()
			for j in range(len(MPPa_genes)):
				try:
					gene_index = gene_names.index(MPPa_genes[j])
					value = high[gene_index]
					if value > OE or value < UE:
						values.append(value)
						labels.append(MPPa_genes[j])
					else:
						print('Debugger: an indicated gene is not actually either UE/OE by value'); flush()
				except ValueError:
   					print('Debugger: an indicated gene cannot be found in UE/OE list'); flush()				

# Make a graph and return the figure and axes. Contains predetermined variables,
# unless otherwise stated.
def make_graph(title, x_axis, y_axis, width, height, 
	scale_width, scale_height, font_title, font_axes, x_ticks, y_ticks,
	pad_title, pad_x, pad_y, shift_x, shift_y):

	fig, (ax1, ax3, ax2) = plt.subplots(1, 3, sharey = False, figsize = (width ,height))
	fig.subplots_adjust(wspace = 0.05)

	ax1.set_title(title, size = font_title, pad = pad_title)
	ax1.set_xlabel(x_axis, size = font_axes, labelpad = pad_x)
	ax1.set_ylabel(y_axis, size = font_axes, labelpad = pad_y)

	plt.setp(ax1.get_xticklabels(), fontsize = x_ticks)
	plt.setp(ax2.get_xticklabels(), fontsize = x_ticks)
	plt.setp(ax3.get_xticklabels(), fontsize = x_ticks)
	plt.setp(ax1.get_yticklabels(), fontsize = y_ticks)

	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	ax3.tick_params(left = False, labelleft = False)

	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax2.tick_params(left = False, labelleft = False)

	box = ax1.get_position()
	ax1.set_position([box.x0 + shift_x * box.width, box.y0 + shift_y * box.width, box.width * scale_width, box.height * scale_height])
	box = ax2.get_position()
	ax2.set_position([box.x0 + shift_x * box.width - box.width * scale_width * 0.5, box.y0 + shift_y * box.width, box.width * scale_width, box.height * scale_height])
	box = ax3.get_position()
	ax3.set_position([box.x0 + shift_x * box.width, box.y0 + shift_y * box.width, box.width * scale_width * 0.5, box.height * scale_height])

	pos = ['data', 0]
	ax3.spines['left'].set_position(pos)

	return fig, ax1, ax2, ax3

# Plot the distribution given
def plot_distribution(i, fold_changes_dist, mu1, mu2, var1, var2, weight1, weight2, UE, OE, value, gene):
	# Vectorizable functions to use for plotting distributions
	def GMM_dist1(x):
		return weight1 * (np.exp(-((x - mu1)**2) / (2 * var1)) / (np.sqrt(2 * PI * var1))) 
	def GMM_dist2(x):
		return weight2 * (np.exp(-((x - mu2)**2) / (2 * var2)) / (np.sqrt(2 * PI * var2)))

	# Define parameters of legend
	x_shift = 0.68 # Increase -> 
	y_shift = 0.53 # Increase -> 
	font_size = 11.75
	tick_rotation = 0
	tick_length = 6.0
	gene_mark_offset = 0
	graph_side_offset = 0.030
	y_max_val = 0.01
	
	# Plot the histogram against split Gaussians & save the figure
	fig, ax1, ax2, ax3 = make_graph(title = 'Protein', 
		x_axis = '', y_axis = 'Normalized Count',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, x_ticks = 8.0, y_ticks = 12.0,
		pad_title = 5.0, pad_x = 6.0, pad_y = 0.0, shift_x = 0.50, shift_y = 0.10) # pad_value = 
	for ax in [ax1, ax2, ax3]:
		n, bins, patches = ax.hist(fold_changes_dist, density = 1, bins = 500, alpha = 1, color = 'silver')
		ax.plot(bins, np.vectorize(GMM_dist1)(bins), 'k', label = 'K = 2/1', color = 'green')
		ax.plot(bins, np.vectorize(GMM_dist2)(bins), 'k', label = 'K = 2/2', color = 'blue')
		ax.plot([UE, UE], [0, 0.3], 'r--', linewidth = 1.0, label = 'UE/OE', color = 'red') # UE cutoff location
		ax.plot([OE, OE], [0, 0.3], 'r--', linewidth = 1.0, color = 'red') # OE cutoff location
		x = np.array([value, value]).flatten() 
		y = np.concatenate((np.zeros(len(value), dtype = float), np.vectorize(GMM_dist2)(value))).flatten()
		for j in range(len(x) / 2):
			ax.plot(x[[j, j + len(x) / 2]] + gene_mark_offset, y[[j, j + len(x) / 2]], linewidth = 1.0, color = 'purple')

	counter = 0
	for ax in [ax1, ax2, ax3]: # ax3 is the middle
		if counter == 2:
			temp_value = np.array(value)
			temp_gene = np.array(gene)
			index = (temp_gene != 'Ssbp2').astype(int)
			index2 = (temp_gene != 'Fdx1').astype(int)
			fin_index = (index * index2).astype(bool)
			ax.xaxis.set_ticks(temp_value[index])
			ax.xaxis.set_ticklabels(temp_gene[index])
		else:
			ax.xaxis.set_ticks(value)
			ax.xaxis.set_ticklabels(gene)
		ax.tick_params(length = tick_length, rotation = tick_rotation)
		for ind, tick in enumerate(ax.xaxis.get_major_ticks()):
			tick.set_pad(MPPa_tick_pads[ind])
		counter += 1
	
	# fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)
	ax1.set_xlim([np.min(bins), UE + graph_side_offset])
	ax1.set_ylim([0, y_max_val])
	ax2.set_xlim([OE - graph_side_offset, np.max(bins)])
	ax2.set_ylim([0, y_max_val])
	ax3.set_xlim([UE - graph_side_offset, OE + 10 * graph_side_offset])
	ax3.set_ylim([0, .3])
	# plt.show()
	plt.savefig(analysis_dir + write_cell_dir + '13_' + cell_order[i] + '_genes_break')
	plt.close(fig)

# Main function
def main():
	# Starting message
	print('Starting [13_final_graph_genes_break]\n'); flush()

	# Reads previously generated GMM parameters and returns the stats for the indicated cell type
	mu1, var1, weight1, mu2, var2, weight2, UE, OE = read_stats('HSC', analysis_dir + data_dir_2)

	# Read the data generated in 3_data_compile.py
	genes, data = read_data_ordered(gen_data_dir + read_cell_dir + data_dir)
	print('Data reading complete...\n'); flush()

	# Find the log 2 of fold changes of cell_order[i] against combinations of all other cells 
	fold_changes_dist, fold_changes_arr, highs_names, highs_genes = comb_fold(cell_order.index('HSC'), data)

	# Get the gene label and the log protein expression value for genes I want to mark
	values = [] # Values of ticks on x-axis
	labels = [] # Labels of ticks on x-axis
	get_gene_values(genes, highs_genes, highs_names, fold_changes_arr, OE, UE, values, labels)

	# Plot the possible distributions against the data histogram
	plot_distribution(cell_order.index('HSC'), fold_changes_dist, mu1, mu2, var1, var2, weight1, weight2, UE, OE, values, labels)

	# Closing message
	print('\n[13_final_graph_genes_break] complete')

main()