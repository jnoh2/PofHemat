# File name: 7_unique_adult.py
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
# Using the compiled data, find the unique proteins expressions for
# different combinations of cells and record them

import sys
import numpy as np
import csv

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
write_cell_dir = 'adult/uniques/'
data_dir = '3_gene_intensity_navg.txt'

# Constant for debugging
num_comparisons = 0
def flush():
	sys.stdout.flush()

#Change this for analysis
cell_order = ['HSC', 'MPPa', 'MPPb', 'MPPc', 'CLP', 'CMP', 'MEP', 'GMP']

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
	
# Recusive function to be used to recurse through gene list.
# remove_list - tracks the indices of cells that must be considered as 'not expressing'
# next_cell - indicates the cell that is to be removed next
# count - count of cells that must be removed more until expressed genes are computed
# data_TF - 
def count_expression(remove_list, next_cell, count, data_TF, genes, w1):
	global num_comparisons
	# count == 0 implies cells have been removed as prescribed at the beginning of this recursive loop
	if count == 0:
		num_comparisons += 1
	# Make the keep list
		keep_list = []
		for i in range(len(cell_order)):
			if i not in remove_list:
				keep_list.append(i)
	# Genes expressed in keep list, not expressed in remove list
		data_TF_keep = data_TF[:, keep_list]
		data_TF_remove = (data_TF[:, remove_list] * (-1) + 1).astype(bool)
	# Find the combination & apply to gene list
		data_TF_comb = np.prod(data_TF_keep, axis = 1) * np.prod(data_TF_remove, axis = 1)
		genes_expressed = genes[data_TF_comb.astype(bool)]
	# List of cells kept
		kept_cells = np.array(cell_order)[keep_list]
	# Write this into file
		w2 = open(analysis_dir + write_cell_dir + '7_unique_' + '_'.join(kept_cells) + '.txt', 'w+')
		for i in range(len(genes_expressed)):
			w2.write(genes_expressed[i] + '\n')
		w2.close()
	# Keep record of genes expressed with different combinations of cells
		if len(genes_expressed) > 0:
			w1.write('_'.join(kept_cells) + ',' + str(len(genes_expressed)) + '\n')
		return(len(genes_expressed))
	# BASE CASE
	elif next_cell >= len(cell_order):
		return 0 
	else:
	# Recurse until final cases are reached
		genes_covered = 0
		for j in range(len(cell_order)):
	# Remove from current cell (j = 0) till the end
			next_cell_rem = next_cell + j
			if next_cell_rem < len(cell_order):
				next_remove_list = remove_list[:]
				next_remove_list.append(next_cell_rem)
				genes_covered += count_expression(next_remove_list, next_cell_rem + 1, count - 1, data_TF, genes, w1)
		return(genes_covered)

# Function purely for debugging; makes sure everything ran OK
def debug(genes, genes_covered, num_comparisons):
	if len(genes) == genes_covered:
		print('Debugger: Final count of genes covered is good'); flush()
	else:
		print('Debugger: Final count of genes covered is NOT GOOD. Check!'); flush()

	if num_comparisons != (2 ** len(cell_order)) - 1: # Sum of all combinations minus 1 (removing all cells)
		print('Debugger: Number of comparisons is NOT good'); flush()
		print('Debugger: comparisons - ' + str(num_comparisons) + ' theoretical - ' + str(2 ** len(cell_order))); flush()
	else:
		print('Debugger: Number of comparisons matches the theoretical count'); flush()

# Main function
def main():
	global num_comparisons
	# Starting message
	print('Starting [7_unique_adult.py]\n'); flush()

	# Read in data
	genes, data = read_data_ordered(gen_data_dir + read_cell_dir + data_dir)
	data_TF = (data > 0)

	w1 = open(analysis_dir + write_cell_dir + '7_unique_counts_adult.csv', 'w+')
	genes_covered = 0

	# For each cell type
	for i in range(len(cell_order)):
	# i indicates the number of cells that don't express a certain list of genes
		if i > 0:
			for j in range(len(cell_order)):
				remove_list = []
				remove_list.append(j)
				genes_covered += count_expression(remove_list, j + 1, i - 1, data_TF, genes, w1)
	# Case when no cells are removed
		else:
			num_comparisons += 1
			data_TF_last = np.prod(data_TF, axis = 1)
			genes_expressed = genes[data_TF_last.astype(bool)]
			genes_covered += len(genes_expressed)
			w1.write('_'.join(cell_order) + ',' + str(len(genes_expressed)) + '\n')
			w2 = open(analysis_dir + write_cell_dir + '7_unique_' + '_'.join(cell_order) + '.txt', 'w+')
			for i in range(len(genes_expressed)):
				w2.write(genes_expressed[i] + '\n')
			w2.close()
	w1.close()

	# Debugging purposes
	debug(genes, genes_covered, num_comparisons)s	

	# Closing message
	print('\n[7_unique_adult] complete')

main()