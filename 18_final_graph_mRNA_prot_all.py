# File name: 18_final_graph_mRNA_prot_all.py
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
# Generate a boxplot of the log2 fold change values between protein
# and mRNA per cell type, only using the genes that are present in
# protein and mRNA for all cell types. Indicate where the specific 
# genes of interest (Hprt1, Ciapin1, Adnp) lie along this boxplot.

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
# dataset, removing cells that are not used; removes genes that are not used
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

# Process the data such that only genes expressed across protein and mRNA
# of all cell types are considered. Normalize within this geneset
def genes_all_norm(data, mRNA_data, genes):
	# Expressed across protein and mRNA of all cell types
	all_express = (np.prod((mRNA_data > 0).astype(int), axis = 1) * np.prod((data > 0).astype(int), axis = 1)).astype(bool)

	# Normalize protein to 1,000,000 
	data = data[all_express, :]
	data = 1000000 * data / np.sum(data, axis = 0)

	# Normalize mRNA to 1,000,000
	mRNA_data = mRNA_data[all_express, :]
	mRNA_data = 1000000 * mRNA_data / np.sum(mRNA_data, axis = 0)

	# Include only those specific genes
	genes = genes[all_express]

	return data, mRNA_data, genes

# Make a graph and return the figure and axes. Contains predetermined variables,
# unless otherwise stated.
def make_graph(title, x_axis, y_axis, width, height, 
	scale_width, scale_height, font_title, font_axes, x_ticks, y_ticks,
	pad_title, pad_x, pad_y, shift_x, shift_y):

	fig, ax = plt.subplots(figsize = (width ,height))

	ax.set_title(title, size = font_title, pad = pad_title)
	ax.set_xlabel(x_axis, size = font_axes, labelpad = pad_x)
	ax.set_ylabel(y_axis, size = font_axes, labelpad = pad_y)

	plt.setp(ax.get_xticklabels(), fontsize = x_ticks)
	plt.setp(ax.get_yticklabels(), fontsize = y_ticks)

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	box = ax.get_position()
	ax.set_position([box.x0 + shift_x * box.width, box.y0 + shift_y * box.width, box.width * scale_width, box.height * scale_height])

	return fig, ax

# Make a boxplot of the log2 of protein over mRNA fold change values per cell type,
# indicating where the specific genes of interest lie
def plot_data(hprt1_data, adnp_data, ciapin1_data, bp_gatherings):
	fig, ax = make_graph(title = 'Simple prot/mRNA', 
			x_axis = 'Cell Type', y_axis = 'Log2 Fold Change Difference',
			width = 4.0, height = 3.2, scale_width = 0.72, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
			font_title = 14.2, font_axes = 12.5, x_ticks = 12.0, y_ticks = 12.0,
			pad_title = 5.0, pad_x = 6.0, pad_y = -2.0, shift_x = 0.08, shift_y = 0.10)

	# Plot the genes of interest
	ax.plot(np.arange(1, 5, 1), hprt1_data, ':', c = '#808080')
	ax.plot(np.arange(1, 5, 1), hprt1_data, 'o', c = '#000000', label = 'Hprt1', markerfacecolor = '#ffffff')
	ax.plot(np.arange(1, 5, 1), adnp_data, ':', c = '#404040')
	ax.plot(np.arange(1, 5, 1), adnp_data, 'o', c = '#000000', label = 'Adnp', markerfacecolor = '#808080')
	ax.plot(np.arange(1, 5, 1), ciapin1_data, ':', c = '#000000')
	ax.plot(np.arange(1, 5, 1), ciapin1_data, 'o', c = '#000000', label = 'Ciapin1', markerfacecolor = '#000000')

	# Define parameters of legend
	x_shift = 0.68 # Increase -> 
	y_shift = 0.53 # Increase -> 
	font_size = 11.75

	# Make the legend only for the genes of interest
	fig.legend(loc = 'center left', bbox_to_anchor = (x_shift, y_shift), frameon = False, fontsize = font_size)

	# Plot the larger data of log2 of protein over mRNA fold change
	bp = ax.boxplot(bp_gatherings, labels = cell_order, showmeans = True, meanline = True, showfliers = False, whis = 0, showcaps = False)

	# Color code the graph 
	for i, patch in enumerate(bp['boxes']):
		patch.set(color = rainbow_color[cell_order[i]])
	for i, patch in enumerate(bp['medians']):
		patch.set(color = rainbow_color[cell_order[i]])
	for i, patch in enumerate(bp['means']):
		patch.set(color = rainbow_color[cell_order[i]])
	for i, patch in enumerate(bp['whiskers']):
		patch.set(color = rainbow_color[cell_order[i / 2]])
	for i, patch in enumerate(bp['caps']):
		patch.set(color = rainbow_color[cell_order[i / 2]])

	# Save the figure
	plt.savefig(analysis_dir + write_cell_dir + '18_mRNA_prot_all')

# Main function
def main():
	# Starting message
	print('Starting [18_final_graph_mRNA_prot_all]\n'); flush()

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

	# Process the data such that only genes expressed across protein and mRNA
	# of all cell types are considered. Normalize within this geneset
	data, mRNA_data, genes = genes_all_norm(data, mRNA_data, genes)

	# Get the expresion log fold change values for the genes of interest
	index_ciapin1 = genes.tolist().index('Ciapin1')
	index_adnp = genes.tolist().index('Adnp')
	index_hprt1 = genes.tolist().index('Hprt1')

	ciapin1_data = np.log2(data[index_ciapin1, :] / mRNA_data[index_ciapin1, :])
	adnp_data = np.log2(data[index_adnp, :] / mRNA_data[index_adnp, :])
	hprt1_data = np.log2(data[index_hprt1, :] / mRNA_data[index_hprt1, :])

	# Determine the genes expressed in both protein and mRNA per cell type; shoud be redundant here
	both_expressed = ((data > 0).astype(int) * (mRNA_data > 0).astype(int)).astype(bool)

	# Find the log fold change for each gene expressed in protein and mRNA for each cell type
	bp_gatherings = []
	for i in range(len(cell_order)):
		log_fold_data = np.log2(data[both_expressed[:, i], i] / mRNA_data[both_expressed[:, i], i])
		bp_gatherings.append(log_fold_data)

	# Plot the gene-specific data, as well as the overall fold change values
	plot_data(hprt1_data, adnp_data, ciapin1_data, bp_gatherings)
	
	# Closing message
	print('\n[18_final_graph_mRNA_prot_all] complete')

main()