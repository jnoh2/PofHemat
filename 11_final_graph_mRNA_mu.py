# File name: 11_final_graph_adult_mu.py
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
# Using the GMM parameters derived from the program 8_mRNA_adult.py,
# create a bar graph of mu2's compared between cell types. Higher
# values of mu2 indicates more protein per mRNA by fold change overall.
# As a default, the option to align mu1's is turned off

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
data_dir = 'adult/mRNA/8_GMM_parameters.csv'

# Change this for analysis
cell_order = ['HSC', 'MPPa', 'MPPb', 'MPPc']
rainbow_color = {'HSC' : 'lightcoral', 'MPPa' : 'orange', 'MPPb' : 'forestgreen', 'MPPc' : 'mediumturquoise'}

# Constant for ease
PI = np.pi
def flush():
	sys.stdout.flush()

# Reads previously generated GMM parameters and returns the stats
def read_stats(directory):
	all_data = []
	with open(directory, 'rb') as r:
		reader = csv.reader(r, delimiter = ',')
		for row in reader:
			I_mu1 = row.index('mu1'); I_mu2 = row.index('mu2'); I_var1 = row.index('var1'); I_var2 = row.index('var2')
			I_weight1 = row.index('weight1'); I_weight2 = row.index('weight2'); I_UE = row.index('UE'); I_OE = row.index('OE')
			break;
		for row in reader:
			all_data.append([row[I_mu1], row[I_var1], row[I_weight1], row[I_mu2], row[I_var2], row[I_weight2], row[I_UE], row[I_OE]])
	
	return np.array(all_data, dtype = float)

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

# Main function
def main():
	# Starting message
	print('Starting [11_final_graph_mRNA_mu]\n'); flush()

	# Reads previously generated GMM parameters and returns the stats
	# mu1, var1, weight1, mu2, var2, weight2, UE, OE
	data = read_stats(analysis_dir + data_dir)

	# Leave the option to adjust the data such that mu1's align
	# data[:, [0, 3]] = data[:,[0, 3]] - np.vstack((data[:, 0], data[:, 0])).T

	# Make a bar graph of mu2's
	color_list = []
	# Make list of colors for bars in bar graph
	for i in range(len(cell_order)):
		color_list.append(rainbow_color[cell_order[i]])

	fig, ax = make_graph(title = 'Log2 Fold Difference Means', 
		x_axis = 'Cell Type', y_axis = 'Distribution 2 Mean',
		width = 4.0, height = 3.2, scale_width = 0.9, scale_height = 0.9, # scale_value represents size of graph wrt entire screen
		font_title = 14.2, font_axes = 12.5, x_ticks = 12.0, y_ticks = 12.0,
		pad_title = 5.0, pad_x = 6.0, pad_y = -2.0, shift_x = 0.08, shift_y = 0.10) # pad_value = 
	ax.bar(cell_order, data[:, 3], color = color_list)
	plt.savefig(analysis_dir + 'adult/mRNA/11_mRNA_mu')
	plt.close(fig)

	# Closing message
	print('\n[11_final_graph_mRNA_mu] complete')

main()