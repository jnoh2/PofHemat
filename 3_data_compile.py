# File name: 3_data_compile.py
# Last Edit: 10/05/2019
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
# Using the mappings created in 2_gene_prot_mapping.py, compile data from all cell types to create
#
# 1) compiled raw data
# 2) compiled data normalized to 1,000,000
# 3) compiled data averaged across equivalent cell types, thn normalized to 1,000,000
# 4) compiled data such that the geometric mean within cell type, ignoring 0's, is normalized to 1,000,000
# 5) compiled data such that the geometric mean within cell type, including 0's, is normalized to 1,000,000, and
# 6) compiled data averaged across the equivalent cell types, with each cell normalized to 1,000,000 first.
# 
# written into six different files with selected UniprotIDs, gene names and the corresponding cell names/types

import csv
import glob
import numpy as np

# Constants used for directory names
#
# raw_data - file with raw data
# gen_data - file with generated data
# imp_data - file with data imported from uniprot.org
# write_cell_dir - additional extension for file to write data to
# cell_type_dirs - blanket vector for additional extensions required
raw_data_dir = './raw_data/'
gen_data_dir = './gen_data/'
imp_data_dir = './imp_data/'
write_cell_dir = 'adult_aged/'
cell_type_dirs = ['adult/', 'aged/']

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

# Generate the following four lists:
#
# cell_list - ordered list of cells including the run detail for each cell type (ex. records 'AgedHSC_Run1t1')
# unique_cell_list - ordered list of cells with only the cell type (ex. records 'AgedHSC' for 'AgedHSC_Run1t1')
# text_list - ordered list of text file names for all cells used
# prot_list - ordered list of unique UniprotIDs
def cell_prot_list(prot_gene_map, gene_prot_map, prot_aliases):
	cell_list = []; unique_cell_list = []; text_list = []; prot_list = []

	# Get order-matched lists of text files, cells and unique cells
	# Get a unique list of all expressed UniprotIDs
	for add in cell_type_dirs:
		for txt_file in glob.glob(raw_data_dir + add + '*.txt'):
	# Add to lists of text file, cells and unique cells
			text_list.append(txt_file)
			cell_list.append(txt_file[txt_file.find('\\') + 1:txt_file.find('Run') - 1] + txt_file[txt_file.find('Run') + 3:-4])
			unique_cell_list.append(txt_file[txt_file.find('\\') + 1:txt_file.find('Run') - 1])
			with open(txt_file, 'rb') as f:
				reader = csv.reader(f, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
	# Get indices for each respective data column
				for row in reader:
					UID = row.index('UniprotID')
					GID = row.index('GeneName')
					IID = row.index('Intensity')
					SID = row.index('NumofSpectra')
					break;
				for row in reader:
					prot = row[UID].strip()
					gene = row[GID].strip()
					intensity = row[IID].strip()
					spectra = row[SID].strip()
	# Only consider genes in the raw file that have some level of expression
					if intensity != '' and spectra != '':
	# Find the UniprotID that is used in prot_gene_map by searching through prot_aliases
						if prot != '' and prot != 'nan':
							FOUND = True
							if prot not in prot_gene_map: # Look through all {UniprotIDs}
								for i in range(len(prot_aliases) / 2): 
									index = 2 * i
									if prot in prot_aliases[index]: # Find the equivalent alias that is used in prot_aliases
										for ele in prot_aliases[index]: 
											if ele in prot_gene_map:
												prot = ele # Correct alias found
												break
										break
	# Debugging purposes
								if prot not in prot_gene_map:
									print('Debugger: Logical error; a proper alias for a UniprotID cannot be found.')
									print('Debugger: UniprotID - ' + prot)
									FOUND =  False
							if FOUND:
								prot_list.append(prot)
	# If UniprotID is invalid, but gene name is valid, find the UniprotID that is used in prot_gene_map by searching through prot_aliases
						elif gene != '' and gene != 'nan' and gene != '2 SV':
							FOUND = True
							if gene not in gene_prot_map:
								for i in range(len(prot_aliases) / 2):
									index = (2 * i) + 1
									if gene in prot_aliases[index]:
										for ele in prot_aliases[index]:
											if ele in gene_prot_map:
												gene = ele
												break
										break
								if gene not in gene_prot_map:
									print('Debugger: Logical error; a proper alias for a gene name cannot be found.')
									print('Debugger: gene name - ' + gene)
									FOUND = False
							if FOUND:
								prot_list.append(gene_prot_map[gene])
	# Create a non-redundant, sorted list of UniprotIDs that appear in the raw files
	prot_list = list(set(prot_list))
	prot_list.sort()
	
	return cell_list, unique_cell_list, text_list, prot_list

# Using the cell, UniprotID and gene name lists generated, 
# return a matrix of intensity values for each protein per cell type
def read_in_data(cell_list, unique_cell_list, text_list, prot_list, prot_gene_map, gene_prot_map, prot_aliases):
	# Matrix of values order-matched to cells and UniprotIDs
	i_value_matrix = np.zeros((len(prot_list), len(cell_list)), dtype = float)

	# Go through all raw files
	for txt_file in text_list:
		with open(txt_file, 'rb') as csvfile:
			reader = csv.reader(csvfile, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
			cell_index = text_list.index(txt_file) # Index of cell file
	# Debugging purposes to ensure data is lined up properly
			if (cell_list[cell_index] != txt_file[txt_file.find('\\') + 1:txt_file.find('Run') - 1] + txt_file[txt_file.find('Run') + 3:-4] 
				or unique_cell_list[cell_index] != txt_file[txt_file.find('\\') + 1:txt_file.find('Run') - 1]):
				print('Debugger: There is a mismatch of data ordering')
				print('Debugger: ' + txt_file)
				print('Debugger: ' + cell_list[cell_index])
				print('Debugger: ' + unique_cell_list[cell_index])
			for row in reader:
				UID = row.index('UniprotID')
				GID = row.index('GeneName')
				IID = row.index('Intensity')
				SID = row.index('NumofSpectra') # Needed to make sure spectra exists
				break;
			for row in reader:
				prot = row[UID].strip()
				gene = row[GID].strip()
				intensity = row[IID].strip()
				spectra = row[SID].strip()
	# Only consider data where intensity and spectra exist & where there is/are valid UniprotID and/or gene name
				if intensity != '' and spectra != '':
					if prot != '' and prot != 'nan':
						FOUND = True
						if prot not in prot_list:
							for i in range(len(prot_aliases) / 2):
								index = 2 * i
								if prot in prot_aliases[index]:
									for ele in prot_aliases[index]:
										if ele in prot_list:
											prot = ele
											break
									break
							if prot not in prot_list:
								FOUND = False
	# If the appropriate UniprotID alias can be matched, record the intensity value
						if FOUND:
							i_value_matrix[prot_list.index(prot)][cell_index] = float(intensity)
					elif gene != '' and gene != 'nan' and gene != '2 SV':
						FOUND = True
						if gene not in gene_prot_map:
							for i in range(len(prot_aliases) / 2):
								index = (2 * i) + 1
								if gene in prot_aliases[index]:
									for ele in prot_aliases[index]:
										if ele in gene_prot_map:
											gene = ele
											break
									break
							if gene not in gene_prot_map:
								FOUND = False
	# If the appropriate gene name alias can be matched, record the intensity value
						if FOUND:
							i_value_matrix[prot_list.index(gene_prot_map[gene])][cell_index] = float(intensity)

	return i_value_matrix

# Process data to four forms:
# 1) Data normalized to 1,000,000 per cell (i_value_matrix_norm)
# 2) Data averaged within cell types normalized to 1,000,000 per cell type (i_value_matrix_comb)
# 3) Data with the geometric norm within cell type, ignoring 0's, normalized to 1,000,000 per cell type (i_value_matrix_geom)
# 4) Data with the geometric norm within cell type, including 0's, normalized to 1,000,000 per cell type (i_value_matrix_zgeo)
# 5) Data averaged within cell type with each individual cell normalized to 1,000,000 (i_value_matrix_navg)
def process_data(i_value_matrix, condensed_unique_cell_list, unique_cell_list):
	# Normalize to 1,000,000
	i_value_matrix_norm = 1000000. * i_value_matrix / np.sum(i_value_matrix, axis = 0)

	# Average the data within the same cell type, then normalize the data to 1,000,000
	i_value_matrix_comb = []
	i_value_matrix_navg = []
	for i in range(len(condensed_unique_cell_list)):
		cell_type = condensed_unique_cell_list[i]
		indices = (np.array(unique_cell_list) == cell_type)
		new_data = np.sum(i_value_matrix[:, indices], axis = 1)
		new_data = 1000000 * new_data / np.sum(new_data)
		new_data_navg = 1000000 * i_value_matrix[:, indices] / np.sum(i_value_matrix, axis = 0)[indices]
		new_data_navg = np.mean(new_data_navg, axis = 1)
		i_value_matrix_comb.append(new_data)
		i_value_matrix_navg.append(new_data_navg)
	i_value_matrix_comb = np.transpose(np.squeeze(i_value_matrix_comb))
	i_value_matrix_navg = np.transpose(np.squeeze(i_value_matrix_navg))

	# Find the geometric mean within the same cell type, ignoring 0's, then normalize the data to 1,000,000
	# Find the geometric mean within the same cell type, including 0's, then normalize the data to 1,000,000
	i_value_matrix_geom = []
	i_value_matrix_zgeo = []
	for i in range(len(condensed_unique_cell_list)):
		cell_type = condensed_unique_cell_list[i]
		indices = (np.array(unique_cell_list) == cell_type)
		new_data = i_value_matrix[:, indices]
		new_data_TF = (new_data > 0)
		new_data_form = np.zeros(len(new_data), dtype = float)
		new_data_form_zero = np.zeros(len(new_data), dtype = float)
		for i in range(len(new_data)):
			row = new_data[i, :]
			row_TF = new_data_TF[i, :]
			if np.sum(row_TF.astype(int)) > 0:
				new_data_form[i] = np.prod(row[row_TF])**(1.0 / np.sum(row_TF.astype(int)))
				new_data_form_zero[i] = np.prod(row[row_TF])**(1.0 / 6.0)
		new_data_form = 1000000 * new_data_form / np.sum(new_data_form)
		new_data_form_zero = 1000000 * new_data_form_zero / np.sum(new_data_form_zero)
		i_value_matrix_geom.append(new_data_form)
		i_value_matrix_zgeo.append(new_data_form_zero)
	i_value_matrix_geom = np.transpose(np.squeeze(i_value_matrix_geom))
	i_value_matrix_zgeo = np.transpose(np.squeeze(i_value_matrix_zgeo))

	return i_value_matrix_norm, i_value_matrix_comb, i_value_matrix_geom, i_value_matrix_zgeo, i_value_matrix_navg

# Write all six data (raw, normalized, combined, geometric, zgeometric, naverage) 
# into their respective files.
def write_to_file(i_value_matrix, i_value_matrix_norm, i_value_matrix_comb, i_value_matrix_geom, i_value_matrix_zgeo, i_value_matrix_navg, 
	cell_list, condensed_unique_cell_list, prot_list, prot_gene_map):
	f1 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_raw.txt', 'w+')
	f2 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_norm.txt', 'w+')
	f3 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_comb.txt', 'w+')
	f4 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_geom.txt', 'w+')
	f5 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_zgeo.txt', 'w+')
	f6 = open(gen_data_dir + write_cell_dir + '3_gene_intensity_navg.txt', 'w+')

	f1.write(',,' + ','.join(cell_list))
	f2.write(',,' + ','.join(cell_list))
	f3.write(',,' + ','.join(condensed_unique_cell_list))
	f4.write(',,' + ','.join(condensed_unique_cell_list))
	f5.write(',,' + ','.join(condensed_unique_cell_list))
	f6.write(',,' + ','.join(condensed_unique_cell_list))

	for i in range(len(prot_list)):
		f1.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix[i])))
		f2.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix_norm[i])))
		f3.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix_comb[i])))
		f4.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix_geom[i])))
		f5.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix_zgeo[i])))
		f6.write("\n" + prot_list[i] + ","  +  prot_gene_map[prot_list[i]] + ',' + ",".join(map(str,i_value_matrix_navg[i])))

	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()

# Main function
def main():
	# Starting message
	print('Starting [3_data_compile]\n')

	# Retrieve mappings generated in 2_gene_prot_mappings.py
	prot_gene_map, gene_prot_map, prot_aliases = retrieve_mappings()

	# Generate the following four lists:
	#
	# cell_list - ordered list of cells including the run detail for each cell type (ex. records 'AgedHSC_Run1t1')
	# unique_cell_list - ordered list of cells with only the cell type (ex. records 'AgedHSC' for 'AgedHSC_Run1t1')
	# text_list - ordered list of text file names for all cells used
	# prot_list - ordered list of unique UniprotIDs
	cell_list, unique_cell_list, text_list, prot_list = cell_prot_list(prot_gene_map, gene_prot_map, prot_aliases)
	# List of cell types present
	condensed_unique_cell_list = np.unique(unique_cell_list)

	# Read in values from raw files
	i_value_matrix = read_in_data(cell_list, unique_cell_list, text_list, prot_list, prot_gene_map, gene_prot_map, prot_aliases)

	# Process data to four forms:
	# 1) Data normalized to 1,000,000 per cell (i_value_matrix_norm)
	# 2) Data averaged within cell types normalized to 1,000,000 per cell type (i_value_matrix_comb)
	# 3) Data with the geometric norm within cell type, ignoring 0's, normalized to 1,000,000 per cell type (i_value_matrix_geom)
	# 4) Data with the geometric norm within cell type, including 0's, normalized to 1,000,000 per cell type (i_value_matrix_zgeo)
	# 5) Data averaged within cell type with each individual cell normalized to 1,000,000 (i_value_matrix_navg)
	(i_value_matrix_norm, i_value_matrix_comb, i_value_matrix_geom, 
		i_value_matrix_zgeo, i_value_matrix_navg) = process_data(i_value_matrix, condensed_unique_cell_list, unique_cell_list)

	# Write the data into file
	write_to_file(i_value_matrix, i_value_matrix_norm, i_value_matrix_comb, i_value_matrix_geom, i_value_matrix_zgeo, i_value_matrix_navg, 
		cell_list, condensed_unique_cell_list, prot_list, prot_gene_map)

	# Closing message
	print('\n[3_data_compile] Complete')

main()