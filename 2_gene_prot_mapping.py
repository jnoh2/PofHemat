# File name: 2_gene_prot_mapping.py
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
# Using the mappings imported from uniprot.org and the raw files,
#
# 1) records all UniprotID <-> gene name mappings that occur
# 2) recursively follows these mappings
# 3) groups genes/proteins that have multiple alises (ex. p53, Tp53, Trp53)
# 4) maps between groups of UniprotIDs and gene names and writes it to a file, and
# 5) selects a single UniprotID alias and gene name alias from each mapping as standard names to be used in downstream analyses.
#
# This method creates a final, singular list that ensures no redundancy & confusion when comparing
# between samples that may have used different aliases for the same gene/protein. Unmapped UniprotIDs
# and/or gene names are given mappings to new aliases (UNM# for UniprotID and Unm# for gene name)

import csv
import glob
import numpy as np

# Constants used for directory names
#
# raw_data - file with raw data
# gen_data - file with generated data
# imp_data - file with data imported from uniprot.org
# cell_type_dirs - blanket vector for additional extensions required

raw_data_dir = './raw_data/'
gen_data_dir = './gen_data/'
imp_data_dir = './imp_data/'
cell_type_dirs = ['adult/', 'aged/']

# Write results into file
# 2_prot_gene_group - the set to set mappings between UniprotIDs and gene names
# 2_prot_gene_single - the UniprotID : gene name mappings created to be used in downstream analyses
# 2_mapping_summary - summary of the mapping process
def write_to_file(prot_set_gene_set, prot_gene_dict, cell_list, UNMs, Unms):
	w1 = open(gen_data_dir + '2_prot_gene_group.txt', 'w+')
	w2 = open(gen_data_dir + '2_prot_gene_single.txt', 'w+')
	w3 = open(gen_data_dir + '2_total_mapping_summary.txt', 'w+')

	# Write all set to set mappings into w1
	for i in range(len(prot_set_gene_set) / 2):
		index = 2 * i
		prots = prot_set_gene_set[index]
		genes = prot_set_gene_set[index + 1]
	# Format- UniprotID1,UniprotID2,UniprotID3,PROT:GENE,gene name1,gene name2,gene name3
		w1.write(','.join(prots) + ',PROT:GENE,' + ','.join(genes) + '\n')
	# Write all UniprotID : gene name mappings into w2
	for key, value in prot_gene_dict.items():
	# Format- UniprotID:gene name
		w2.write(key + ',' + value + '\n')

	# Summarize mapping process into w3
	cell_list = np.array(cell_list, dtype = str) # Converted to numpy object for ease of manipulation
	w3.write('Cell types used to compile mappings: ' + ', '.join(np.unique(cell_list)) + '\n\n')
	w3.write('Counts of each cell type: ' + '\n')
	for ele in np.unique(cell_list):
		w3.write(ele + ' - ' + str(np.sum((cell_list == ele).astype(int))) + '\n')
	w3.write('\nTotal unique UniprotID set : Gene name set - ' + str(len(prot_set_gene_set) / 2))
	w3.write('\n\nTotal UNM\'s (new UniprotIDs) - ' + str(len(UNMs)))
	w3.write('\n\n' + ', '.join(UNMs))
	w3.write('\n\nTotal Unm\'s (new gene names) - ' + str(Unms))
	w3.write('\n\n' + ', '.join(Unms))

	# Close all files
	w1.close()
	w2.close()
	w3.close()

# Recursively defined function to create non-redundant mappings between
# UniprotID aliases with gene name alises.
# prot/gene - contains the next UniprotID/gene name to begin another node from
# used_prot/used_gene - tracks the list of prot/gene previously used to prevent redundant recursion
# prot_set/gene_set - tracks the UniprotID/gene name alises of the same protein/gene name
# prot_gene_set/gene_prot_set - the dictionaries to recurse off of
def recurse_gene_prot(prot, gene, used_prot, used_gene, prot_set, gene_set, prot_gene_set, gene_prot_set):
	if prot == -1: # Recurse from a gene node
		if gene not in used_gene: # Base case
			used_gene.add(gene)
			gene_set.add(gene) # Adds the gene name to the group of aliases
			proteins = gene_prot_set[gene] # Searches for the protein nodes to recurse towards
			for ele in proteins:
				recurse_gene_prot(ele, -1, used_prot, used_gene, prot_set, gene_set, prot_gene_set, gene_prot_set)
	else: # Recurse from a protein node
		if prot not in used_prot: # Base case
			used_prot.add(prot)
			prot_set.add(prot) # Adds the UniprotID to the group of aliases
			genes = prot_gene_set[prot] # Searches for the gene nodes to recurse towards
			for ele in genes:
				recurse_gene_prot(-1, ele, used_prot, used_gene, prot_set, gene_set, prot_gene_set, gene_prot_set)

# Creates {UniprotIDs} : {gene names} as well as 
# UniprotID : gene name to be used in the downstream analyses
def create_set_mappings(prot_gene_set, gene_prot_set, prot_set_gene_set, prot_gene_dict):
	# Sets tracking the UniprotIDs/gene names used throughout recursion
	used_prot = set()
	used_gene = set()

	# Iterate through all prot_gene_set mappings to recursively find
	# mappings between UniprotID aliases and gene name alises
	for key, value in prot_gene_set.items():
		if key not in used_prot:
	# prot_set/gene_set keeps track of aliases that represent the same protein/gene
			prot_set = set()
			gene_set = set()
			first_prot = key # first_prot is selected as the UniprotID that will be used in downstream analyses
			used_prot.add(key) # Important for maintaining base case of recursion
			prot_set.add(key) # Updates the set
			for ele in value:
				first_gene = ele # first_gene is selected as the gene name that will be used in downstream analyses
				break;
	# For a given UniprotID, recurse through every gene name present
			for ele in value:
				recurse_gene_prot(-1, ele, used_prot, used_gene, prot_set, gene_set, prot_gene_set, gene_prot_set)
	# Record the {UniprotIDs} : {gene names} as well as UniprotID : gene name
			prot_set_gene_set.append(prot_set)
			prot_set_gene_set.append(gene_set)
			prot_gene_dict[first_prot] = first_gene

	# Double check logic with gene list to ensure all gene names were utilized once
	ALL_USED = True
	for key in gene_prot_set:
		if key not in used_gene:
			print('Debugger: ' + key + ' has not been mapped by the recursive function.')
			print('Debugger: To ensure that no data was lost, look into why this is the case')
			ALL_USED = False
	if ALL_USED:
		print('Debugger: All gene names were mapped at least once; no data has been lost throughout recursion')

# For all UniprotIDs and gene names that have not been mapped by uniprot.org,
# if they have not yet been mapped based on the raw files, assign them to new
# UniprotIDs/gene names (UNMs/Unms)
def name_unmapped(prot_gene_set, gene_prot_set, unmap_list, UNMs, Unms):
	counter1 = 0 # Tracks the Unm number
	with open(imp_data_dir + '1_prot_gene_unmapped.txt', 'rb') as r:
		reader = csv.reader(r, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
	# Go through all UniprotIDs not mapped by uniprot.org
		for row in reader:
			prot = row[0]
	# If the UniprotID has not been mapped before
			if prot not in prot_gene_set:
				counter1 += 1
				gene = 'Unm' + str(counter1) # New gene name to be mapped to
	# Create mappings with new gene name
				prot_gene_set[prot] = set((gene,))
				gene_prot_set[gene] = set((prot,))
				Unms.append(prot)
	# unmap_list is maintained for debugging purposes
				while prot in unmap_list:
					unmap_list.remove(prot)
				
	counter2 = 0 # Tracks the UNM number
	with open(imp_data_dir + '1_gene_prot_unmapped.txt', 'rb') as r:
		reader = csv.reader(r, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
	# Go through all gene names not mapped by uniprot.org
		for row in reader:
			gene = row[0]
	# If the gene name has not been mapped before
			if gene not in gene_prot_set:
				counter2 += 1
				prot = 'UNM' + str(counter2) # New UniprotID to be mapped to
				prot_gene_set[prot] = set((gene,))
				gene_prot_set[gene] = set((prot,))
				UNMs.append(gene)
	# unmap_list is maintained for debugging purposes
				while gene in unmap_list:
					unmap_list.remove(gene)

	# Remove all mapped UniprotIDs and gene names from unmap_list, if still present
	for key in prot_gene_set:
		while key in unmap_list:
			unmap_list.remove(key)
	for key in gene_prot_set:
		while key in unmap_list:
			unmap_list.remove(key)

	# This message is for debugging purposes
	print('Debugger: List of unmapped genes and proteins - ' + str(unmap_list))
	print('Debugger: If the list above is empty, then the mapping has occurred successfully')

# Adds UniprotID : gene name & gene name : UniprotID as identified in the raw data files
# Any defined UniprotID mapped to '2 SV', 'nan' or '' are put into unmap_list
# Any defined gene name mapped to 'nan' or '' are put into unmap_list
def find_mappings(prot_gene_set, gene_prot_set, unmap_list, cell_list):
	# Read files from both adult and aged
	for add in cell_type_dirs:
		for txt_file in glob.glob(raw_data_dir + add + '*.txt'):
	# Record the cell used
			cell_list.append(txt_file[txt_file.find('\\') + 1:txt_file.find('Run') - 1]) # Ex. \AgedHSC_Run2t1.txt -> AgedHSC
			with open(txt_file, 'rb') as csvfile:
				reader = csv.reader(csvfile, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
	# Identify the indices for UniprotID and gene name
				for row in reader:
					UID = row.index('UniprotID')
					GID = row.index('GeneName')
					break;
				for row in reader:
					prot = row[UID].strip()
					gene = row[GID].strip()
	# If both UniprotID and gene name are valid, record in the mapping
					if prot != '' and prot != 'nan' and gene != '' and gene != 'nan' and gene != '2 SV':
						if prot in prot_gene_set: # This UniprotID has been encountered previously
							prot_gene_set[prot].add(gene)
						else:
							prot_gene_set[prot] = set((gene,))
						if gene in gene_prot_set: # This gene name has been encountered previously
							gene_prot_set[gene].add(prot)
						else:
							gene_prot_set[gene] = set((prot,))
	# If only one of the two is invalid, record the valid one in the unmap_list
					else:
						if prot != '' and prot != 'nan' and prot not in prot_gene_set:
							unmap_list.append(prot)
						elif gene != '' and gene != 'nan' and gene != '2 SV' and gene not in gene_prot_set:
							unmap_list.append(gene)	

# Adds UniprotID : {gene name set} as imported from uniprot.org
# Adds gene name : {UniprotID set} as imported from uniprot.org
def import_mappings(prot_gene_set, gene_prot_set):
	# UniprotID : {gene name set}
	with open(imp_data_dir + '1_prot_gene_mapped.txt', 'rb') as r:
		reader = csv.reader(r, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
		for row in reader:
			prot = row[0]
			gene = row[1]
	# Adds mapping in both directions
			if prot in prot_gene_set: # UniprotID has already been mapped before
				prot_gene_set[prot].add(gene)
			else:
				prot_gene_set[prot] = set((gene,))
			if gene in gene_prot_set: # Gene name has already been mapped before
				gene_prot_set[gene].add(prot)
			else:
				gene_prot_set[gene] = set((prot,))

	# gene name : {UniprotID set}
	# Format in data file - gene name1, gene name2, gene name 3\tUniprotID
	with open(imp_data_dir + '1_gene_prot_mapped.txt', 'rb') as r:
		reader = csv.reader(r, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_NONE)
		for row in reader:
	# Find and record all gene names in the row
			gene = row[0]
			genes = []
			if gene.find(',') != -1: # If more than one gene name exists
				while True: # Repeat until all gene names have been added
					index = gene.find(',')
					if index >= 0: # -1 if ',' is no longer contained; what's leftover is the final gene name
						genes.append(gene[:index])
						gene = gene[index + 1:]
					else:
						genes.append(gene)
						break
			else:
				genes = [gene]
	# Find and record the mapped UniprotID
			prot = row[1]
	# Adds mapping in both directions
			if prot in prot_gene_set: # UniprotID has already been mapped before
				for gene in genes:
					prot_gene_set[prot].add(gene)
			else: 
				for gene in genes:
					prot_gene_set[prot] = set((gene,))
					break
				for gene in genes:
					prot_gene_set[prot].add(gene)
			for gene in genes:
				if gene in gene_prot_set: # Gene name has already been mapped before
					gene_prot_set[gene].add(prot)
				else:
					gene_prot_set[gene] = set((prot,))

# Main function
def main():
	# Starting message
	print('Starting [2_gene_prot_mapping]\n')

	# prot_gene - set of all UniprotID : {gene name set} mappings
	# gene_prot - set of all gene name : {UniprotID set} mappings
	prot_gene_set = {}
	gene_prot_set = {}

	# Add mappings imported from uniprot.org
	import_mappings(prot_gene_set, gene_prot_set)

	# Potential list of unmapped UniprotIDs and gene names as identified in the raw files
	# This list is only maintained for debugging purposes
	unmap_list = []
	# Keeps track of cells (raw files) used to compile mappings (does not include run #)
 	cell_list = []

	# Finds mappings (and lack of mappings) between UniprotIDs and gene names as seen in the raw files
	find_mappings(prot_gene_set, gene_prot_set, unmap_list, cell_list)
	
	# List of UniprotIDs mapped to new gene names (Unms) and
	# list of gene names mapped to new UniprotIDs (UNMs)
	UNMs = []
	Unms = []

	# For all UniprotIDs and gene names that have not been mapped by uniprot.org
	# or by the raw files, creae new UniprotIDs or gene names and use those for mapping
	name_unmapped(prot_gene_set, gene_prot_set, unmap_list, UNMs, Unms)

	# Vector of UniprotID alias sets and gene name alias sets
	# Given len(prot_set_gene_set) = n, where n is even
	# the (2k)th UniprotID alias set is mapped to (2k + 1)th gene name alias set
	# where 0 <= k <= (n / 2)
	prot_set_gene_set = []
	# Dictionary of UniprotID : gene name, such that a single UniprotID and gene name
	# are selected for each set to set mapping determined in prot_set_gene_set
	prot_gene_dict = {}

	# Creates {UniprotIDs} : {gene names} as well as 
	# UniprotID : gene name to be used in the downstream analyses
	create_set_mappings(prot_gene_set, gene_prot_set, prot_set_gene_set, prot_gene_dict)

	# Write results into file
	# 2_prot_gene_group - the set to set mappings between UniprotIDs and gene names
	# 2_prot_gene_single - the UniprotID : gene name mappings created to be used in downstream analyses
	# 2_mapping_summary - summary of the mapping process
	write_to_file(prot_set_gene_set, prot_gene_dict, cell_list, UNMs, Unms)
			
	# Closing message
	print('\n[2_gene_prot_mapping] Complete')

main()