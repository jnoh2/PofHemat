# File name: 1_raw_gene_prot.py
# Last Edit: 10/22/2019
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
# Reads all raw files and extracts the UniprotIDs and gene names that appear

import csv
import glob

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

# Main function
def main():
	# Starting message
	print('Starting [1_raw_gene_prot]\n')

	# Contains the UniprotIDs and gene names
	prot_list = set()
	gene_list = set()

	# w1 - UniprotID list path
	# w2 - gene name list path
	w1 = open(gen_data_dir + '1_raw_prot_list.txt', 'w+')
	w2 = open(gen_data_dir + '1_raw_gene_list.txt', 'w+')

	# Read from all raw data files
	for add in cell_type_dirs:
		for txt_file in glob.glob(raw_data_dir + add + '*.txt'):
			with open(txt_file, 'rb') as f:
				reader = csv.reader(f, delimiter = '\t')

	# Get index of UniprotID and gene name
				for row in reader:
					UID = row.index('UniprotID')
					GID = row.index('GeneName')
					break

	# Add each UniprotID and gene name to their respective lists
	# 'nan', '' and '2 SV' are skipped, as they are not real UniprotIDs or gene names
				for row in reader:
					prot = row[UID].strip()
					gene = row[GID].strip()
					if prot != 'nan' and prot != '' and prot not in prot_list:
						prot_list.add(prot)
						w1.write(prot + '\n')
					if gene != 'nan' and gene != '' and gene != '2 SV' and gene not in gene_list:
						gene_list.add(gene)
						w2.write(gene + '\n')

	# Close all files
	w1.close()
	w2.close()

	# Closing message
	print('\n[1_raw_gene_prot] Complete')	

main()
