#!/usr/bin/python

'''
This script uses text files of upregulated genes to create a count table for
which genes are upregulated above a set threshold
'''

import argparse
from collections import defaultdict
import numpy
import csv

# -----------------------------------------------------
# Step 1
# Import variables, load input files & create list of genes
# If using a different number of files
# arguments & appending to list of genes will need to be changed
# -----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True, type=str, help='text file of genes \
with fold change information')
ap.add_argument('--out_dir', required=True, type=str, help='the tsv file where \
the count table is output to')
conf = ap.parse_args()

inp1_dict = defaultdict(list)
with open(conf.input) as f1:
    inp1_lines = f1.readlines()[1:]
    genes_list = []
    inp1 = []
    for x in inp1_lines:
        genes_list.append(x.split('\t')[0])
        inp1.append(x.split('\t')[0])
        gene_name = x.split('\t')[0]
        value = float(x.split('\t')[2])
        inp1_dict[gene_name].append(value)

# -----------------------------------------------------
# Step 2
# Load gene names to a numpy array and create new columns
# If doing with a different number of files
# change the number in the numpy.reshape() command
# -----------------------------------------------------

a = numpy.array(["Gene_Name", "Fold_Change"])

for x in genes_list:
    to_add = []
    to_add.append(x)
    try:
        b = inp1.index(x)
    except ValueError:
        to_add.append('0')
    else:
        for y in inp1_dict[x]:
            test = abs(y)
            if test > 1:
                to_add.append('1')
            else:
                to_add.append('0')
    a = numpy.append(a, to_add, axis=0)

z = len(genes_list) + 1
a = numpy.reshape(a, (z, 2))

outfile = str(conf.out_dir)
with open(outfile, "w") as o:
    writer = csv.writer(o, delimiter='\t')
    writer.writerows(a)
