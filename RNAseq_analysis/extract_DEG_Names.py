#!/usr/bin/python

'''
This script uses the tsv file produced by parsing RNA-Seq with a previous
python script to extract only DEGs above a certain threshold
'''

import argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True, type=str,
                help='tsv file of binary presence table of DEGs')
ap.add_argument('--output', required=True, type=str,
                help='Output text file of gene names')
conf = ap.parse_args()

# -----------------------------------------------------
# Step 1
# Load .tsv file into a dictionary and create gene list to iterate over
# -----------------------------------------------------

DEG_dict = defaultdict(list)
gene_list = []
with open(conf.input) as f:
    DEG_lines = f.readlines()
    for line in DEG_lines:
        if line.startswith('Gene_Name'):
            continue
        else:
            split_line = line.split()
            gene_name = split_line[0]
            gene_list.append(gene_name)
            time_a = split_line[1]
            time_b = split_line[2]
            time_c = split_line[3]
            DEG_dict[gene_name].extend([time_a, time_b, time_c])

# -----------------------------------------------------
# Step 2
# Test for presence of a '1' in the dictionary entry for each gene,
# denotating it's above the threshold value
# -----------------------------------------------------

positive_hits = []
for x in gene_list:
    if '1' in DEG_dict[x]:
        positive_hits.append(x)

# -----------------------------------------------------
# Step 3
# Write list of positive hits to text file
# -----------------------------------------------------

positive_hits = map(lambda x: x+'\n', positive_hits)
outfile = str(conf.output)
with open(outfile, "w") as o:
    o.writelines(positive_hits)
    o.close()
