#!/usr/bin/python
#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import math
import numpy as np
from sets import Set
from collections import defaultdict
from collections import Counter
from operator import itemgetter
from Bio import SeqIO

ap = argparse.ArgumentParser()

ap.add_argument('--sig', required=True,type=str,help='Gregs Deseq2 output file showing LFC and P values in tsv format')
conf = ap.parse_args()

with open(conf.sig) as f:
    DEG_lines = f.readlines()

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

problem_genes = []

for line in DEG_lines[0:1]:
    line = line.rstrip()
    split_line = line.split("\t")
    condition_list = list(itemgetter(12,15,18,21,24,27,30,33)(split_line))
    # print condition_list
    condition_list = [x.replace('LFC_', '') for x in condition_list]

for line in DEG_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    if len(split_line) < 16:
        continue
    # print line
    LFCs = itemgetter(12,15,18,21,24,27,30,33)(split_line)
    Pvals = itemgetter(13,16,19,22,25,28,31,34)(split_line)
    Pvals = [p or 1 for p in Pvals]
    # if (
    #     any([float(p) < 0.05 for p in Pvals]) and
    #     any([abs(float(x)) > 2 for x in LFCs])
    #     ):
    #     print line
    DEG_conditions = []
    # print LFCs
    for condition, lfc, p in zip(condition_list, LFCs, Pvals):
        # print "\t".join([gene_id, condition, lfc, p])
        if (
            float(p) < 0.05 and
            abs(float(lfc)) > 2
            ):
                DEG_conditions.append(condition)
    if DEG_conditions:
        print "\t".join([gene_id, ";".join(DEG_conditions)])
