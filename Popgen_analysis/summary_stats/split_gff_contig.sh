#!/bin/bash

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
#1st argument: Whole-genome GFF file.
gff=$1
filename=$(basename "$gff")
#First sort the gff file by chromosomes
sort -n -k 1 $gff >${filename%.gff*}_sorted.gff
#And then split
python $scripts/summary_stats/split_gff.py ${filename%.gff*}_sorted.gff
