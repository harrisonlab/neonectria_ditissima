#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace11.blacklace

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis

/home/adamst/prog/R/R-3.2.5/bin/Rscript --vanilla $scripts/Export_Network_Cytoscape.R --out_dir $1 --module $2
