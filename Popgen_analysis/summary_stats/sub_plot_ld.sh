#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/summary_stats
input=$1

/home/adamst/prog/R/R-3.2.5/bin/Rscript --vanilla $scripts/plot_ld.R $input
