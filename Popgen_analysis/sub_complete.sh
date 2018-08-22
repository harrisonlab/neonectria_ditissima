#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h=blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

input_dir=$1
n=$2
four_Ner_max=$3
no_points=$4
theta=$5

cd $input_dir

/home/gomeza/prog/LDhat-master/complete -n $n -rhomax $four_Ner_max -n_pts $no_points -theta $theta
