#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=16G
#$ -l h=blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

filename=$1
config=$2

pgdspid=/home/gomeza/bin/PGDSpider_2.1.0.3
java -Xmx15240m -Xms5120m -jar $pgdspid/PGDSpider2-cli.jar -inputfile $filename -inputformat VCF -outputfile ${filename%.vcf}.geno -outputformat GESTE_BAYE_SCAN -spid $config
