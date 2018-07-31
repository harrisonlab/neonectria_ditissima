#!/bin/bash
input=analysis/popgen/SNP_calling
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp

input_hap=analysis/popgen/SNP_calling/N.ditissima_contigs_unmasked.vcf

#Create a directory for each individual BayeScan runS
mkdir -p $input/outliers/N.ditissima

#Copy input files for the analysis
cp -r $input_hap $input/outliers/N.ditissima
cd $input/outliers

#################### 1) Filter SNPs to retain only biallelic SNPs, otherwise not compatible with PGDSpider
#and Bayescan. Furthermore, keep only SNPs with max 5% missing genotypes.
### It may be necessary to also pre-filter for samples which were poorly sequenced/
#aligned beforehand to avoid removal of too many potentially informative SNPs.
vcftools=/home/sobczm/bin/vcftools/bin
for filename in $input/outliers/N.ditissima/*.vcf
do
$vcftools/vcftools --vcf $filename --max-missing 0.95 --mac 1 --min-alleles 2 --max-alleles 2 --recode --out ${filename%.vcf}_bi_filtered
done

############# 2) PGDSpider conversion from VCF to Bayescan format
#First need to prepare a simple file with custom assignment of each individual sample
#name to one of the populations identified using PCA, STRUCTURE analysis etc.
#Format: Sample_name Population_name
#But first: double check what sample names are present in the following VCF file:
filename=$input/outliers/N.ditissima/N.ditissima_contigs_unmasked_bi_filtered.recode.vcf
grep "#CHROM" $filename | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${filename%.vcf}.lst"

#Resulting example file population assignment file:
cat $input/outliers/N.ditissima/N.ditissima_contigs_unmasked_bi_filtered.recode.lst

#Now, need to prepare the configuration file. Copy the conversion script
pgdspid=/home/gomeza/bin/PGDSpider_2.1.0.3
#For haploid:
config=vcf_to_bayescan_haploid.spid~
cp $pgdspid/$config ./

#Now, need to edit the conversion script to include the population assignment file
pop_assignment_file=$input/outliers/N.ditissima/N.ditissima_contigs_unmasked_bi_filtered.recode.lst
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$pop_assignment_file"',' $config

#Convert to Bayescan input format with PGDSpider.
#For bigger input VCF file (>50 MB) use the script to convert to BayeScan input (will take a couple of hours):
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/
qsub $scripts/sub_pgdspider.sh $filename $config
dos2unix ${filename%.vcf}.geno

#This is the command to run it on the head node. But requires memory, so it is better qsub it.
#java -jar -Xmx1024m -Xms512m $pgdspid/PGDSpider2-cli.jar -inputfile $filename -inputformat VCF -outputfile ${filename%.vcf}.geno -outputformat GESTE_BAYE_SCAN -spid $config
#dos2unix ${filename%.vcf}.geno

############# 3) Bayescan analysis
#Bayescan run and plot the results in R.
qsub -R y $scripts/sub_bayescan.sh ${filename%.vcf}.geno
