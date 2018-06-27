vcftools=/home/sobczm/bin/vcftools/bin

#Filter vcf outputs, only retain biallelic high-quality SNPS with no missing data for genetic analyses.

```bash
for vcf in $(ls *_contigs_unmasked.vcf)
do
    echo $vcf
    script=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/sub_vcf_parser.sh
    qsub $script $vcf
done
```

Only one of these can be run at a time!

#General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
cd ../
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/N.ditissima_contigs_unmasked.vcf > SNP_calling/N.ditissima_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/N.ditissima_contigs_unmasked_filtered.vcf > SNP_calling/N.ditissima_contigs_unmasked_filtered.stat
```

#Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls SNP_calling/*_filtered.vcf)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

#Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.

```bash
for vcf in $(ls SNP_calling/*_filtered.vcf)
do
    echo $vcf
    out=$(basename $vcf .vcf)
    echo $out
    $vcftools/vcftools --vcf $vcf --mac 1 --recode --out SNP_calling/$out
done
```

#Create custom SnpEff genome database
```bash
snpeff=/home/gomeza/prog/snpEff
nano $snpeff/snpEff.config
```

##Add the following lines to the section with databases:

#---
# EMR Databases
#----
# Hg199 genome
Hg199_minion.genome: Hg199

#Collect input files

```bash
mkdir -p $snpeff/data/Hg199_minion
cp Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa $snpeff/data/Hg199_minion
cp gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3 $snpeff/data/Hg199_minion
```
#Rename input files
```bash
cd $snpeff/data/Hg199_minion
mv final_genes_appended.gff3 genes.gff
mv N.ditissima_contigs_unmasked.fa sequences.fa
```
#Build database using GFF3 annotation

java -jar $snpeff/snpEff.jar build -gff3 -v Hg199_minion

#Annotate VCF files

```bash
input=/data/scratch/gomeza/analysis/popgen/
cd $input
for a in SNP_calling/*recode.vcf
do
    echo $a
    filename=$(basename "$a")
    java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 Hg199_minion $a > ${filename%.vcf}_annotated.vcf
    mv snpEff_genes.txt SNP_calling/snpEff_genes_${filename%.vcf}.txt
    mv snpEff_summary.html  SNP_calling/snpEff_summary__${filename%.vcf}.html
done
```

#Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls analysis/popgen/SNP_calling/*distance.log)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  Rscript --vanilla $scripts/distance_matrix.R $log
done
```
#Carry out PCA and plot the results

```bash
for vcf in $(ls analysis/popgen/SNP_calling/*filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    out=$(basename $vcf contigs_unmasked_filtered.vcf)
    echo $out
    Rscript --vanilla $scripts/pca.R $vcf $out
done
```
#Calculate an NJ tree based on all the SNPs. Outputs a basic diplay of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

```bash
cd analysis/popgen/SNP_calling
for vcf in $(ls *filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    $scripts/nj_tree.sh $vcf 1
done
```
This script isn't happy with the low levels of variation I have. TODO: Investigate other tree methods - Michelle mentioned some possibilities