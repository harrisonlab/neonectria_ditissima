#Analysis of the summary stats section from Maria's github

##Copy input for the analysis into a new directory

```bash
input=/data/scratch/gomeza/analysis/summary_stats2
snpeff=/home/gomeza/bin/snpEff
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
```

##All individuals

```bash
cp /data/scratch/gomeza/analysis/popgen/SNP_calling3/Hg199_contigs_unmasked.vcf $input
cp /data/scratch/gomeza/analysis/popgen/SNP_calling3/Hg199_contigs_unmasked_filtered.vcf  $input
cp /data/scratch/gomeza/repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa $input
cp /data/scratch/gomeza/gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3 $input
cd $input
```

##Create additional subsets of VCF files with reduced number of individuals

```bash
vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin

#All N.ditissima
#First argument: unfiltered input VCF file with all SNPs
#Subsequent arguments: Sample names of individuals to be removed
$vcflib/vcfremovesamples Hg199_contigs_unmasked.vcf NMaj > Hg199_contigs_unmasked_allNd.vcf
#Filter the SNPs
$scripts/snp/vcf_parser_haploid.py --i Hg199_contigs_unmasked_allNd.vcf

#Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.
$vcftools/vcftools --vcf Hg199_contigs_unmasked_allNd_filtered.vcf --mac 1 --recode --out N.ditissima_contigs_unmasked_allNd_filtered

After filtering, kept 88126 out of a possible 493419 Sites

#Only Northwestern countries pathogens
$vcflib/vcfremovesamples N.ditissima_contigs_unmasked.vcf BGV344 ND8 ND9 OPC304 P112 R68-17 > N.ditissima_contigs_unmasked_NorthWest.vcf
$scripts/vcf_parser_haploid.py --i N.ditissima_contigs_unmasked_NorthWest.vcf
$vcftools/vcftools --vcf N.ditissima_contigs_unmasked_NorthWest_filtered.vcf --mac 1 --recode --out N.ditissima_contigs_unmasked_NorthWest_filtered

#Only highly pathogenic tested
$vcflib/vcfremovesamples N.ditissima_contigs_unmasked.vcf Ag04 Ag05 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 ND8 ND9 OPC304 P112 R0905 R37-15 R68-17 > N.ditissima_contigs_unmasked_patho.vcf
$scripts/vcf_parser_haploid.py --i N.ditissima_contigs_unmasked_patho.vcf
$vcftools/vcftools --vcf N.ditissima_contigs_unmasked_patho_filtered.vcf --mac 1 --recode --out N.ditissima_contigs_unmasked_patho_filtered

#Only non-pathogens tested
$vcflib/vcfremovesamples N.ditissima_contigs_unmasked.vcf Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 ND8 ND9 Hg199 OPC304 P112 R0905 R39-15 R41-15 R42-15 R45-15 R6-17-3 R68-17 > N.ditissima_contigs_unmasked_nopatho.vcf
$scripts/vcf_parser_haploid.py --i N.ditissima_contigs_unmasked_nopatho.vcf
$vcftools/vcftools --vcf N.ditissima_contigs_unmasked_nopatho_filtered.vcf --mac 1 --recode --out N.ditissima_contigs_unmasked_nopatho_filtered
```

##Create custom SnpEff genome database

```bash
$scripts/summary_stats/build_genome_database.sh N.ditissima_contigs_unmasked.fa final_genes_appended.gff3 Hg199_minion
```

#Annotate VCF files

```bash
for a in *recode.vcf
do
$scripts/summary_stats/annotate_snps_genome.sh $a Hg199_minion
done
```
```bash
###Create FASTA alignment files containing only select subsets of SNPs. Required
### for analyses in the fus_popgenome_analysis.sh script. Here, using option 1 as haploid organism, but for diploid organisms use
### typically option 2 (for Popgenome analysis) or 3.
### From now onwards, analysing the dataset without R68.
cd $input/NoR68
ref_genome=/data/scratch/gomeza/analysis/summary_stats/N.ditissima_contigs_unmasked.fa
#All
python $scripts/summary_stats/vcf_to_fasta.py N.ditissima_contigs_unmasked_noR68_filtered.recode_annotated.vcf $ref_genome 1
#Moving each subset of FASTA files into a separate dir.
mkdir all
mv *.fasta ./all
##coding
python $scripts/summary_stats/vcf_to_fasta.py N.ditissima_contigs_unmasked_noR68_filtered.recode_coding.vcf $ref_genome 1
mkdir coding
mv *.fasta ./coding
##silent(four-fold degenerate)
python $scripts/summary_stats/vcf_to_fasta.py N.ditissima_contigs_unmasked_noR68_filtered.recode_syn_4fd.vcf $ref_genome 1
mkdir silent
mv *.fasta ./silent
```
