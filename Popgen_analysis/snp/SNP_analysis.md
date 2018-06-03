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
