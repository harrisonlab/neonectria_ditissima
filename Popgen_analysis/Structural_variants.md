#Looks for structural variations between the genomes
##Structural variants can include: duplications, deletions, inversions & translocations. Uses read-pair configuration, split-reads & read-depth

###First, run BWA-mem to align reads to Reference genome assembly before SV calling

####Run BWA-mem

```bash
Reference=repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 ND8 ND9 OPC304 P112 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj R0905_all Hg199; do
for Reads in $(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/qc_dna/paired/N.*/$Strain)
do
  Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    FRead=$Reads/F/*.fq.gz
    RRead=$Reads/R/*.fq.gz
    OutDir=alignment/bwa
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa.sh $Strain $Reference $FRead $RRead $OutDir
done
done
```

R0905 genome as reference

```bash
Reference=R0905_good/R0905_contigs_unmasked.fa
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 ND8 ND9 OPC304 P112 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 NMaj R0905_all Hg199; do
for Reads in $(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/qc_dna/paired/N.*/$Strain)
do
  Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    FRead=$Reads/F/*.fq.gz
    RRead=$Reads/R/*.fq.gz
    OutDir=alignment/bwa_vsR0905
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa.sh $Strain $Reference $FRead $RRead $OutDir
done
done
```

####Run svaba

```bash
Prefix=Nd_svaba
Reference=repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa
AlignDir=alignment/bwa
OutDir=svaba
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
qsub $ProgDir/sub_svaba.sh $Prefix $Reference $AlignDir $OutDir
```

### Analysis of files produced by svaba

#### Analysis of indel file

This file contains smaller indels

##### Set initial variables

```bash
scripts=/home/sobczm/bin/popgen/summary_stats
input=/data/scratch/gomeza/analysis/sv_calling2/svaba
```

##### Create a cut-down VCF and filter it

```bash
cd $input

vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Nd_svaba_sv.svaba.indel.vcf NMaj_sorted.bam > Nd_svaba_sv.svaba.indel_cut.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Nd_svaba_sv.svaba.indel_cut.vcf  --max-missing 0.95 --recode --out Nd_svaba_sv.svaba.indel_cut_filtered
```

No sites removed by filtering

#### Ancestral variants

##### Set UK isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Nd_svaba_sv.svaba.indel_cut_filtered.recode.vcf --out Nd_svaba_sv.svaba.indel_filtered_UKvsIT.vcf --ply 1 --pop1 Ag02_sorted.bam,,Ag04_sorted.bam,,Ag05_sorted.bam,,Ag06_sorted.bam,,Hg199_sorted.bam,,R0905_sorted.bam,,R6-17-2_sorted.bam,,R6-17-3_sorted.bam --pop2 R68-17-C2_sorted.bam,,R68-17-C3_sorted.bam --thr 0.95
```

```bash
python $scripts/vcf_find_difference_pop.py --vcf Nd_svaba_sv.svaba.indel_cut_filtered.recode.vcf --out Nd_svaba_sv.svaba.indel_filtered_Pop1vsPop2.vcf --ply 1 --pop1 Ag02_sorted.bam,,Ag04_sorted.bam,,Ag05_sorted.bam,,Ag06_sorted.bam,,Ag08_sorted.bam,,Ag09_A_sorted.bam,,Ag11_A_sorted.bam,,Ag11_B_sorted.bam,,Ag11_C_sorted.bam,,BGV344_sorted.bam,,Hg199_sorted.bam,,ND8_sorted.bam,,ND9_sorted.bam,,OPC304_sorted.bam,,P112_sorted.bam,,R0905_sorted.bam,,R37-15_sorted.bam,,R39-15_sorted.bam,,R41-15_sorted.bam,,R42-15_sorted.bam,,R45-15_sorted.bam,,R6-17-2_sorted.bam,,R6-17-3_sorted.bam --pop2 R68-17-C2_sorted.bam,,R68-17-C3_sorted.bam,,SVK1_sorted.bam,,SVK2_sorted.bam --thr 0.95
```

## lumpy: a general probabilistic framework for structural variant discovery

```bash
input_hap=/home/groups/harrisonlab/project_files/neonectria_ditissima/qc_dna/paired/N.ditissima
input_hap_assembly=../../../Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa

scripts=/home/sobczm/bin/popgen/snp
input=/data/scratch/gomeza/analysis/sv_calling/lumpy

##Warning!!!!. In some cases, the forward and reverse read files are corrupted (reads do not match in the two files)
#and bwa-mem will complain about it, and exit prematurely. For those samples, one needs to first fix the input reads with
#$scripts/sub_pairfq.sh  

cd $input/
for sample in $input_hap/*
do
reads_forward=$sample/F/*trim.fq.gz
#Sample name is the first part of the filename with reads, until the first underscore (_) encountered.
sname=$(echo $(basename "$reads_forward") | cut -d "_" -f1)
qsub $scripts/sub_bwa_mem.sh Illumina $sname $input_hap_assembly $reads_forward $reads_reverse
done

#Warning: the last step in the script - genotype calling - takes days.
qsub $scripts/sub_lumpy.sh Nditissima_struc_variants
```
