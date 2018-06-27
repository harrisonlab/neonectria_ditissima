#Looks for structural variations between the genomes
##Structural variants can include: duplications, deletions, inversions & translocations. Uses read-pair configuration, split-reads & read-depth

###First, run BWA-mem to align reads to Reference genome assembly before SV calling

####Run BWA-mem

```bash
Reference=Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa
for Strain in Ag11_C BGV344 ND9 OPC304 P112 Ag06 Ag09_A Ag11_A R39-15 R42-15 R68-17 Ag11_B R41-15 R6-17-2 R6-17-3 Ag02 Ag05 ND8 R37-15 Ag04 R45-15 R0905_canu_2017_v2 Hg199; do
for Reads in $(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/qc_dna/paired/N.ditissima/*)
do
  Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    FRead=$Reads/F/*.fq.gz
    RRead=$Reads/R/*.fq.gz
    OutDir=alignment
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa.sh $Strain $Reference $FRead $RRead $OutDir
done
```

####Run svaba

```bash
Prefix=Nd_svaba
Reference=Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa
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
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/sv_calling
```

##### Create a cut-down VCF and filter it

```bash
cd $input

vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Pfrag_svaba_sv.svaba.indel.vcf SCRP245_sorted.bam ONT3_sorted.bam Nov77_sorted.bam Bc23_sorted.bam > Pfrag_svaba_sv.svaba.indel_cut.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Pfrag_svaba_sv.svaba.indel_cut.vcf  --max-missing 0.95 --recode --out Pfrag_svaba_sv.svaba.indel_cut_filtered
```

No sites removed by filtering

#### Ancestral variants

##### For UK2, set UK2 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK2.vcf --ply 2 --pop1 Bc16,,A4,,SCRP249,,SCRP324,,SCRP333 --pop2 Nov5,,Bc1,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK1 based analysis, set UK1 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK1.vcf --ply 2 --pop1 Bc1,,Nov5,,SCRP249,,SCRP324,,SCRP333 --pop2 A4,,Bc16,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK3 based analysis, set UK3 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK3.vcf --ply 2 --pop1 Nov9,,Nov27,,Nov71,,SCRP249,,SCRP324,,SCRP333 --pop2 A4,,Bc16,,Nov5,,Bc1 --thr 0.95
```

```
Nothing found
```

#### Just private variants, without addressing ancestral state

##### Create cut down VCF and filter it

```bash
vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Pfrag_svaba_sv.svaba.indel.vcf SCRP245_sorted.bam ONT3_sorted.bam Nov77_sorted.bam Bc23_sorted.bam SCRP249_sorted.bam SCRP324_sorted.bam SCRP333_sorted.bam > Pfrag_svaba_sv.svaba.indel_cut_UK123.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Pfrag_svaba_sv.svaba.indel_cut_UK123.vcf  --max-missing 0.95 --recode --out Pfrag_svaba_sv.svaba.indel_cut_UK123_filtered
```

##### For UK2, set UK2 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_UK123_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK2.vcf --ply 2 --pop1 Bc16,,A4 --pop2 Nov5,,Bc1,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK1 based analysis, set UK1 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_UK123_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK1.vcf --ply 2 --pop1 Bc1,,Nov5 --pop2 A4,,Bc16,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK3 based analysis, set UK3 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.indel_cut_UK123_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.indel_cut_filtered_fixed_UK3.vcf --ply 2 --pop1 Nov9,,Nov27,,Nov71 --pop2 A4,,Bc16,,Nov5,,Bc1 --thr 0.95
```

```
Nothing found
```

### Analysis of files produced by svaba

#### Analysis of sv file

This file contains larger indels

##### Set initial variables

```bash
scripts=/home/sobczm/bin/popgen/summary_stats
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/sv_calling
```

##### Create a cut-down VCF and filter it

```bash
cd $input

vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Pfrag_svaba_sv.svaba.sv.vcf SCRP245_sorted.bam ONT3_sorted.bam Nov77_sorted.bam Bc23_sorted.bam > Pfrag_svaba_sv.svaba.sv_cut.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Pfrag_svaba_sv.svaba.sv_cut.vcf  --max-missing 0.95 --recode --out Pfrag_svaba_sv.svaba.sv_cut_filtered
```

#### Ancestral variants

##### For UK2, set UK2 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK2.vcf --ply 2 --pop1 Bc16,,A4,,SCRP249,,SCRP324,,SCRP333 --pop2 Nov5,,Bc1,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK1 based analysis, set UK1 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK1.vcf --ply 2 --pop1 Bc1,,Nov5,,SCRP249,,SCRP324,,SCRP333 --pop2 A4,,Bc16,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK3 based analysis, set UK3 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK3.vcf --ply 2 --pop1 Nov9,,Nov27,,Nov71,,SCRP249,,SCRP324,,SCRP333 --pop2 A4,,Bc16,,Nov5,,Bc1 --thr 0.95
```

```
Nothing found
```

#### Just private variants, without assessing ancestral state

##### Create cut down VCF and filter it

```bash
vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Pfrag_svaba_sv.svaba.sv.vcf SCRP245_sorted.bam ONT3_sorted.bam Nov77_sorted.bam Bc23_sorted.bam SCRP249_sorted.bam SCRP324_sorted.bam SCRP333_sorted.bam > Pfrag_svaba_sv.svaba.sv_cut_UK123.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Pfrag_svaba_sv.svaba.sv_cut_UK123.vcf  --max-missing 0.95 --recode --out Pfrag_svaba_sv.svaba.sv_cut_UK123_filtered
```

##### For UK2, set UK2 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK2.vcf --ply 2 --pop1 Bc16,,A4 --pop2 Nov5,,Bc1,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK1 based analysis, set UK1 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK1.vcf --ply 2 --pop1 Bc1,,Nov5 --pop2 A4,,Bc16,,Nov9,,Nov27,,Nov71 --thr 0.95
```

```
Nothing found
```

##### UK3 based analysis, set UK3 isolates and P. rubi isolates as pop1

```bash
python $scripts/vcf_find_difference_pop.py --vcf Pfrag_svaba_sv.svaba.sv_cut_filtered.recode.vcf --out Pfrag_svaba_sv.svaba.sv_cut_filtered_fixed_UK3.vcf --ply 2 --pop1 Nov9,,Nov27,,Nov71 --pop2 A4,,Bc16,,Nov5,,Bc1 --thr 0.95
```

```
Nothing found
