# Neonectria_Reference_Genome_Assembly

==========

This document details the commands used to assemble and annotate the Hg199 Neonectria genome.

# Data management

## First MinION run date is 20170717. Data stored in 20170717_1504_Neonectria_Hg199 and 20170717_1452_Neonectria_Hg199

Albacore 1.1.1\. was used to rerun basecalling. Data stored in data/seq_data/minion/2017/20170717_recalled_Neonectria_Hg199

New version of Albacore available. v.2.0.2.

```bash
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.0.2-cp34-cp34m-manylinux1_x86_64.whl
pip3 install --user ont_albacore-2.0.2-cp34-cp34m-manylinux1_x86_64.whl
```

Running albacore 2.0.2.

```bash
screen -S gomeza

~/.local/bin/read_fast5_basecaller.py --flowcell FLO-MIN106 --kit SQK-LSK108 --input /data/seq_data/minion/2017/20170717_1504_Neonectria_Hg199 --recursive --worker_threads 12 --save_path /home/nanopore/Neonectria_Hg199_2_0_2 --output_format fastq,fast5 --reads_per_fastq_batch 4000
```

```bash
cat  *.fastq | gzip > fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz
cat  *.fastq | gzip > fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_fail.fastq.gz
tar -czf fast5_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.tar.gz [0-9]*
tar -czf fast5_runid_5832f037a56936787d17e66d1e3b8ac05572199f_fail.tar.gz [0-9]*

rm -rf *.fastq [0-9]*
```

Rob command to transfer the date from nanopore node to head node. Permission needed scp -r ./Neonectria_Hg199_2_0_2 miseq_data@192.168.1.200:/data/seq_data/minion/2017/20171025_Neonectria_Hg199_2_0_2

## Second MinION run date 20171203. Data stored in 20171203_Hg199.

Running albacore 2.0.2.

```bash
screen -S gomeza

~/.local/bin/read_fast5_basecaller.py \
--flowcell FLO-MIN106 \
--kit SQK-LSK108 \
--input /data/seq_data/minion/2017/20171203_Hg199/Hg199/GA50000/reads/ \
--recursive \
--worker_threads 12 \
--save_path /home/nanopore/20171203_Hg199_2_0_2 \
--output_format fastq,fast5 \
--reads_per_fastq_batch 4000
```

Data stored in /data/seq_data/minion/2017/20171203_Hg199

# Building of directory structure

```bash
  screen -a
  RawDatDir=/data/seq_data/minion/2017/20171025_Neonectria_Hg199_2_0_2/workspace
  Organism=N.ditissima
  Strain=Hg199
  Date=25-10-17
  mkdir -p raw_dna/minion/$Organism/$Strain/$Date
  cat /data/seq_data/minion/2017/20170717_recalled_Neonectria_Hg199/workspace/*fastq | gzip -cf > raw_dna/minion/$Organism/$Strain/$Date/"$Strain"_"$Date".fastq.gz
```

```bash
  RawDatDir=/data/seq_data/minion/2017/20170717_recalled_Neonectria_Hg199/workspace
  Organism=N.ditissima
  Strain=Hg199
  Date=23-10-17
  for Fast5Dir in $(ls -d $RawDatDir/*); do
    poretools fastq $Fast5Dir | gzip -cf
  done > raw_dna/minion/$Organism/$Strain/"$Strain"_"$Date"_fast5.fastq.gz
```

```bash
screen -a
  RawDatDir=/data/seq_data/minion/2017/20171203_Hg199/Hg199/GA50000
  Organism=N.ditissima
  Strain=Hg199
  Date=03-12-17
  mkdir -p raw_dna/minion/$Organism/$Strain/$Date/
  for Fast5Dir in $(ls -d $RawDatDir/*.fastq); do
      poretools fastq $Fast5Dir | gzip -cf
    done > raw_dna/minion/$Organism/$Strain/$Date/"$Strain"_"$Date".fastq.gz
```

Both, 2.0.2 version recalled runs were merged together

```bash
cat raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/fastq_runid_298a8dbc00c3db453901232f1ad01b11fd094980_pass.fastq.gz > Hg199_fastq_allfiles.fastq.gz
```

The following processes were applied to Neonectria genomes prior to analysis: Data qc Genome assembly Repeatmasking Gene prediction Functional annotation

Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/N.ditissima

# Identify sequencing coverage

For Minion data:

```bash
    for RawData in $(ls raw_dna/minion/*/*/*q.gz); do
        echo $RawData;
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
        GenomeSz=45
        OutDir=$(dirname $RawData)
        qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
    done
```

```bash
    for RawData in $(ls raw_dna/minion/*/*/03-12-17/*/*/*q.gz); do
        echo $RawData;
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
        GenomeSz=45
        OutDir=$(dirname $RawData)
        qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
    done
```

```bash
  for StrainDir in $(ls -d raw_dna/minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

Allfiles MinION coverage was: 94.46

03-12-17 MinION coverage was: 71.39

For Miseq data:

```bash
    for RawData in $(ls qc_dna/paired/*/Hg199/*/*q.gz); do
        echo $RawData;
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
        qsub $ProgDir/run_fastqc.sh $RawData;
        GenomeSz=45
        OutDir=$(dirname $RawData)
        qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
    done
```

```bash
    for StrainDir in $(ls -d qc_dna/paired/*/Hg199); do
        Strain=$(basename $StrainDir)
        printf "$Strain\t"
        for File in $(ls $StrainDir/*/*.txt); do
            echo $(basename $File);
            cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
        done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
    done
```

Miseq coverage was: 89.92

## Assembly

### Assembly will with full trimming of reads:

Splitting reads and trimming adapters using porechop

```bash
    RawReads=raw_dna/minion/N.ditissima/Hg199/*allfiles.fastq.gz
    OutDir=qc_dna/minion/N.ditissima/Hg199
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
```

Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu_minion/N.ditissima/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 45m $Strain $OutDir
done
```

Assembly using Canu

```bash
for CorrectedReads in $(ls assembly/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu_minion2/N.ditissima/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_assembly_only.sh $CorrectedReads 45m $Strain $OutDir
done
```

Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix=$Strain
OutDir=assembly/SMARTdenovo/N.ditissima/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

Quast for the SMARTdenovo assembly:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
# printf "Organism\tStrain\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/*/assembly/*/short_summary_*.txt | grep 'Hg199'); do
# echo $File;
# Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
# Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
# printf "$Organism\t$Strain\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/*10.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/*10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly/racon_10
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

# Assembly correction using nanopolish

```bash
ScratchDir=/data/scratch/gomeza
mkdir -p $ScratchDir
cp raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/*.tar.gz $ScratchDir/.
cp raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/*.tar.gz $ScratchDir/.
for Tar in $(ls $ScratchDir/*.tar.gz); do
  tar -zxvf $Tar -C $ScratchDir
done
```

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Event information would have been used from all of the runs, however MinKnow doesnt
# produce event-level information and therefore just the albacore data was used.
ReadsFq1=$(ls raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/fastq_runid_298a8dbc00c3db453901232f1ad01b11fd094980_pass.fastq.gz)
ReadsFq2=$(ls raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz)
cat $ReadsFq1 $ReadsFq2 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq
/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
Fast5Dir1=$(ls -d /data/scratch/gomeza/Nd_Hg199_20171203)
Fast5Dir2=$(ls -d /data/scratch/gomeza/Nd_Hg199_20171025)
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
#if [ -d $ReadDir ]; then
#    echo "reads already extracted"
#else
#    echo "extracting reads"
#    mkdir -p $ReadDir
#    CurDir=$PWD
#    cd $ReadDir
#    nanopolish extract -r $Fast5Dir1 $Fast5Dir2 | gzip -cf > "$Strain"_reads.fa.gz
#    cd $CurDir
#fi
#RawReads=$(ls $ReadDir/"$Strain"_reads.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
done
```

Split the assembly into 50Kb fragments and submit each to the cluster for nanopolish correction

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt
Ploidy=1
echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
```

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py assembly/SMARTdenovo/$Organism/$Strain/racon_10/*/*.fa > $OutDir/"$Strain"_nanoplish.fa
# for File in $(ls assembly/SMARTdenovo/F.*/*/racon2/*/*.fa | grep ":*-"); do
#     # echo $File
#     cat $File >> $OutDir/"$Strain"_nanoplish.fa
# done
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/nanopolish/Hg199_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/nanopolish/Hg199_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/*/assembly/*/short_summary_*.txt | grep 'Hg199'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

short_summary_Hg199_nanoplish_min_500bp_renamed.txt	3317	11	160	248	3725


## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
    for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep 'Hg199'); do
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
        TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
        TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
        OutDir=$(dirname $Assembly)/../pilon
        Iterations=5
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
        qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done
```

Contigs were renamed

```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/N.d*/Hg199/pilon/*.fasta | grep 'pilon_5')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.di*/Hg199/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.di*/Hg199/pilon/*.fasta ); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/*/assembly/*/short_summary_*.txt | grep 'Hg199'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
Filename    Complete    Duplicated    Fragmented    Missing    Total
short_summary_contigs_min_500bp.txt	3673	15	24	28	3725
short_summary_Hg199.dmo.lay.txt	408	0	438	2879	3725
short_summary_Hg199.dmo.lay.txt	789	1	661	2275	3725
short_summary_Hg199_nanoplish_min_500bp_renamed.txt	3317	11	160	248	3725
short_summary_pilon_1.txt	3581	17	22	122	3725
short_summary_pilon_2.txt	3583	17	22	120	3725
short_summary_pilon_3.txt	3586	17	20	119	3725
short_summary_pilon_4.txt	3586	17	20	119	3725
short_summary_pilon_5.txt	3586	17	20	119	3725
short_summary_pilon_min_500bp_renamed.txt	3586	17	20	119	3725
```

# Hybrid Assembly

## Spades Assembly

```bash
for TrimReads in $(ls qc_dna/minion/*/Hg199/*allfiles_trim.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=assembly/spades_minion/$Organism/"$Strain"
echo $TrimF1_Read
echo $TrimR1_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
qsub $ProgDir/sub_spades_minion.sh $TrimReads $TrimF1_Read $TrimR1_Read $OutDir
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```

Quast and busco were run to assess the quality of hybrid assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/spades_minion/*/Hg199/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly_spades
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Merging MinION and hybrid assemblies

Note - the anchor length is the starting point for contigs to be merged - only contigs larger than these will be extended.

```bash
for MinIONAssembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
Organism=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $MinIONAssembly | rev | cut -f3 -d '/' | rev)
HybridAssembly=$(ls assembly/spades_minion/$Organism/$Strain/filtered_contigs/contigs_min_500bp.fasta)
# QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
# N50=$(cat $QuastReport | grep 'N50' | cut -f2)
# AnchorLength=$N50

AnchorLength=20000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_first_20k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid__first_20k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
AnchorLength=5000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_5k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_5k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
done
```

Quast and busco were run to assess the quality of hybrid assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*_minion_5k/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly2
#OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
short_summary_merged.txt        3563    107     23      139     3725

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*_minion_first_20k/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly2
#OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
short_summary_merged.txt        2435    29      478     812     3725

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/Hg199/assembly/*/short_summary_*.txt | grep 'Hg199'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Minion_5K genome
3563 Complete BUSCOs (C)
INFO    3456 Complete and single-copy BUSCOs (S)
INFO    107 Complete and duplicated BUSCOs (D)
INFO    23 Fragmented BUSCOs (F)
INFO    139 Missing BUSCOs (M)

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
  for Assembly in $(ls assembly/merged_canu_spades/*/*_minion_5k/merged.fasta | grep 'Hg199'); do
      Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
      IlluminaDir=$(ls -d qc_dna/paired/*/Hg199)
      TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
      TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
      OutDir=$(dirname $Assembly)/../pilon
      Iterations=5
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
      qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```

Contigs were renamed

```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/N.d*/Hg199/pilon/*.fasta | grep 'pilon_5')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

# Isolate R09/05 Reference

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/R0905/*_contigs/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Contigs=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain - $Contigs"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/$Contigs
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/R0905/*_contigs/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Filename	Complete	Duplicated	Fragmented	Missing	Total
short_summary_R0905_canu_contigs_modified.txt	3510	93	26	189	3725
short_summary_R0905_contigs_renamed.txt	3510	94	26	189	3725

# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*_minion_5k/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/Ref_Genomes/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/Ref_Genomes/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```
Number of masked bases:
3832437
```

## Gene Prediction

Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

```bash

for Assembly in $(ls repeat_masked/N.ditissima/Ref_Genomes/Hg199*/*/*unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


short_summary_Hg199_contigs_unmasked.txt
				3670    Complete BUSCOs (C)
        3655    Complete and single-copy BUSCOs (S)
        15      Complete and duplicated BUSCOs (D)
        25      Fragmented BUSCOs (F)
        30      Missing BUSCOs (M)
        3725    Total BUSCO groups searched
short_summary_contigs_min_500bp.txt Hg199
				3673    Complete BUSCOs (C)
        3658    Complete and single-copy BUSCOs (S)
        15      Complete and duplicated BUSCOs (D)
        24      Fragmented BUSCOs (F)
        28      Missing BUSCOs (M)
        3725    Total BUSCO groups searched
short_summary_R0905_contigs_unmasked.txt R0905_merged_assembly
				3510    Complete BUSCOs (C)
        3417    Complete and single-copy BUSCOs (S)
        93      Complete and duplicated BUSCOs (D)
        26      Fragmented BUSCOs (F)
        189     Missing BUSCOs (M)
        3725    Total BUSCO groups searched
short_summary_R0905_contigs_unmasked.txt R0905 R0905_merged_2017
				3651    Complete BUSCOs (C)
        3635    Complete and single-copy BUSCOs (S)
        16      Complete and duplicated BUSCOs (D)
        20      Fragmented BUSCOs (F)
        54      Missing BUSCOs (M)
        3725    Total BUSCO groups searched
short_summary_R0905_contigs_unmasked.txt R0905_canu_2017_v2
				3668    Complete BUSCOs (C)
        3652    Complete and single-copy BUSCOs (S)
        16      Complete and duplicated BUSCOs (D)
        21      Fragmented BUSCOs (F)
        36      Missing BUSCOs (M)
        3725    Total BUSCO groups searched


# Gene model training

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
	for Assembly in $(ls repeat_masked/*/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
	  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
	  Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
	  echo "$Organism - $Strain"
	  for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199); do
	    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
	    echo "$Timepoint"
	    FileF=$(ls $RNADir/F/*_trim.fq.gz)
	    FileR=$(ls $RNADir/R/*_trim.fq.gz)
	    OutDir=alignment/$Organism/$Strain/$Timepoint
	    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
	    qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
done
```

Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

screen -a

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
for Assembly in $(ls repeat_masked/*/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
# AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
OutDir=gene_pred/cufflinks/$Organism/$Strain/vs_Hg199reads/concatenated_prelim
echo "$Organism - $Strain"
mkdir -p $OutDir
cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
done
```
I have aligned every isolate with the Hg199 RNA reads.

Output from stdout included:


The Estimated Mean: 219.68 allowed calculation of of the mean insert gap to be
-140bp 182-(180*2) where 180? was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
for Assembly in $(ls repeat_masked/N*/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199 | grep -v -e '_rep'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=alignment/$Organism/$Strain/$Timepoint
InsertGap='-140'
InsertStdDev='43'
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
done
printf "\n"
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
done
done
done
```

#### Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
  for Assembly in $(ls repeat_masked/N.ditissima/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		#for Assembly in $(ls repeat_masked/N.ditissima/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
    done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/$Strain
		#OutDir=gene_pred/braker/$Organism/"$Strain"_braker_Nov2017
    AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker_Jan2018
		#GeneModelName="$Organism"_"$Strain"_braker_Nov2017
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_Jan2018
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
  for File in $(ls gene_pred/braker/N.*/*/*/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```

## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
    for Assembly in $(ls repeat_masked/N*/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
    done
	done
```

Secondly, genes were predicted using CodingQuary:

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
    for Assembly in $(ls repeat_masked/N*/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
    done
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for Strain in Ag02 Ag05 ND8 R37-15; do
for BrakerGff in $(ls gene_pred/braker/$Organism/$Strain/*/augustus.gff3); do
#Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/codingquary/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

done
done
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/codingquary/N.*/*/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
