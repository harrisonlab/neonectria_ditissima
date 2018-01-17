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
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stocks4'); do
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
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly/grep4)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/*.fasta); do
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
ScratchDir=/data/scratch/nanopore_tmp_data/Neonectria_temp
mkdir -p $ScratchDir
cp raw_dna/minion/N.ditissima/Hg199/pass/*.tar.gz /data/scratch/nanopore_tmp_data/Neonectria_temp

for Tar in $(ls $ScratchDir/*.tar.gz); do
  tar -zxvf $Tar -C $ScratchDir
done
```

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
Fast5Dir=$(ls -d raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/fast5_runid_298a8dbc00c3db453901232f1ad01b11fd094980_pass.tar.gz)
ReadDir=raw_dna/nanopolish/$Organism/$Strain
if [ -d $ReadDir ]; then
    echo "reads already extracted"
else
    echo "extracting reads"
    mkdir -p $ReadDir
    CurDir=$PWD
    cd $ReadDir
    nanopolish extract -r $Fast5Dir | gzip -cf > "$Strain"_reads.fa.gz
    cd $CurDir
fi

RawReads=$(ls $ReadDir/"$Strain"_reads.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $RawReads $OutDir/nanopolish
```

Split the assembly into 50Kb fragments and submit each to the cluster for nanopolish correction

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/Hg199_reads.fa.gz)
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
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/nanopolish/Stocks4_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/nanopolish/Stocks4_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stocks4'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

## Pilon error correction

Assemblies were polished using Pilon

```bash
    for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep 'Stocks4'); do
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
        TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
        TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
        OutDir=$(dirname $Assembly)/../pilon
        Iterations=5
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
        qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done
```

Contigs were renamed

```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/pilon/*.fasta | grep 'pilon_5')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/pilon/*.fasta ); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stocks4'); do
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
short_summary_Stocks4_SMARTdenovo.txt    365    0    400    2960    3725
short_summary_wtasm_racon_round_1.txt    1940    11    873    912    3725
short_summary_wtasm_racon_round_2.txt    2239    13    744    742    3725
short_summary_wtasm_racon_round_3.txt    2308    16    717    700    3725
short_summary_wtasm_racon_round_4.txt    2281    17    750    694    3725
short_summary_wtasm_racon_round_5.txt    2312    17    709    704    3725
short_summary_wtasm_racon_round_6.txt    2286    14    749    690    3725
short_summary_wtasm_racon_round_7.txt    2306    17    741    678    3725
short_summary_wtasm_racon_round_8.txt    2314    12    740    671    3725
short_summary_wtasm_racon_round_9.txt    2314    19    736    675    3725
short_summary_wtasm_racon_round_10.txt    2350    18    730    645    3725
short_summary_Stocks4_nanoplish_min_500bp_renamed.txt    3326    27    211    188    3725
short_summary_pilon_1.txt    3669    33    21    35    3725
short_summary_pilon_2.txt    3679    34    18    28    3725
short_summary_pilon_3.txt    3684    34    15    26    3725
short_summary_pilon_4.txt    3684    34    15    26    3725
short_summary_pilon_5.txt    3684    34    15    26    3725
short_summary_pilon_min_500bp_renamed.txt    3684    34    15    26    3725
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
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/hybrid_assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

## Merging MinION and hybrid assemblies

Note - the anchor length is the starting point for contigs to be merged - only contigs larger than these will be extended.

```bash
for MinIONAssembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
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
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_first_20k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
AnchorLength=200000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_first_200k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_first_200k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
done
```

Quast and busco were run to assess the quality of hybrid assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'Stocks4'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'Stocks4'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
# OutDir=gene_pred/busco/$Organism/$Strain/assembly
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

-->

# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta | grep 'Stocks4'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'Stocks4'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'Stocks4'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```
Number of masked bases:
8321901
```
