# Neonectria Reference Genome Assembly

======================================

This document details the commands used to assemble and annotate a reference genome of the Hg199 isolate.

Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/N.ditissima or /data/scratch/gomeza

## Data management

### First run

First MinION run date is 20170717. Data stored in 20170717_1504_Neonectria_Hg199 and 20170717_1452_Neonectria_Hg199

Albacore 1.1.1\. was used to rerun basecalling. Data stored in data/seq_data/minion/2017/20170717_recalled_Neonectria_Hg199

New version of Albacore is available, v.2.0.2.

```bash
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.0.2-cp34-cp34m-manylinux1_x86_64.whl
pip3 install --user ont_albacore-2.0.2-cp34-cp34m-manylinux1_x86_64.whl
```

Running albacore 2.0.2.

```bash
screen -S gomeza

~/.local/bin/read_fast5_basecaller.py \
--flowcell FLO-MIN106
--kit SQK-LSK108 \
--input /data/seq_data/minion/2017/20170717_1504_Neonectria_Hg199 \
--recursive \
--worker_threads 12 \
--save_path /home/nanopore/Neonectria_Hg199_2_0_2 \
--output_format fastq,fast5 \
--reads_per_fastq_batch 4000
```

```bash
cat  *.fastq | gzip > fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz
cat  *.fastq | gzip > fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_fail.fastq.gz
tar -czf fast5_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.tar.gz [0-9]*
tar -czf fast5_runid_5832f037a56936787d17e66d1e3b8ac05572199f_fail.tar.gz [0-9]*

rm -rf *.fastq [0-9]*
```

Rob's command to transfer the date from nanopore node to head node.
Permission needed scp -r ./Neonectria_Hg199_2_0_2 miseq_data@192.168.1.200:/data/seq_data/minion/2017/20171025_Neonectria_Hg199_2_0_2

### Second run

Second MinION run date 20171203. Data stored in 20171203_Hg199.

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

## Building of directory structure

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

## Identify sequencing coverage

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

## Splitting reads and trimming adapters using porechop

```bash
    RawReads=raw_dna/minion/N.ditissima/Hg199/*allfiles.fastq.gz
    OutDir=qc_dna/minion/N.ditissima/Hg199_2ndround
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
```

## Identify sequencing coverage after porechop for the Hg199 assembly


For Minion data:

```bash
  for RawData in $(ls qc_dna/minion/*/*round/*q.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    GenomeSz=45
    OutDir=$(dirname $RawData)
    qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
  done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/*round); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

Allfiles MinION coverage was: 94.14

# CANU

## Assembly using CANU. Read correction and assembly

Read correction using Canu

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=Hg199
    OutDir=assembly_vAG/canu_2steps/canu_minion/N.ditissima/"$Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction_nanopore.sh $TrimReads 45m $Strain $OutDir
  done
```
```bash
  for CorrectedReads in $(ls assembly_vAG/canu_2steps/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
    Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly_vAG/canu_minion/N.ditissima/"$Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_assembly_only_nanopore.sh $CorrectedReads 45m $Strain $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu_2steps/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/canu_2steps/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
for Assembly in $(ls assembly_vAG/canu_2steps/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Error correction using racon:

```bash
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/racon_10/*10.fasta); do
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/canu_2step/canu_minion/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```
## Assembly using CANU. Default run.

```bash
	for Reads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/*allfiles_trim.fastq.gz); do
  	GenomeSz="45m"
  	Strain=Hg199
  	Organism=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir=assembly_vAG/canu/$Organism/"$Strain"
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  	qsub $ProgDir/submit_canu_minion.sh $Reads $GenomeSz $Prefix $OutDir
  done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Error correction using racon: (no need to do this with canu)

```bash
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Pilon error correction (this should be done without racon)

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

# SMARTdenovo

## Assembly using SMARTdenovo

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=Hg199
    Prefix="$Strain"_smartdenovo
    OutDir=assembly_vAG/SMARTdenovo/N.ditissima/"$Strain"
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
    qsub $ProgDir/sub_SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Error correction using racon:

```bash
  for Assembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.txt
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/Hg199_racon10_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/SM*/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

# Miniasm

## Assembly of uncorrected reads. Racon required after.

```bash
# If reads have same names or same part splitted by space, fix them using rename.sh from bbtools.
# This can only be done in blacklace11, since it needs a specific library for blasr.
ssh blacklace11.blacklace
/home/gomeza/prog/bbtools/bbmap/rename.sh in=Hg199_fastq_allfiles_trim.fastq.gz out=trimmed_renamed.fasta prefix=Hg199

# Fast all-against-all overlap of raw reads
# Overlap for MinION reads (or use "-x ava-pb" for Pacbio read overlapping)
/home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > Hg199_fastq_allfiles.paf.gz

# Concatenate pieces of read sequences to generate the final sequences.
# Thus the per-base error rate is similar to the raw input reads. Make sure you correct your reads.
# Layout
/home/gomeza/prog/miniasm/miniasm -f trimmed_renamed.fasta Hg199_fastq_allfiles.paf.gz > reads.gfa

# Convert gfa file to fasta file.
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > Hg199_miniasm.fa
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/Hg199_miniasm.fa); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/Hg199_miniasm.fa); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

#### Error correction using racon:

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/Hg199_miniasm.fa); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.txt
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/Hg199_racon10_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
mv final assemblies to assembly_vAG folder

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.txt
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/pilon_10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```
Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/*10.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

# Merging MinION and hybrid assemblies

Note - the anchor length is the starting point for contigs to be merged - only contigs larger than these will be extended.

Note- The Python wrapper was written for MUMmer version 3.x command line. Hash out version 4 from profile while this is run.

### Canu

```bash
  for MinIONAssembly in $(ls assembly_vAG/canu/N.ditissima/Hg199/pilon/pilon_10.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_minion/$Organism/$Strain/filtered_contigs/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=700000
    OutDir=assembly_vAG/merged_assemblies/canu_spades/$Organism/"$Strain"_minion_700k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/canu_spades/$Organism/"$Strain"_hybrid_700k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

### Miniasm

```bash
  for MinIONAssembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/Hg199_pilon10_renamed.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_minion/$Organism/$Strain/filtered_contigs/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=600000
    OutDir=assembly_vAG/merged_assemblies/miniasm_spades/$Organism/"$Strain"_minion_600k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/miniasm_spades/$Organism/"$Strain"_hybrid_600k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

###SMARTdenovo

```bash
  for MinIONAssembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/Hg199/racon_10/pilon/pilon_10.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_minion/$Organism/$Strain/filtered_contigs/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=700000
    OutDir=assembly_vAG/merged_assemblies/SMARTdenovo_spades/$Organism/"$Strain"_minion_700k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/SMARTdenovo_spades/$Organism/"$Strain"_hybrid_700k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

## Quast and busco were run to assess the quality of hybrid assemblies:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/*/merged.fasta); do
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
Quast and busco were run to assess the quality of hybrid assemblies:

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/*/pilon/*_10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/*/pilon/*_10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Repeat masking. Best assembly for further analysis.

Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/pilon_10.fasta); do
    OutDir=$(dirname $Assembly)/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_pilon10_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv  
```
```bash
  for Assembly in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/Hg199_pilon10_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)/repeat_masked
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```
```bash
for RepDir in $(ls -d assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked); do
    Strain=$(echo $RepDir | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f6 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/Hg199/pilon/pilon_10.fasta); do
    OutDir=$(dirname $Assembly)/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_pilon10_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv  
```
```bash
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/Hg199/pilon/filtered_contigs/Hg199_pilon10_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)/repeat_masked
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```
```bash
for RepDir in $(ls -d assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked); do
    Strain=$(echo $RepDir | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f6 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly_vAG/merged_assemblies/miniasm_spades/N.ditissima/Hg199_minion_600k/pilon/pilon_10.fasta); do
    OutDir=$(dirname $Assembly)/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_pilon10_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv  
```
```bash
  for Assembly in $(ls assembly_vAG/merged_assemblies/miniasm_spades/N.ditissima/Hg199_minion_600k/pilon/filtered_contigs/Hg199_pilon10_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)/repeat_masked
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```
```bash
for RepDir in $(ls -d assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked); do
    Strain=$(echo $RepDir | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f6 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```
