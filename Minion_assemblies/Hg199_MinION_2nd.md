# Neonectria_Reference_Genome_Assembly

======================================

This document details the commands used to assemble and annotate the Hg199 Neonectria genome and a 3rd version of the R0905 genome, using same tools and comparing between them.

Here, I repeat the genome assembly of the Reference Genome.

## Data extraction

R0905 PacBio

```bash
  	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  	RawDatDir=/home/harrir/projects/pacbio_test/n_dit/
  	mkdir -p raw_dna/pacbio/N.ditissima/R0905
  	cp -r $RawDatDir/D08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/E08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/F08_1 raw_dna/pacbio/N.ditissima/R0905/.
    OutDir=raw_dna/pacbio/N.ditissima/R0905/extracted
  	mkdir -p $OutDir
  	cat raw_dna/pacbio/N.ditissima/R0905/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
```

## Identify sequencing coverage

## Hg199

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

## R0905

For PacBio data

```bash
    for RawData in $(ls raw_dna/pacbio/*/*/extracted/*.fastq); do
        echo $RawData;
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
        GenomeSz=45
        OutDir=$(dirname $RawData)
        qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
    done
```

```bash
  for StrainDir in $(ls -d raw_dna/pacbio/*/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

Pacbio extracted coverage was: 91.44

For Miseq data:

```bash
    for RawData in $(ls qc_dna/paired/*/R0905/*/*q.gz); do
        echo $RawData;
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
        qsub $ProgDir/run_fastqc.sh $RawData;
        GenomeSz=45
        OutDir=$(dirname $RawData)
        qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
    done
```

```bash
    for StrainDir in $(ls -d qc_dna/paired/*/R0905); do
        Strain=$(basename $StrainDir)
        printf "$Strain\t"
        for File in $(ls $StrainDir/*/*.txt); do
            echo $(basename $File);
            cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
        done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
    done
```

Miseq coverage was:	134.1

## Identify sequencing coverage after porechop for the Hg199 assembly

cp /home/.../neonectria_ditissima/qc_dna/minion/N.ditisima/Hg199/Hg199_fastq_allfiles_trim.fastq.gz /data/scratch/gomeza/qc_dna/N.ditisima/Hg199

For Minion data:

```bash
  for RawData in $(ls qc_dna/minion/*/*/*q.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    GenomeSz=45
    OutDir=$(dirname $RawData)
    qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
  done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

Allfiles MinION coverage was: 94.14

## Read correction and assembly

Hg199 - Assembly using CANU

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

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_minion/N.ditissima/"$Strain"_40m
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction.sh $TrimReads 40m $Strain $OutDir
  done
```
```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_minion/N.ditissima/"$Strain"_50m
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction.sh $TrimReads 50m $Strain $OutDir
  done
```

```bash
  for CorrectedReads in $(ls assembly/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
    Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_minion/N.ditissima/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_assembly_only.sh $CorrectedReads 45m $Strain $OutDir
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    165 contigs

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

Hg199 - Assembly using CANU 1 step

```bash
	for Reads in $(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz); do
  	GenomeSz="45m"
  	Strain=$(echo $Reads | rev | cut -f2 -d '/' | rev)
  	Organism=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir=assembly/canu_minion_1step/$Organism/"$Strain"
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  	qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
  done    
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu_minion_1step/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    178 contigs

```bash
for Assembly in $(ls assembly/canu_minion_1step/N.ditissima/Hg199/Hg199_canu.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly_canu_1step
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

Hg199 - Assembly using SMARTdenovo

```bash
  for CorrectedReads in $(ls assembly/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
    Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/N.ditissima/"$Strain"
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
    qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    137 contigs

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

R09/05 - I will test different assembly methods. Canu in 2 steps and in 1 step.

R09/05 - Assembbly using CANU

```bash
  for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu_pacbio/N.ditissima/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_correction.sh $TrimReads 45m $Strain $OutDir
  done
```

```bash
  for CorrectedReads in $(ls assembly/canu_minion/N.d*/Hg199/*.trimmedReads.fasta.gz); do
    Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_pacbio/N.ditissima/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_assembly_only.sh $CorrectedReads 45m $Strain $OutDir
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/R0905.contigs.fasta); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    541 contigs. BAD, no Busco done.

R0905 - Assembly using CANU 1 step

```bash
	for Reads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
  	GenomeSz="45m"
  	Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir=assembly/canu2/$Organism/"$Strain"
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  	qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
  done    
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    50 contigs

```bash
for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

R09/05 - Assembbly using SMARTdenovo

```bash
  for CorrectedReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $CorrectedReads | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/N.ditissima/"$Strain"
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
    qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
  done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

    135 contigs

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

    Bad prediction

## Error correction using racon:

Hg199

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
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

```bash
for Assembly in $(ls assembly/canu_minion/N.ditissima/Hg199/Hg199.contigs.fasta); do
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
for Assembly in $(ls assembly/*/N.ditissima/Hg199/racon_10/*10.fasta); do
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
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_SMARTdenovo/racon_10
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls assembly/canu_minion/N.ditissima/Hg199/racon_10/*10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_canu_minion/racon_10
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

I will run busco in all of the racon iterations

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Iterations=$(echo $Assembly | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_canu_minion/racon_all/$Iterations
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

I will repeat racon with 20 interations

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz)
Iterations=20
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_20/*20.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_20/*20.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_SMARTdenovo/racon_20
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

R09/05

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/*/N.ditissima/R0905/Original_v3/racon_10/*10.fasta); do
  Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/R0905/racon_10/*10.fasta); do
  Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/*/N.ditissima/R0905/Original_v3/racon_10/*10.fasta); do
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_canu_minion/racon_10
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/R0905/racon_10/*10.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_SMARTdenovo/racon_10
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Hg199 Assembly correction using nanopolish

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_smartdenovo_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Event information would have been used from all of the runs, however MinKnow doesnt
# produce event-level information and therefore just the albacore data was used.
ReadsFq1=$(ls raw_dna/minion/Hg199_fastq_allfiles.fastq.gz)
#ReadsFq2=$(ls raw_dna/minion/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz)
#ReadsFq1=$(ls raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/fastq_runid_298a8dbc00c3db453901232f1ad01b11fd094980_pass.fastq.gz)
#ReadsFq2=$(ls raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz)
cat $ReadsFq1 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq
/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
Fast5Dir1=$(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/Nd_Hg199_20171203)
Fast5Dir2=$(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/Nd_Hg199_20171025)
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
```

```bash
Assembly=$(ls assembly/canu_minion/N.ditissima/Hg199/racon_10/*10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadDir=raw_dna/nanopolish/$Organism/$Strain
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
```

Split the assembly into 50Kb fragments and submit each to the cluster for nanopolish correction

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_smartdenovo_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)
NanoPolishDir=/home/gomeza/prog/nanopolish/scripts
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
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_smartdenovo_racon_round_10.fasta)
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

## R0905 Assembly correction using nanopolish

This part is not needed.

```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/R0905/racon_10/R0905_smartdenovo_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $Reads $OutDir/nanopolish
```

```bash
Assembly=$(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/racon_10/R0905_canu_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $Reads $OutDir/nanopolish
```

## Hg199 Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep 'Hg199'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/../pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

Contigs were renamed

```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/N.d*/Hg199/pilon/*.fasta | grep 'pilon_10')
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

    Filename    Complete    Duplicated    Fragmented    Missing    Total
    short_summary_Hg199_contigs_unmasked.txt	3670	15	25	30	3725 (MiSeq data)
    short_summary_Hg199_nanoplish_min_500bp_renamed.txt	3314	12	162	249	3725
    short_summary_pilon_10.txt	3586	17	20	119	3725
    short_summary_pilon_1.txt	3582	17	21	122	3725
    short_summary_pilon_2.txt	3584	17	21	120	3725
    short_summary_pilon_3.txt	3586	17	20	119	3725
    short_summary_pilon_4.txt	3586	17	20	119	3725
    short_summary_pilon_5.txt	3586	17	20	119	3725
    short_summary_pilon_6.txt	3586	17	20	119	3725
    short_summary_pilon_7.txt	3586	17	20	119	3725
    short_summary_pilon_8.txt	3586	17	20	119	3725
    short_summary_pilon_9.txt	3586	17	20	119	3725
    short_summary_pilon_min_500bp_renamed.txt	3586	17	20	119	3725

## Hybrid Assembly

### Spades Assembly

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

Quast and busco were run to assess the quality of hybrid assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades_minion/*/*/contigs.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
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
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
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
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_20k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_20k
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
AnchorLength=5000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_5k
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_5k
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
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly_merged
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

Hybrid assemblies are more fragmented with bigger number of busco genes. Minion are less fragmented but with less busco genes presnt.

I will run again quickmerge using N50 values for AnchorLength. Lowering this value may lead to more merging but may increase the probability of mis-joins.

```bash
for MinIONAssembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
Organism=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $MinIONAssembly | rev | cut -f3 -d '/' | rev)
HybridAssembly=$(ls assembly/spades_minion/$Organism/$Strain/filtered_contigs/contigs_min_500bp.fasta)
QuastReport=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/pilon/report.tsv)
N50=$(cat $QuastReport | grep 'N50' | cut -f2)
AnchorLength=$N50
#AnchorLength=600000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_600k
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_600k
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength

QuastReport=$(ls assembly/spades_minion/N.ditissima/Hg199/filtered_contigs/report.tsv)
N50=$(cat $QuastReport | grep 'N50' | cut -f2)
AnchorLength=$N50
#AnchorLength=500000
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_minion_500k
qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_hybrid_500k
qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
done
```

Quast and busco were run to assess the quality of hybrid assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*00k/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*00k/merged.fasta | grep 'Hg199'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_merged
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
AnchorLength=N50 does not improve the completeness of the genome. 





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
      OutDir=$(dirname $Assembly)/pilon
      Iterations=5
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
      qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```

Contigs were renamed

```bash
echo "" > tmp.txt
Assembly=$(ls assembly/merged_canu_spades/N.ditissima/Hg199_minion_5k/pilon/*.fasta | grep 'pilon_5')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/Hg199_minion_5k/pilon/pilon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/Hg199_minion_5k/pilon/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly2
#OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/*/assembly2/*/short_summary_*.txt | grep 'Hg199'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

## Isolate R09/05 Reference

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

## Repeat Masking

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

    Number of masked bases:
    3832437

## R0905 Assemblies were polished using Pilon

```bash
	for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
	Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905)
  TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
  TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
  Iterations=5
  OutDir=assembly/canu_pacbio/$Organism/$Strain/Original_v3/polished
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
  qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done
```

```bash
	for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/racon_10/R0905_canu_racon_round_10.fasta); do
	Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905)
  TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
  TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
  Iterations=5
  OutDir=assembly/canu_pacbio/$Organism/$Strain/Original_v3/racon_10/polished
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
  qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu_pacbio/N.di*/R0905/Original_v3/polished/*5.fasta); do
  Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done

for Assembly in $(ls assembly/canu_pacbio/N.di*/R0905/Original_v3/racon_10/polished/*5.fasta); do
  Strain=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/canu_pacbio/N.di*/R0905/Original_v3/polished/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls assembly/canu_pacbio/N.di*/R0905/Original_v3/racon_10/polished/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/assembly_afterracon
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/N*/Ref_Genomes/*/assembly/*/short_summary_*.txt | grep 'R0905'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

### Spades Assembly

```bash
  	for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $StrainPath
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/R0905)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    OutDir=assembly/spades_pacbio/$Organism/$Strain
    qsub -R y $ProgDir/sub_spades_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $OutDir 15
  	done
```

cat contigs.fasta | grep 'NODE' | wc -l
641

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/spades/*/*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

Contigs shorter than 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs_min_500bp
    FilterDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs_min_500bp/contigs_min_500bp.fasta
  done
```

cat contigs_min_500bp.fasta | grep 'NODE' | wc -l
364

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

## Merging pacbio and hybrid assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
  # for PacBioAssembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/polished/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    # HybridAssembly=$(ls assembly/spades_pacbio/$Organism/Fus2/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
    # OutDir=assembly/pacbio_test/$Organism/$Strain
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir
  done
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

This merged assembly was polished using Pilon

```bash
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2/merged.fasta); do
  # for Assembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    # IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    # OutDir=assembly/pacbio_test/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

\#Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  # for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
  for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/R0905/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    # OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Assembly stats were collected using quast

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # for Assembly in $(ls assembly/merged_canu_spades/*/*/filtered_contigs/R0905_contigs_renamed.fasta); do
  for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

\#Canu assembly contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/canu/N.ditissima/R0905/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

## Renaming assemblies - temporarily

Fus2 was temporarily renamed for preliminary analysis

```bash
  cp -r assembly/canu/F.oxysporum_fsp_cepae/Fus2 assembly/canu/F.oxysporum_fsp_cepae/Fus2_pacbio_test_canu
  cp -r assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2 assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_pacbio_test_merged
  cp -r assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged_richards
```

# Repeatmasking assemblies

```bash
  R0905_pacbio_merged=$(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
  # for Assembly in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  for Assembly in $(ls $R0905_pacbio_merged); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

```bash
  R0905_pacbio_canu=$(ls assembly/canu/*/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
  # for Assembly in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  for Assembly in $(ls $R0905_pacbio_canu); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

\*\* % bases masked by repeatmasker: 11.57% (bases masked:5290972 bp)

\*\* % bases masked by transposon psi: 10.29% (bases masked:4704506 bp)

Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done
```

repeat_masked/N.ditissima/R0905_pacbio_canu/filtered_contigs_repmask/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:5398957
