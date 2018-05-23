# Neonectria_Reference_Genome_Assembly

======================================

This document details the commands used to assemble and annotate the Hg199 Neonectria genome.

Here, I repeat the genome assembly of the Reference Genome.

## Identify sequencing coverage after porechop

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
Prefix=$Strain_smartdenovo
OutDir=assembly/SMARTdenovo/N.ditissima/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```
