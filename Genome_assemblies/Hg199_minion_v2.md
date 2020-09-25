 # Neonectria ditissima Hg199 isolate. 
======================================

## Genome assembly

Genome assembly of minion reads using new assemblers

### Flye

```bash
    # This work imply testing new tools, so it will be performed in the gomez_WD folder
    for TrimReads in $(ls ../qc_dna/minion/N.ditissima/Hg199/Hg199_fastq_allfiles_trim.fastq.gz); do
        Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
        Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) 
        Prefix="$Strain"_flye
        TypeSeq=nanoraw
        OutDir=assembly/flye/$Organism/$Strain
        mkdir -p $OutDir
        Size=45m
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq
    done
  done
```
```bash
    for Assembly in $(ls assembly/flye/N.ditissima/Hg199/assembly.fasta); do
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
        sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done
```
```bash
    # This work imply testing new tools, so it will be performed in the gomez_WD folder
    for TrimReads in $(ls ../qc_dna/minion/N.ditissima/Hg199/Hg199_fastq_allfiles_trim.fastq.gz); do
        Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
        Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) 
        Prefix="$Strain"_flye
        TypeSeq=nanocorrected
        OutDir=assembly/flye_corrected/$Organism/$Strain
        mkdir -p $OutDir
        Size=45m
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq
    done
```
```bash
    for Assembly in $(ls assembly/flye_corrected/N.ditissima/Hg199/assembly.fasta); do
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
        sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done
```

## Medaka 

```bash
conda activate medaka
  for Assembly in $(ls assembly/fly*/N.ditissima/Hg199/assembly.fasta); do
    ReadsFq=$(ls ../qc_dna/minion/N.ditissima/Hg199/Hg199_fastq_allfiles_trim.fastq.gz)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done
```

## Racon

Consensus module for raw de novo DNA assembly

```bash
for Assembly in $(ls assembly/fly*/N.ditissima/Hg199/assembly.fasta); do
ReadsFq=$(ls ../qc_dna/minion/N.ditissima/Hg199/Hg199_fastq_allfiles_trim.fastq.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
# You might rename your contigs at this point using remove_contaminants.py
```
