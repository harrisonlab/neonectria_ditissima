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
```
