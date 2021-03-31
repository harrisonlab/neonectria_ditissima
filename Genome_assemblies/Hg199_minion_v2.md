 # Neonectria ditissima Hg199 isolate. 
======================================

## Genome assembly

Genome assembly of minion reads using new assemblers

### Flye

```bash
    # This work imply testing new tools, so it will be performed in the gomez_WD folder
    for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/Hg199_fastq_allfiles_trim.fastq.gz); do
        Organism=N.ditissima
        Strain=Hg199
        Prefix="$Strain"_flye
        TypeSeq=nanoraw
        OutDir=assembly_VP/flye/$Organism/$Strain
        mkdir -p $OutDir
        Size=45m
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq
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
    for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
        Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev)
        Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) 
        Prefix="$Strain"_flye
        TypeSeq=pacbioraw
        OutDir=assembly_VP/flye/$Organism/$Strain
        mkdir -p $OutDir
        Size=45m
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq
    done
```


```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
for Assembly in $(ls assembly_VP/flye/*/*/assembly.fasta); do
OutDir=$(dirname $Assembly)
sbatch $ProgDir/quast.sh $Assembly $OutDir
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

### SMARTdenovo

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/Hg199_fastq_allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=Hg199
    Prefix="$Strain"_smartdenovo
    OutDir=assembly_VP/SMARTdenovo/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
```

```bash
  for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Prefix="$Strain"_smartdenovo
    OutDir=assembly_VP/SMARTdenovo/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
```

### Miniasm


```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/Hg199_fastq_allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) 
    Strain=Hg199
    Prefix="$Strain"_miniasm
    OutDir=assembly_VP/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch -p himem $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done
```

```bash
  for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev) 
    Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Prefix="$Strain"_miniasm
    OutDir=assembly_VP/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch -p himem $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done
```

### Canu


```bash
for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/Hg199_fastq_allfiles_trim.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) 
Strain=Hg199
Size=45m
Prefix="$Strain"_canu
type=nanopore
OutDir=assembly_VP/canu/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch -p long $ProgDir/canu.sh $TrimReads $Size $Prefix $type $OutDir
done
```

```bash
for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev) 
Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Size=45m
Prefix="$Strain"_canu
type=pacbio
OutDir=assembly_VP/canu/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
sbatch -p himem $ProgDir/canu.sh $TrimReads $Size $Prefix $type $OutDir
done
```




#### Error correction using racon:

```bash
  for Assembly in $(ls assembly_VP/miniasm/N.ditissima/Hg199/Hg199_miniasm.fa); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/miniasm/N.ditissima/R0905/*_miniasm.fa); do
    Strain=R0905
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon2.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/Hg199/Hg199_smartdenovo.dmo.lay.utg); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
    Strain=R0905
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon2.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_VP/flye/N.ditissima/Hg199/assembly.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/flye/N.ditissima/R0905/assembly.fasta); do
    Strain=R0905
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon2.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```


## Quast and BUSCO 

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/*/racon_10/*_racon_round_10.fasta); do
  OutDir=$(dirname $Assembly)
  sbatch $ProgDir/quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/*/racon_10/*_racon_round_10.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/canu/N.ditissima/*/*_canu.contigs.fasta); do
  OutDir=$(dirname $Assembly)
  sbatch $ProgDir/quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/canu/N.ditissima/*/*_canu.contigs.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch -p long $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```

## Medaka post racon

```bash
conda activate medaka
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/Hg199/racon_10/*_racon_round_10.fasta); do
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done
```

```bash
conda activate medaka
  for Assembly in $(ls assembly_VP/SMARTdenovo/N.ditissima/R0905/racon_10/*_racon_round_10.fasta); do
    ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done
```
```bash
conda activate medaka
  for Assembly in $(ls assembly_VP/canu/N.ditissima/Hg199/*_canu.contigs.fasta); do
    ReadsFq=$(ls qc_dna/minion/N.ditissima/*round/*allfiles_trim.fastq.gz)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done

  for Assembly in $(ls assembly_VP/canu/N.ditissima/R0905/*_canu.contigs.fasta); do
    ReadsFq=$(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq)
    OutDir=$(dirname $Assembly)/medaka
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
  done
```

## Quast and Busco

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/*/N.ditissima/*/racon_10/medaka/medaka/consensus.fasta); do
      OutDir=$(dirname $Assembly)
      sbatch $ProgDir/quast.sh $Assembly $OutDir
  done

  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/canu/N.ditissima/*/medaka/medaka/consensus.fasta); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/*/N.ditissima/*/racon_10/medaka/medaka/consensus.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done

   for Assembly in $(ls assembly_VP/canu/N.ditissima/*/medaka/medaka/consensus.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```

## Pilon

```bash
  for Assembly in $(ls assembly_VP/*/N.ditissima/Hg199/racon_10/medaka/medaka/consensus.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done

  for Assembly in $(ls assembly_VP/canu/N.ditissima/Hg199/medaka/medaka/consensus.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```


```bash
  for Assembly in $(ls assembly_VP/*/N.ditissima/R0905/racon_10/medaka/medaka/consensus.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done

  for Assembly in $(ls assembly_VP/canu/N.ditissima/R0905/medaka/medaka/consensus.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```

## Quast and Busco

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/*/N.ditissima/*/racon_10/medaka/medaka/pilon/pilon_10.fasta); do
      OutDir=$(dirname $Assembly)
      sbatch $ProgDir/quast.sh $Assembly $OutDir
  done

  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/canu/N.ditissima/*/medaka/medaka/pilon/pilon_10.fasta); do
    OutDir=$(dirname $Assembly)
    sbatch $ProgDir/quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly_VP/*/N.ditissima/*/racon_10/medaka/medaka/pilon/pilon_10.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done

   for Assembly in $(ls assembly_VP/canu/N.ditissima/*/medaka/medaka/pilon/pilon_10.fasta); do
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
  done
```








## Rename contigs

```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    # If split or remove contigs is needed, provide FCSreport file by NCBI.
    touch tmp.txt
    for Assembly in $(ls assembly/flye/N.ditissima/Hg199/racon_10/assembly_racon_round_10.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
```


#

