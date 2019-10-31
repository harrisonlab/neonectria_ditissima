# R0905 isolate. Neonectria Reference_Genome_Assembly

======================================

This document details the commands used to assemble and annotate the R0905 genome.

Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/N.ditissima or /data/scratch/gomeza

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


# Genome assembly

## CANU

```bash
  for TrimReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $TrimReads | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    OutDir=assembly2/canu_pacbio/N.ditissima/"$Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction_pacbio.sh $TrimReads 45m $Strain $OutDir
  done
```
```bash
  for CorrectedReads in $(ls assembly2/canu_pacbio/N.d*/R0905/*.trimmedReads.fasta.gz); do
    Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly2/canu_pacbio/N.ditissima/"$Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_assembly_only_pacbio.sh $CorrectedReads 45m $Strain $OutDir
  done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/canu_pacbio/N.ditissima/R0905/R0905.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly2/canu_pacbio/N.ditissima/R0905/R0905.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

### Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
	for Assembly in $(ls assembly2/canu_pacbio/N.ditissima/R0905/R0905.contigs.fasta); do
  	Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    Iterations=10
    OutDir=assembly2/canu_pacbio/$Organism/$Strain/pilon
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done
```
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/canu_pacbio/N.ditissima/R0905/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls assembly2/canu_pacbio/N.ditissima/R0905/pilon/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## SMARTdenovo

```bash
  for RawReads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $RawReads | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f3 -d '/' | rev)
    Prefix="$Strain"_smartdenovo
    OutDir=assembly2/SMARTdenovo/N.ditissima/"$Strain"
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
    qsub $ProgDir/sub_SMARTdenovo.sh $RawReads $Prefix $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
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
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/R0905_smartdenovo.dmo.lay.utg); do
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

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/*10.fasta); do
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
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/*10.fasta); do
  OutDir=$(dirname $Assembly)
  $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/R0905_racon10_renamed.fasta); do
    Strain=R0905
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
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
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/pilon/*10.fasta); do
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
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.txt
  for Assembly in $(ls assembly2/SMARTdenovo/N.ditissima/R0905/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_pilon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```


## Miniasm

#### Assembly of uncorrected reads. Racon required after.

```bash
# If reads have same names or same part splited by space, fix them using rename.sh from bbtools.
# This can only be done in blacklace11, since it needs a specific library for blasr.
ssh blacklace11.blacklace
/home/gomeza/prog/bbtools/bbmap/rename.sh in=concatenated_pacbio.fastq.gz out=trimmed_renamed.fasta prefix=R0905

# Fast all-against-all overlap of raw reads
# Overlap for Pacbio reads (or use "-x ava-ont" for Minion read overlapping)
/home/gomeza/prog/minimap2/minimap2 -x ava-pb -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > R0905_fastq_allfiles.paf.gz

# Concatenate pieces of read sequences to generate the final sequences.
# Thus the per-base error rate is similar to the raw input reads. Make sure you correct your reads.
# Layout
/home/gomeza/prog/miniasm/miniasm -f trimmed_renamed.fasta R0905_fastq_allfiles.paf.gz > reads.gfa

# Convert gfa file to fasta file.
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > R0905_miniasm.fa
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls R0905_miniasm/R0905_miniasm.fa); do
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls R0905_miniasm/R0905_miniasm.fa); do
Strain=R0905
Organism=N.ditissima
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=R0905_miniasm/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Error correction using racon:

```bash
  for Assembly in $(ls R0905_miniasm/R0905_miniasm.fa); do
    Strain=R0905
    Organism=N.ditissima
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
  for Assembly in $(ls R0905_miniasm/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls R0905_miniasm/racon_10/*10.fasta); do
    Strain=R0905
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
  for Assembly in $(ls R0905_miniasm/racon_10/*10.fasta); do
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_racon10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
  done
  rm tmp.txt
```

## Pilon error correction

Assemblies were polished using Pilon
Note: qsub -R y 'Book blacklace11 avoiding more job in this node. Pilon requires a lot of memory'

```bash
for Assembly in $(ls R0905_miniasm/racon_10/R0905_racon10_renamed.fasta); do
    Strain=R0905
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
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
  for Assembly in $(ls R0905_miniasm/racon_10/pilon/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
```bash
  for Assembly in $(ls R0905_miniasm/racon_10/pilon/*10.fasta); do
    Strain=R0905
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

## Assembly using CANU. Default run.

```bash
	for Reads in $(ls raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq); do
  	GenomeSz="45m"
  	Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir=assembly_vAG/canu_1step/$Organism/"$Strain"
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/canu
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

### R0905 Assemblies were polished using Pilon

```bash
	for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
  	Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    Iterations=10
    OutDir=assembly2/canu_pacbio/$Organism/$Strain/Original_v3/polished
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done

  for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/R0905_canu.contigs.fasta); do
  	Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_v2)
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    Iterations=10
    OutDir=assembly2/canu_pacbio/$Organism/$Strain/Original_v3/polished_v2reads
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/canu_pacbio/N.di*/R0905/Original_v3/polished/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done

  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly2/canu_pacbio/N.di*/R0905/Original_v3/polished_v2reads/*10.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

```bash
  for Assembly in $(ls assembly2/canu_pacbio/N.di*/R0905/Original_v3/polished/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done

  for Assembly in $(ls assembly2/canu_pacbio/N.di*/R0905/Original_v3/polished_v2reads/*10.fasta); do
    Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
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
  for Assembly in $(ls R0905_good/pilon_5.fasta); do
    OutDir=R0905_good/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```
```bash
  for Assembly in $(ls R0905_good/filtered_contigs/R0905_renamed.fasta); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```
```bash
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/R0905_pilon10_renamed.fasta); do
    Strain=R0905
    Organism=N.ditissima
    echo "$Organism - $Strain"
    OutDir=$(dirname $Assembly)/repeat_masked/filtered_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```
The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash
  for File in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked.fa); do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    echo "$OutFile"
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
  done
```

```
Number of masked bases:
5386283
```

```bash
  for RepDir in $(ls -d assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs); do
    Strain=$(echo $RepDir | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f5 -d '/' | rev)
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

Identify Telomere repeats:
Telomeric repeats were identified in assemblies

```bash
  for Assembly in $(ls repeat_masked/Ref_Genomes/*/*/filtered_contigs/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/telomere/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
    $ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
  done
  cat $OutDir/telomere_hits.txt | sort -nr -k5 | less
```


## Merging Pacbio and hybrid assemblies

Note - the anchor length is the starting point for contigs to be merged - only contigs larger than these will be extended.

Note- The Python wrapper was written for MUMmer version 3.x command line. Hash out version 4 from profile while this is run.

### Canu

```bash
  for MinIONAssembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/pilon_10.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_pacbio/N.ditissima/R0905_all/filtered_contigs_min_500bp/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=900000
    OutDir=assembly_vAG/merged_assemblies/canu_spades/$Organism/"$Strain"_pacbio_900k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/canu_spades/$Organism/"$Strain"_hybrid_2100k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

### Miniasm

```bash
  for MinIONAssembly in $(ls assembly_vAG/miniasm/N.ditissima/R0905/racon_10/pilon/pilon_10.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_pacbio/$Organism/R0905_all/filtered_contigs_min_500bp/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=1300000
    OutDir=assembly_vAG/merged_assemblies/miniasm_spades/$Organism/"$Strain"_minion_1300k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/miniasm_spades/$Organism/"$Strain"_hybrid_1300k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

### SMARTdenovo

```bash
  for MinIONAssembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/R0905/racon_10/pilon/pilon_10.fasta); do
    Organism=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $MinIONAssembly | rev | cut -f4 -d '/' | rev)
    HybridAssembly=$(ls assembly_vAG/hybrid_spades/spades_pacbio/$Organism/R0905_all/filtered_contigs_min_500bp/contigs_min_500bp.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=900000
    OutDir=assembly_vAG/merged_assemblies/SMARTdenovo_spades/$Organism/"$Strain"_minion_900k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    OutDir=assembly_vAG/merged_assemblies/SMARTdenovo_spades/$Organism/"$Strain"_hybrid_900k
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $MinIONAssembly $OutDir $AnchorLength
  done
```

Quast and busco were run to assess the quality of hybrid assemblies:

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
  for Assembly in $(ls assembly_vAG/merged_assemblies/*/*/R0905_pac*/merged.fasta); do
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
