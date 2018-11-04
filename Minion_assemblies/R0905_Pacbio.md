# R0905 isolate. Neonectria_Reference_Genome_Assembly

======================================

This document details the commands used to assemble and annotate the R0905 genome.

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


## Genome assembly

R09/05 - I will test different assembly methods. Canu in 2 steps and in 1 step.

#### R09/05 - Assembly using CANU

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

#### R0905 - Assembly using CANU 1 step

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

#### R09/05 - Assembbly using SMARTdenovo

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

#### Racon correction

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

#### R0905 Assembly correction using nanopolish

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

### R0905 Assemblies were polished using Pilon

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

### R0905 Hybrid Spades Assembly

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
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/spades_pacbio
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
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

Quast

```bash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs_min_500bp/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
```
```bash
    for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs_min_500bp/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/spades_pacbio
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

### Merging pacbio and hybrid assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/polished/pilon_5.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp/contigs_min_500bp.fasta)
    AnchorLength=5000
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/"$Strain"_pacbio_5k
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```
```bash
  for PacBioAssembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/racon_10/polished/pilon_5.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f6 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f5 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp/contigs_min_500bp.fasta)
    AnchorLength=5000
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/"$Strain"_pacbio_afterracon_5k
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Merged=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/SMARTdenovo/$Organism/$Strain/$Merged/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/merged_canu_spades/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Merged=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/merged_canu_spades/$Merged
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

This merged assembly was polished using Pilon.

```bash
  for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/merged.fasta); do
    IlluminaDir=$(ls -d qc_dna/paired/*/R0905)
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    OutDir=$(dirname $Assembly)/pilon
    Iterations=5
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/pilon/*5.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/pilon/*5.fasta); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```
Same statistics before pilon

The "best" genome was obtained after polishing the canu genome with pilon.


### Improving R0905 assembly.

### R0905 Hybrid Spades Assembly

```bash
  	for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $StrainPath
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/R0905_all)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    OutDir=assembly/spades_pacbio/$Organism/"$Strain"_all
    qsub -R y $ProgDir/sub_spades_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $OutDir 15
  	done
```
```bash
  	for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $StrainPath
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/R0905_v2)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    OutDir=assembly/spades_pacbio/$Organism/"$Strain"_v2
    qsub -R y $ProgDir/sub_spades_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $OutDir 15
  	done
```
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/spades_pacbio/*/R0905_*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/spades_pacbio/*/R0905_*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/spades_pacbio/$Organism/$Strain/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

Contigs shorter than 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_pacbio/*/R0905_*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs_min_500bp
    FilterDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs_min_500bp/contigs_min_500bp.fasta
  done
```

Quast

```bash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    for Assembly in $(ls assembly/spades_pacbio/*/R0905_*/filtered_contigs_min_500bp/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
```
```bash
    for Assembly in $(ls assembly/spades_pacbio/*/R0905_*/filtered_contigs_min_500bp/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

### Merging pacbio and hybrid assemblies

I used the pilon_5 genome because busco prediction was slightly better than racon genome.

```bash
for PacBioAssembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/polished/pilon_5.fasta); do
  HybridAssembly=assembly/spades_pacbio/*/R0905_all/filtered_contigs_min_500bp/contigs_min_500bp.fasta
  Organism=$(echo $HybridAssembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $HybridAssembly | rev | cut -f3 -d '/' | rev)
  AnchorLength=500000
  OutDir=assembly/merged_canu_spades/$Organism/$Strain/"$Strain"_pacbio_500k
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
  qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
done
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Merged=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/SMARTdenovo/$Organism/$Strain/$Merged/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/merged_canu_spades/*/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Merged=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/Ref_Genomes/$Strain/merged_canu_spades/$Merged
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

This merged assembly was polished using Pilon.

```bash
  for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/merged.fasta); do
    IlluminaDir=$(ls -d qc_dna/paired/*/R0905)
    TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
    OutDir=$(dirname $Assembly)/pilon
    Iterations=5
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
  done
```
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/pilon/*5.fasta); do
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/merged_canu_spades/N.dit*/R0905/R0905_pacbio_5k/pilon/*5.fasta); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=$(dirname $Assembly)/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```


After pilon, only two additional gene was predicted in the R0905 genome.
Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/CSAR/Scaffold_v2/N.ditissima/*/polished/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/CSAR/Scaffold_v2/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

## Repeat masking of pilon_5. The best R0905 genome.


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
for Assembly in $(ls R0905_good/filtered_contigs/R0905_renamed.fasta); do
Strain=R0905
Organism=N.ditissima
echo "$Organism - $Strain"
OutDir=R0905_good/repeat_masked/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```
The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.









```bash

for File in $(ls repeat_masked/Ref_Genomes/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
```
```
repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4060878
repeat_masked/Ref_Genomes/N.ditissima/R0905/filtered_contigs/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5516834
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
