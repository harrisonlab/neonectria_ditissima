Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

The following is a summary of the work presented in this Readme:

Draft Genome assembly Hg199
  * Data qc
  * Genome assembly
  * Repeatmasking

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
  Read_F=raw_dna/paired/N.ditissima/*/F/*.fastq.gz
  Read_R=raw_dna/paired/N.ditissima/*/R/*.fastq.gz
  IluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
  qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/N.ditissima/*/*/*.fq.gz); do
  echo $RawData;
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc;
  qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
	Trim_F=qc_dna/paired/N.ditissima/*/F/*.fq.gz
	Trim_R=qc_dna/paired/N.ditissima/*/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

mode kmer abundance prior to error correction was reported using the following commands:
```bash
    for File in $(ls qc_dna/kmc/N.ditissima/*/*_true_kmer_summary.txt); do
        basename $File;
        cat $File | grep -e 'abundance' -e 'size'
    done
```

** Estimated Genome Size is:

** Esimated Coverage is:


#Assembly
Assembly was performed using: Velvet / Abyss / Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
	F_Read=qc_dna/paired/N.ditissima/*/F/*.fq.gz
	R_Read=qc_dna/paired/N.ditissima/*/R/*.fq.gz
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Outdir=assembly/spades/N.ditissima/Hg199/
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $Outdir correct
```

cat contigs.fasta | grep '>' | wc -l

cat contigs_min_500bp.fasta | grep '>' | wc -l


```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/N.ditissima/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:
  * N80:
  * N20:
  * Longest contig:
  **


# Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss=/assembly/spades/N.ditissima/*/filtered_contigs/contigs_min_500bp.fasta
	qsub $ProgDir/rep_modeling.sh $BestAss
	qsub $ProgDir/transposonPSI.sh $BestAss
 ```

** % bases masked by repeatmasker: % (bases masked:  bp)

** % bases masked by transposon psi: % (bases masked:  bp)

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

    repeat_masked/N.ditissima/*/filtered_contigs_repmask/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
    Number of masked bases:  
```

## Pacbio

### Canu assembly

```bash
  	Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq)
  	GenomeSz="46m"
  	Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir="assembly/canu/$Organism/$Strain"
  	ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  	qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

Assemblies were polished using Pilon

```bash
  	for Assembly in $(ls assembly/canu/*/R0905/*_canu.contigs.fasta); do
  	Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/NG-R0905)
    TrimF1_Read=$(ls $IlluminaDir/F/NG-R0905_qc_F.fastq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/NG-R0905_qc_R.fastq.gz);
    OutDir=assembly/canu/$Organism/$Strain/polished
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  	done
```
cat pilon.fasta | grep 'tig' | wc -l
125

### Spades Assembly


For N. ditissima

```bash
  	for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $StrainPath
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/NG-R0905)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/NG-R0905_qc_F.fastq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/NG-R0905_qc_R.fastq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    OutDir=assembly/spades_pacbio/$Organism/$Strain
    qsub $ProgDir/sub_spades_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $OutDir 15
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
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
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
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

#Contigs were renamed in accordance with ncbi recomendations.


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

#Canu assembly contigs were renamed in accordance with ncbi recomendations.


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

** % bases masked by repeatmasker: % (bases masked:  bp)

** % bases masked by transposon psi: % (bases masked:  bp)


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
Number of masked bases:

# Preliminary analysis

## Checking PacBio coverage against Neonectria contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.


Merged canu spades assembly

```bash
    Assembly=assembly/merged_canu_spades/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/merged_canu_spades/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```

Canu assembly

```bash
    Assembly=assembly/canu/N.ditissima/R0905/polished/pilon.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/canu/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```





```bash
  Reference=$(ls repeat_masked/N.*/R0905_pacbio_canu/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/NG-R0905); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_NG-R0905_pacbio2
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
  ```
