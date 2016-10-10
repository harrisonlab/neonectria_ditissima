Scripts used during analysis of neonectria_ditissima genomes.

neonectria_ditissima
====================

Commands used during analysis of the neonectria_ditissima genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

The following is a summary of the work presented in this Readme:

Draft Genome assembly Hg199
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Pacbio assembly
  * Merged hybrid assemblies
  * Rename contigs
  * Preliminary analysis

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

** Estimated Coverage is:


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

## Canu assembly contigs were renamed in accordance with ncbi recomendations (without repair contigs).


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

## Repair merged assemblies and rename.

Using IGV we can identify those missassembled regions and trim them. This will correct our assemblies.


```bash
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    touch tmp.csv
    # printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly" > tmp.csv
    printf \
    "Trim:
    Sequence name,\tlength,\tspan(s),\tapparent source
    contig_2\t3685016\t2630642..2631051\tlow coverage
    contig_7\t2122932\t14290..15776\tlow coverage
    contig_15\t1318801\t1303809..1306436\tlow coverage
    contig_17\t1243348\t1098076..1099539\tlow coverage
    contig_21\t602725\t120117..121484,162018..164688,289092..289649\tlow coverage
    " \
    > tmp.csv
    for Assembly in $(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/*_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
    done
    rm tmp.csv
```

# Repeatmasking assemblies

```bash
R0905_pacbio_merged=$(ls assembly/merged_canu_spades/*/R0905/edited_contigs/*_contigs_modified.fasta)
for Assembly in $(ls $R0905_pacbio_merged); do
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
OutDir=repeat_masking2/filtered_contigs
qsub $ProgDir/rep_modeling.sh $Assembly
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
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
## Repair merged assemblies

Using IGV we can identify those missassembled regions and trim them. This will correct our assemblies.


```bash
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    touch tmp.csv
    # printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly" > tmp.csv
    printf \
    "Trim:
    Sequence name,\tlength,\tspan(s),\tapparent source
    contig_2\t3685016\t2630642..2631051\tlow coverage
    contig_7\t2122932\t14290..15776\tlow coverage
    contig_15\t1318801\t1303809..1306436\tlow coverage
    contig_17\t1243348\t1098076..1099539\tlow coverage
    contig_21\t602725\t120117..121484,162018..164688,289092..289649\tlow coverage
    " \
    > tmp.csv
    for Assembly in $(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/*_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
    done
    rm tmp.csv
```

## Alignment of raw reads vs the Neonectria genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Neonectria genome


```bash
  Reference=$(ls repeat_masked/N.*/R0905_pacbio_canu/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/*); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_pacbio
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
  ```


# Gene Prediction


  Gene prediction followed three steps:
  	Pre-gene prediction
  		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
  	Gene model training
  		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
  	Gene prediction
  		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

  Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

  ```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
  	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  	for Genome in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  		echo $Genome;
  		qsub $ProgDir/sub_cegma.sh $Genome dna;
  	done
  ```
  ** Number of cegma genes present and complete: 94.35%
  ** Number of cegma genes present and partial: 95.97%


  Outputs were summarised using the commands:
  ```bash
  	for File in $(ls gene_pred/cegma/N.*/R0905_canu_assembly/*_dna_cegma.completeness_report); do
  		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
  		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
  		printf "$Species\t$Strain\n";
  		cat $File | head -n18 | tail -n+4;printf "\n";
  	done > gene_pred/cegma/cegma_results_dna_summary.txt

  	less gene_pred/cegma/cegma_results_dna_summary.txt
  ```


## Gene prediction

  Gene prediction was performed for Fusarium genomes. Two gene prediction
  approaches were used:

  Gene prediction using Braker1
  Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
  approach.


  ## Gene prediction 1 - Braker1 gene model training and prediction

  Gene prediction was performed using Braker1.

  First, RNAseq data was aligned to Fusarium genomes.
  * qc of RNA seq data is detailed below:

  ```bash
   for Folder in $(ls -d raw_rna/paired/N.ditissima/*); do
  	 FolderName=$(echo $Folder | rev | cut -f1 -d '/' | rev);
  	 echo $FolderName;
  	 ls $Folder/F;
  	 ls $Folder/R;
  	done
  ```
  This contained the following data:

  ```
  Hg199
  Hg199_S1_L001_R1_001_fastqc  Hg199_S1_L001_R1_001.fastq.gz
  Hg199_S1_L001_R2_001_fastqc  Hg199_S1_L001_R2_001.fastq.gz
  R0905
  R09-05_S2_L001_R1_001_fastqc  R09-05_S2_L001_R1_001.fastq.gz
  R09-05_S2_L001_R2_001_fastqc  R09-05_S2_L001_R2_001.fastq.gz
  ```

  Perform qc of RNAseq data

  ```bash
    for FilePath in $(ls -d raw_rna/paired/N.*/R0905); do
      echo $FilePath;
      FileF=$(ls $FilePath/F/*.gz);
      FileR=$(ls $FilePath/R/*.gz);
      IlluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa; ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc;
      qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA;
    done
  ```

  Data quality was visualised using fastqc:
  ```bash
  	for RawData in $(ls qc_rna/paired/N.ditissima/R0905/*/*.fq.gz); do
  		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/dna_qc
  		echo $RawData;
  		qsub $ProgDir/run_fastqc.sh $RawData
  	done
  ```

#### Aligning

  Insert sizes of the RNA seq library were unknown until a draft alignment could
  be made. To do this tophat and cufflinks were run, aligning the reads against a
  single genome. The fragment length and stdev were printed to stdout while
  cufflinks was running.

  ```bash
  	for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  		echo "$Organism - $Strain"
  		for RNADir in $(ls -d qc_rna/paired/N.ditissima/R0905); do
  			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
  			echo "$Timepoint"
  			FileF=$(ls $RNADir/F/*_trim.fq.gz)
  			FileR=$(ls $RNADir/R/*_trim.fq.gz)
  			OutDir=alignment/$Organism/$Strain/$Timepoint
  			ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
  		done
  	done
  ```

  76.4% overall read mapping rate.
  66.8% concordant pair alignment rate.

  Alignments were concatenated prior to running cufflinks:
  Cufflinks was run to produce the fragment length and stdev statistics:

  ```bash
  for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  AcceptedHits=alignment/$Organism/$Strain/R0905/accepted_hits.bam
  OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
  echo "$Organism - $Strain"
  mkdir -p $OutDir
  cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
  done
  ```

  Output from stdout included:
  ```
    Processed 19254 loci.                        [*************************] 100%
    Map Properties:
    Normalized Map Mass: 12059221.57
    Raw Map Mass: 12059221.57
    Fragment Length Distribution: Empirical (learned)
                  Estimated Mean: 219.68
  	            Estimated Std Dev: 39.56
  [15:04:46] Assembling transcripts and estimating abundances.
    Processed 19333 loci.                        [*************************] 100%
  ```

  The Estimated Mean: 219.68 allowed calculation of of the mean insert gap to be
  -140bp 182-(180*2) where 180? was the mean read length. This was provided to tophat
  on a second run (as the -r option) along with the fragment length stdev to
  increase the accuracy of mapping.


  Then Rnaseq data was aligned to each genome assembly:

  ```bash
  for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/paired/N.ditissima/R0905 | grep -v -e '_rep'); do
  Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
  echo "$Timepoint"
  FileF=$(ls $RNADir/F/*_trim.fq.gz)
  FileR=$(ls $RNADir/R/*_trim.fq.gz)
  OutDir=alignment/$Organism/$Strain/$Timepoint
  InsertGap='-140'
  InsertStdDev='40'
  Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 10
  printf "."
  Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
  done
  printf "\n"
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
  qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
  done
  done

    cd alignment/N.ditissima/R0905_pacbio_canu/
    mkdir R0905_accurate
    mv -r R0905/* R0905_accurate/
  ```
  76.4% overall read mapping rate.
  66.7% concordant pair alignment rate.

#### Braker prediction

  Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

  ```bash
  ls ~/.gm_key
  cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
  ```
  ```bash
    for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
      Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
      while [ $Jobs -gt 1 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
      done
      printf "\n"
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      OutDir=gene_pred/braker/$Organism/"$Strain"_braker_new
      AcceptedHits=alignment/$Organism/R0905_merged_assembly/R0905/accepted_hits.bam
      GeneModelName="$Organism"_"$Strain"_braker_new
      rm -r /home/gomeza/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
      qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
    done
  ```

  Fasta and gff files were extracted from Braker1 output.

  ```bash
    for File in $(ls gene_pred/braker/N.*/R0905_braker_new/*/augustus.gff); do
      getAnnoFasta.pl $File
      OutDir=$(dirname $File)
      echo "##gff-version 3" > $OutDir/augustus_extracted.gff
      cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
    done
  ```

  The relationship between gene models and aligned reads was investigated. To do
  this aligned reads needed to be sorted and indexed:

  Note - IGV was used to view aligned reads against the Fus2 genome on my local
  machine.


## Supplimenting Braker gene models with CodingQuary genes

  Additional genes were added to Braker gene predictions, using CodingQuary in
  pathogen mode to predict additional regions.

  Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

  Note - cufflinks doesn't always predict direction of a transcript and
  therefore features can not be restricted by strand when they are intersected.

  ```bash
  		for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  			Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  			Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  			echo "$Organism - $Strain"
  			OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
  			mkdir -p $OutDir
  			AcceptedHits=alignment/$Organism/R0905_merged_assembly/*/accepted_hits.bam
  			ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
  			qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  		done
  ```

  Secondly, genes were predicted using CodingQuary:

  ```bash
  		for Assembly in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  		echo "$Organism - $Strain"
  		OutDir=gene_pred/codingquary/$Organism/$Strain
  		CufflinksGTF=gene_pred/cufflinks/$Organism/R0905_merged_assembly/concatenated_prelim/cufflinks/transcripts.gtf
  		ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary
  		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  		done
  ```

  Then, additional transcripts were added to Braker gene models, when CodingQuary
  genes were predicted in regions of the genome, not containing Braker gene
  models:

  ```bash
  	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
  for BrakerGff in $(ls gene_pred/braker/N.*/*_braker_new/*/augustus.gff3); do
  Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g')
  Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Assembly=$(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
  PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
  AddDir=gene_pred/codingquary/$Organism/$Strain/additional
  FinalDir=gene_pred/codingquary/$Organism/$Strain/final
  AddGenesList=$AddDir/additional_genes.txt
  AddGenesGff=$AddDir/additional_genes.gff
  FinalGff=$AddDir/combined_genes.gff
  mkdir -p $AddDir
  mkdir -p $FinalDir

  bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
  bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
  $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary

  $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
  $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
  cp $BrakerGff $FinalDir/final_genes_Braker.gff3
  $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
  cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
  cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
  cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
  cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

  GffBraker=$FinalDir/final_genes_CodingQuary.gff3
  GffQuary=$FinalDir/final_genes_Braker.gff3
  GffAppended=$FinalDir/final_genes_appended.gff3
  cat $GffBraker $GffQuary > $GffAppended

  # cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
  done
  ```

  The final number of genes per isolate was observed using:
  ```bash
  for DirPath in $(ls -d gene_pred/codingquary/N.*/*/final); do
  echo $DirPath;
  cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
  cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
  cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
  echo "";
  done
  ```

  gene_pred/codingquary/N.ditissima/R0905_merged_assembly/final
  12686
  804
  13490


## ORF finder

  The genome was searched in six reading frames for any start codon and following
  translated identification of a start codon translating sequence until a stop
  codon was found. This is based upon the atg.pl script used in paper describing
  the P. infestans genome. Additional functionality was added to this script by
  also printing ORFs in .gff format.


  ```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  	for Genome in $(ls repeat_masked/*/*/edited_*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  	qsub $ProgDir/run_ORF_finder.sh $Genome
  	done
  ```

  The Gff files from the the ORF finder are not in true Gff3 format. These were
  corrected using the following commands:

  ```bash
  	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  	for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
  		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_corrected.gff3/g')
  		echo ""
  		echo "Correcting the following file:"
  		echo $ORF_Gff
  		echo "Redirecting to:"
  		echo $ORF_Gff_mod
  		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
  	done
  ```

#Functional annotation

## A) Interproscan

  Interproscan was used to give gene models functional annotations.
  Annotation was run using the commands below:

  Note: This is a long-running script. As such, these commands were run using
  'screen' to allow jobs to be submitted and monitored in the background.
  This allows the session to be disconnected and reconnected over time.

  Screen ouput detailing the progress of submission of interporscan jobs
  was redirected to a temporary output file named interproscan_submission.log .

  ```bash
  	screen -a
  	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  	for Genes in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
  	echo $Genes
  	$ProgDir/sub_interproscan.sh $Genes
  	done 2>&1 | tee -a interproscan_submisison.log
  ```

  Following interproscan annotation split files were combined using the following
  commands:

  ```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  	for Proteins in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
  		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
  		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
  		echo "$Organism - $Strain"
  		echo $Strain
  		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
  		$ProgDir/append_interpro.sh $Proteins $InterProRaw
  	done
  ```


## B) SwissProt

  ```bash
  	for Proteome in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
  		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  		OutDir=gene_pred/swissprot/$Organism/$Strain
  		SwissDbDir=../../uniprot/swissprot
  		SwissDbName=uniprot_sprot
  		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
  		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  	done
  ```

  ```bash
  	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl); do
  		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
  		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
  		echo "$Organism - $Strain"
  		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
  		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
  		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
  	done
  ```


#Genomic analysis

  ## Comparison to FoL 4287

  BLast searches were used to identify which genes had homologs on which
  chromosomes of the Fusarium lycopersici genome.

  ```bash
  	FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
  	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
  	# for Proteome in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta); do
  	# for Proteome in $(ls gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa); do
  		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  		echo "$Organism - $Strain"
  		OutDir=analysis/blast_homology/$Organism/$Strain
  		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  		qsub $ProgDir/run_blast2csv.sh $Proteome protein $FoLGenomeFa $OutDir
  	done
  ```

## Effector genes

  Putative pathogenicity and effector related genes were identified within Braker
  gene models using a number of approaches:

   * A) From Augustus gene models - Identifying secreted proteins
   * B) From Augustus gene models - Effector identification using EffectorP
   * C) CAZY proteins


### A) From Augustus gene models - Identifying secreted proteins

   Required programs:
    * SignalP-4.1
    * TMHMM

   Proteins that were predicted to contain signal peptides were identified using
   the following commands:

   ```bash
   	SplitfileDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   	CurPath=$PWD
   	for Proteome in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
   		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   		SplitDir=gene_pred/final_genes_split/$Organism/$Strain
   		mkdir -p $SplitDir
   		BaseName="$Organism""_$Strain"_final_preds
   		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
   		for File in $(ls $SplitDir/*_final_preds_*); do
   			Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
   			while [ $Jobs -gt 1 ]; do
   				sleep 10
   				printf "."
   				Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
   			done
   			printf "\n"
   			echo $File
   			qsub $ProgDir/pred_sigP.sh $File signalp-4.1
   		done
   	done
   ```

   The batch files of predicted secreted proteins needed to be combined into a
   single file for each strain. This was done with the following commands:
   ```bash
   	for SplitDir in $(ls -d gene_pred/final_genes_split/N.*/R0905_merged_assembly); do
   		Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
   		Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
   		InStringAA=''
   		InStringNeg=''
   		InStringTab=''
   		InStringTxt=''
   		SigpDir=final_genes_signalp-4.1
   		for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
   			InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
   			InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
   			InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
   			InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
   		done
   		cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
   		cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
   		tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
   		cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
   	done
   ```

   Some proteins that are incorporated into the cell membrane require secretion.
   Therefore proteins with a transmembrane domain are not likely to represent
   cytoplasmic or apoplastic effectors.

   Proteins containing a transmembrane domain were identified:

   ```bash
   	for Proteome in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
   		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
   		qsub $ProgDir/submit_TMHMM.sh $Proteome
   	done
   ```

### B) From Augustus gene models - Effector identification using EffectorP

  Required programs:
   * EffectorP.py

  ```bash
  	for Proteome in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
  		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  		BaseName="$Organism"_"$Strain"_EffectorP
  		OutDir=analysis/effectorP/$Organism/$Strain
  		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
  		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  	done
  ```

### C) CAZY proteins

  Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
  at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

  ```bash
  	for Proteome in $(ls gene_pred/codingquary/N.*/R0905_merged_assembly/*/final_genes_combined.pep.fasta); do
  		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  		OutDir=gene_pred/CAZY/$Organism/$Strain
  		mkdir -p $OutDir
  		Prefix="$Strain"_CAZY
  		CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
  		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
  		qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  	done
  ```

  The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

  Those proteins with a signal peptide were extracted from the list and gff files
  representing these proteins made.

  ```bash
    for File in $(ls gene_pred/CAZY/N.*/R0905_pacbio_canu/*CAZY.out.dm); do
      Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
      Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
      OutDir=$(dirname $File)
      echo "$Organism - $Strain"
      ProgDir=/home/groups/harrisonlab/dbCAN
      $ProgDir/hmmscan-parser.sh $OutDir/R0905_merged_assembly_CAZY.out.dm > $OutDir/R0905_merged_assembly_new_CAZY.out.dm.ps
      SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
      SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
      cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
      Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
      CazyGff=$OutDir/R0905_merged_assembly_new_CAZY.gff
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $Gff CAZyme ID > $CazyGff
    done
  ```

  Note - the CAZY genes identified may need further filtering based on e value and
  cuttoff length - see below:

  Cols in yourfile.out.dm.ps:
  1. Family HMM
  2. HMM length
  3. Query ID
  4. Query length
  5. E-value (how similar to the family HMM)
  6. HMM start
  7. HMM end
  8. Query start
  9. Query end
  10. Coverage

  * For fungi, use E-value < 1e-17 and coverage > 0.45

  * The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)
