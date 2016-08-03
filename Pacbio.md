
## Data extraction


for N. ditissima data:
```bash
  	cd /home/groups/harrisonlab/project_files/idris
  	RawDatDir=/home/harrir/projects/pacbio_test/n_dit/
  	mkdir -p raw_dna/pacbio/N.ditissima/R0905
  	cp -r $RawDatDir/D08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/E08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/F08_1 raw_dna/pacbio/N.ditissima/R0905/.
    OutDir=raw_dna/pacbio/N.ditissima/R0905/extracted
  	mkdir -p $OutDir
  	cat raw_dna/pacbio/N.ditissima/R0905/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
```


## Assembly


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

# Preliminary analysis

## Checking PacBio coverage against Fus2 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.


Merged canu spades assembly

```bash
    Assembly=assembly/merged_canu_spades/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```

Canu assembly

```bash
    Assembly=assembly/canu/N.ditissima/R0905/polished/pilon.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
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

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2' -e 'fo47'); do
	#for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w 'fo47'); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

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


#Gene prediction

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
	for Assembly in $(ls repeat_masked/*/R0905_pacbio_canu/*/*_contigs_unmasked.fa); do
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
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/R0905_pacbio_canu/*/*_contigs_softmasked.fa); do
	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	# AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
  AcceptedHits=alignment/$Organism/$Strain/R0905/accepted_hits.bam
  OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
	echo "$Organism - $Strain"
	mkdir -p $OutDir
	# samtools merge -f $AcceptedHits \
	# alignment/$Organism/$Strain/R0905/accepted_hits.bam \
  # alignment/$Organism/$Strain/R0905/accepted_hits.bam
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
	for Assembly in $(ls repeat_masked/*/R0905_pacbio_canu/*/*_contigs_unmasked.fa); do
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
```


#### Braker prediction

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
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
	mkdir -p alignment/$Organism/$Strain/concatenated
	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
	alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
	OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
	AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
	GeneModelName="$Organism"_"$Strain"_braker_new
	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
	qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
	done
```

Fasta and gff files were extracted from Braker1 output.

```bash
	for File in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff); do
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

```bash
	InBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated.bam
	ViewBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_view.bam
	SortBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_sorted
	samtools view -b $InBam > $ViewBam
	samtools sort $ViewBam $SortBam
	samtools index $SortBam.bam
```


Cufflinks was run to compare the predicted genes to assembled transcripts:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		done
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```

<!--
The number of Fo47 genes was determined for comparison to number predicted by Braker (16269):
```bash
	fo47_transcripts=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf
	cat $fo47_transcripts | grep 'gene_id' | cut -f2 -d '"' | sort | uniq | wc -l
	# 18191
```
 -->



## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		mkdir -p $OutDir
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
	for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_pacbio/*/augustus.gff3 | grep -e 'Fus2'); do
		Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
		Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
		Assembly=$(ls repeat_masked/$Organism/$Strain/*/"$Strain"_contigs_softmasked.fa)
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
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
		$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

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
for DirPath in $(ls -d gene_pred/codingquary/F.*/*/final | grep -w -e'Fus2'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
<!--
## Suplimenting gene models with known genes

Additional gene models were consrtucted in braker / transfered from other isolates
and were exported to the following locations:

```bash
	Fus2Six9=assembly/spades/F.oxysporum_fsp_cepae/Fus2_edited_v2/filtered_contigs/Fus2_edited_v2_contig_1090_six9.gff
```

These gene models were then edited with nano to give names and IDs to these genes. -->

## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
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
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w 'Fus2'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```


## B) SwissProt
<!--
```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
    mkdir -p $OutDir
    blastp \
    -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
    -query $ProjDir/$Proteome \
    -out $OutDir/swissprot_v2015_10_hits.tbl \
    -evalue 1e-100 \
    -outfmt 6 \
    -num_threads 16 \
    -num_alignments 10
  done
``` -->


```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred/swissprot/$Organism/$Strain
		SwissDbDir=../../uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl | grep -w -e 'Fus2'); do
	# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
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

Convert top blast hits into gff annotations

```bash
	for BlastHitsCsv in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv | grep -w -e 'Fus2'); do
		Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
		Column2="$Strain"_gene_homolog
		NumHits=1
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
	done
```

#### Intersecting blast hits with genes from FoL

```bash
	for HitsGff in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.gff | grep -w -e 'Fus2'); do
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $HitsGff| rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsDir=$(dirname $HitsGff)
		FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
		FoLIntersect=$HitsDir/4287_chromosomal_final_genes_combined_intersect.bed
		bedtools intersect -wo -a $HitsGff -b $FoLGenes > $FoLIntersect
	done
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP
 * D) From ORF fragments - Signal peptide & RxLR motif
 * E) From ORF fragments - Hmm evidence of WY domains
 * F) From ORF fragments - Hmm evidence of RxLR effectors


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
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
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/Fus2); do
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
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
		qsub $ProgDir/submit_TMHMM.sh $Proteome
	done
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```

### C) Identification of MIMP-flanking genes

```bash
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
		Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
		Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
		BrakerGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_CodingQuary.gff3)
		QuaryGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_Braker.gff3)
		OutDir=analysis/mimps/$Organism/$Strain
		mkdir -p "$OutDir"
		echo "$Organism - $Strain"
		ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
		$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
		$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
		echo "The number of mimps identified:"
		cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
		bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		bedtools intersect -u -a $QuaryGff -b $OutDir/"$Strain"_mimps_exp.gff >> $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		echo "The following transcripts intersect mimps:"
		MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
		cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' > $MimpGenesTxt
		cat $MimpGenesTxt | wc -l
		echo ""
	done
```

```
	F.oxysporum_fsp_cepae - Fus2
	The number of mimps identified:
	151
	The following transcripts intersect mimps:
	155
```

```bash
	cat gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp.aa | grep '>' | sed 's/>//g' | cut -f1 -d '-' | sed -r 's/\s//g' > gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp_headers.txt
	cat analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP.txt | grep 'Effector' | cut -f1 | sed 's/|//g' > analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP_headers.txt
	cat gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos.txt | cut -f1 > gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos_headers.txt
```



# 4. Genomic analysis
<!--
## 4.1 Identifcation of protospacers

To facilitate CriprCas editing of the Fusarium oxysporum genome target sites
known as protospacers were identified.

This was done using a published program OPTIMuS as well as a parser script. The
commands to do this were as follows:

```bash
	for GeneSeq in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.codingseq | grep -v 'old'); do
	  Organism=$(echo $GeneSeq | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $GeneSeq | rev | cut -f2 -d '/' | rev)
		echo $Organism
		echo $Strain
	  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
	  OutDir=analysis/protospacers/$Organism/$Strain
	  mkdir -p $OutDir
	  $ProgDir/OPTIMuS_EMR.pl $GeneSeq "threshold" 1 > $OutDir/"$Strain"_protospacer_sites.txt
	  $ProgDir/Optimus2csv.py --inp $OutDir/"$Strain"_protospacer_sites.txt  --out $OutDir/"$Strain"_protospacer_by_gene.csv
	done
``` -->


## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |



<!--
### 4.2.a) Orthology between pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogenic_orthology.md


### 4.2.b) Orthology between non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_non-pathogenic_orthology.md
 -->

### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md


<!-- ### 4.2.d) Orthology between all isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_isolates.md -->


## 5.  Searches

## 5.1.A) Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemblies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2'); do
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv | grep -w 'Fus2'); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=SIX_homolog
		NumHits=1
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	done
```

	The blast hits were summarised in a single table for all the genomes. The top
	identity of the top blast hit in relation to the enquire query sequence was
	presented for each blast hit.

```bash
	OutFile=analysis/blast_homology/Fo_path_genes_CRX_summary.tab
	echo "Organism" > tmp2.tab
	cat analysis/blast_homology/F.proliferatum/A8/A8_Fo_path_genes_CRX.fa_homologs.csv | cut -f1 >> tmp2.tab
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		echo "$Organism" > tmp.tab
		echo "$Strain" >> tmp.tab
		cat $BlastHits | cut -f10  | tail -n +2 >> tmp.tab
		paste tmp2.tab tmp.tab > $OutFile
		cp $OutFile tmp2.tab
	done
	rm tmp.tab
	rm tmp2.tab
```

```bash
	for HitsGff in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.gff | grep -v 'trinity' | grep -w 'Fus2'); do
		Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		# GffBraker=gene_pred/codingquary/$Organism/$Strain/final/final_genes_Braker.gff3
		# GffQuary=gene_pred/codingquary/$Organism/$Strain/final/final_genes_CodingQuary.gff3
		GffAppended=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		OutDir=$(dirname $HitsGff)
		SixIntersect=$OutDir/"$Strain"_Fo_path_genes_CRX.fa_hit_genes.bed
		# bedtools intersect -wo -a $HitsGff -b $GffBraker > $SixIntersect
		# bedtools intersect -wo -a $HitsGff -b $GffQuary >> $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended > $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended | cut -f9,18 | grep -v 'Parent'
		echo ""
	done > analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab
```


## 5.1.B) Identifying FTF genes

Previously published FTF genes from Sanchez et al 2016 were blasted against
Fusarium genomes.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'cepae' | grep -v -e '_edited_v2' -e 'HB17'); do
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
		OutDir=analysis/FTF/$Organism/$Strain
		mkdir -p $OutDir
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
	done
```

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
	for BlastHits in $(ls analysis/FTF/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_hits.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/FTF/$Organism/$Strain
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=FTF_homolog
		NumHits=1
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

		GffAppended=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
		bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_FTF_hits_intersected.bed
	done
```



## 5.2 Identifying PHIbase homologs

The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.
<!--
```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	printf "header\n" > $PhibaseDir/PHI_headers.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/\r//g' >> $PhibaseDir/PHI_headers.csv
	printf "effect\n" > .$PhibaseDir/PHI_virulence.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '|' | rev  >> $PhibaseDir/PHI_virulence.csv
```


```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	PhibaseHeaders=$PhibaseDir/PHI_headers.csv
	PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
	for BlastCSV in $(ls analysis/blast_homology/F*/*/*_PHI_36_accessions.fa_homologs.csv); do
		Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
		echo "$Strain"
		OutDir=$(dirname $BlastCSV)
		paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
		cat $OutDir/"$Strain"_PHIbase_virulence.csv | grep 'NODE_' | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
	done
```
results were:

```
	125
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     66 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    131 Reduced virulence
	     87 Unaffected pathogenicity
	55
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     58 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    121 Reduced virulence
	     83 Unaffected pathogenicity
	A23
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     79 Unaffected pathogenicity
	A28
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     58 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    125 Reduced virulence
	     79 Unaffected pathogenicity
	D2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     62 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     81 Unaffected pathogenicity
	Fus2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     59 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    120 Reduced virulence
	     79 Unaffected pathogenicity
	HB17
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     84 Unaffected pathogenicity
	PG
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     61 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    126 Reduced virulence
	     78 Unaffected pathogenicity
	N139
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     72 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    149 Reduced virulence
	     93 Unaffected pathogenicity
	PG18
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     72 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    153 Reduced virulence
	     97 Unaffected pathogenicity
	PG3
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     65 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    139 Reduced virulence
	     93 Unaffected pathogenicity
	A8
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      4 Increased virulence (Hypervirulence)
	     15 Lethal
	     71 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    152 Reduced virulence
	     82 Unaffected pathogenicity

```
-->

<!-- The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Proteins in $(ls gene_pred/braker/F.oxysporum_fsp_cepae/*/*/*.aa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteins protein ../../phibase/v3.8/PHI_accessions.fa
	done
``` -->

### Expressed Genes

As a preliminary analysis of the RNAseq data, highly expressed genes at 72hrs
post infection were identified in Fus2.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	# samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam
	for AcceptedHits in $(ls alignment/*/*/concatenated/concatenated.bam | grep -v -e 'Fus2_edited_v2'); do
		Strain=$(echo $AcceptedHits | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $AcceptedHits | rev | cut -f4 -d '/' | rev)
		OutDir=$(dirname $AcceptedHits)
		echo "$Organism - $Strain"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

```bash
	for Transcripts in $(ls alignment/F.*/*/concatenated/transcripts.gtf); do
		Strain=$(echo $Transcripts | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Transcripts | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		GeneGff=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		ExpressedGenes=alignment/$Organism/$Strain/concatenated/expressed_genes.bed
		bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
	done
```


<!--
```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=Fus2
	samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam

	cufflinks -o timecourse/2016_genes/Fus2/72hrs/cufflinks -p 16 --max-intron-length 4000 alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam

	Transcripts=timecourse/2016_genes/Fus2/72hrs/cufflinks/transcripts.gtf
	GeneGff=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_appended.gff3
	ExpressedGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes.gff
	bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
```
-->

<!--
It was noted that not all the top expressed transcripts had gene models. The
transcripts without gene models were identified:

```bash
	ExpressedTranscriptsTxt=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | cut -f2 -d'"' | uniq | head -n20 > $ExpressedTranscriptsTxt
	ExpressedTranscriptsInter=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts_intersect.gff
	cat $ExpressedGenes | grep -v -P "\ttranscript_id" | sort -r -n -t'"' -k6 | grep -w -f $ExpressedTranscriptsTxt > $ExpressedTranscriptsInter
	ExpressedTranscriptsGff=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.gff
	cat $ExpressedTranscriptsInter | cut -f1-9 | uniq> $ExpressedTranscriptsGff
	cat $ExpressedTranscriptsInter | grep -v -w 'exon' | cut -f2,18 -d'"' --output-delimiter " - " | uniq | less
```

The top 20 expressed transcripts are shown below with genes they correspond to.
9 of the 20 top expressed transcripts have not been predicted as genes. This was
determined, in large, to be a result of transposon being expressed. It was also
noted that cufflinks doesn't always predict direction of a transcript and
therefore when features are intersected, they can not be restricted by strand.

```
	CUFF.4396 - 33 bp - A-rich
	CUFF.2444 - 69 bp - T-rich - repmasked
	CUFF.540 - 60 bp - C-Rich - repmasked
	CUFF.11337 - 92 bp - T-rich - repmasked
	CUFF.14001 - g11792 - Six5
	CUFF.5007 - repmasked
	CUFF.4291 - repmasked
	CUFF.12587 - g11716 - Six3
	CUFF.4090 - repmasked
	CUFF.12214 - g4762
	CUFF.8542 - 103 bp - A-rich - repmasked
	CUFF.13998 - g12630
	CUFF.13635 - g10978
	CUFF.4290 - no gene - Six9
	CUFF.2922 - g12344
	CUFF.9876 - g10974
	CUFF.13435 - g10975
	CUFF.9559 - g10973
	CUFF.12530 - g4720
	CUFF.6539 - g10346
```
-->


<!--
```bash
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | grep -w 'gene' | head -n 20 | cut -f18 > $TopGenes
	InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
	TopFunctions=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.tsv
	cat $InterPro | grep -w -f $TopGenes > $TopFunctions

	IsolateAbrv=FoC_path_vs_non_path
	WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
	Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
	PathDir=$WorkDir/path_orthogroups
	echo "The number of pathogen common genes in the top 20 expressed Fus2 genes at 72hr"
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $TopGenes | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt
	# secreted
	PathOrthogroupsFus2Secreted=$PathDir/Fus2_path_orthogroup_secreted.txt
	cat $PathOrthogroupsFus2Secreted | grep -w -f $TopGenes
	# near mimps
	MimpGenesTxt=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_genes_in_2kb_mimp.txt
  cat $MimpGenesTxt | grep -w -f $TopGenes
	# Intersecting SIX genes

``` -->


## 6. Summarising the Fusarium Proteome

```bash# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/fusarium
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2' -e 'fo47'); do
	#for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w 'fo47'); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/F*/Fus2/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```


#Gene prediction

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
 for Folder in $(ls -d raw_rna/paired/F.oxysporum_fsp_cepae/*); do
	 FolderName=$(echo $Folder | rev | cut -f1 -d '/' | rev);
	 echo $FolderName;
	 ls $Folder/F;
	 ls $Folder/R;
	done
```
This contained the following data:
```
	55_72hrs_rep1
	sample013_1.combined.fastq.gz
	sample013_2.combined.fastq.gz
	55_72hrs_rep2
	sample014_1.combined.fastq.gz
	sample014_2.combined.fastq.gz
	55_72hrs_rep3
	sample015_1.combined.fastq.gz
	sample015_2.combined.fastq.gz
	control_72hrs_rep1
	sample002_1.combined.fastq.gz
	sample002_2.combined.fastq.gz
	control_72hrs_rep2
	sample004_1.combined.fastq.gz
	sample004_2.combined.fastq.gz
	control_72hrs_rep3
	sample005_1.combined.fastq.gz
	sample005_2.combined.fastq.gz
	FO47_72hrs_rep1
	sample006_1.combined.fastq.gz
	sample006_2.combined.fastq.gz
	FO47_72hrs_rep2
	sample007_1.combined.fastq.gz
	sample007_2.combined.fastq.gz
	FO47_72hrs_rep3
	sample012_1.combined.fastq.gz
	sample012_2.combined.fastq.gz
	Fus2_0hrs_prelim
	1_S1_L001_R1_001.fastq.gz
	1_S1_L001_R2_001.fastq.gz
	Fus2_16hrs_prelim
	3_S2_L001_R1_001.fastq.gz
	3_S2_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep1
	4_S3_L001_R1_001.fastq.gz
	4_S3_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep2
	Fus2_36hrs_prelim
	36hr-root_S10_L001_R1_001.fastq.gz
	36hr-root_S10_L001_R2_001.fastq.gz
	Fus2_48hrs_prelim
	6_S4_L001_R1_001.fastq.gz
	6_S4_L001_R2_001.fastq.gz
	Fus2_4hrs_prelim
	4hr-root_S7_L001_R1_001.fastq.gz
	4hr-root_S7_L001_R2_001.fastq.gz
	Fus2_72hrs_prelim
	7_S5_L001_R1_001.fastq.gz
	7_S5_L001_R2_001.fastq.gz
	Fus2_72hrs_rep1
	sample016_1.combined.fastq.gz
	sample016_2.combined.fastq.gz
	Fus2_72hrs_rep2
	sample018_1.combined.fastq.gz
	sample018_2.combined.fastq.gz
	Fus2_72hrs_rep3
	sample019_1.combined.fastq.gz
	sample019_2.combined.fastq.gz
	Fus2_8hrs_prelim
	8hr-root_S8_L001_R1_001.fastq.gz
	8hr-root_S8_L001_R2_001.fastq.gz
	Fus2_96hrs_prelim
	8_S6_L001_R1_001.fastq.gz
	8_S6_L001_R2_001.fastq.gz
	Fus2_CzapekDox
	6_S2_L001_R1_001_fastqc  6_S2_L001_R1_001.fastq.gz
	6_S2_L001_R2_001_fastqc  6_S2_L001_R2_001.fastq.gz
	Fus2_GlucosePeptone
	7_S3_L001_R1_001_fastqc  7_S3_L001_R1_001.fastq.gz
	7_S3_L001_R2_001_fastqc  7_S3_L001_R2_001.fastq.gz
	Fus2_PDA
	9_S4_L001_R1_001_fastqc  9_S4_L001_R1_001.fastq.gz
	9_S4_L001_R2_001_fastqc  9_S4_L001_R2_001.fastq.gz
	Fus2_PDB
	4_S1_L001_R1_001_fastqc  4_S1_L001_R1_001.fastq.gz
	4_S1_L001_R2_001_fastqc  4_S1_L001_R2_001.fastq.gz
```

Perform qc of RNAseq timecourse data
```bash
	for num in $(qstat | grep 'rna_qc' | cut -f1 -d ' '); do
		Treatment=$(head -n1 rna_qc_fastq-mcf.sh.o$num | cut -f10 -d '/')
		for FilePath in $(raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1); do
			echo $FilePath
			FileF=$(ls $FilePath/F/*.gz)
			FileR=$(ls $FilePath/R/*.gz)
			IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
			qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
			sleep 10
		done
	done
```
```
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1
Your job 6436212 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep2
Your job 6436213 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep3
Your job 6436214 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep2
Your job 6436215 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep3
Your job 6436216 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep1
Your job 6436217 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep2
Your job 6436218 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep3
Your job 6436219 ("rna_qc_fastq-mcf.sh") has been submitted
```

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls qc_rna/paired/F.oxysporum_fsp_cepae/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
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
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'FOP2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/*); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
```
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/Fo47/*/*_contigs_softmasked.fa | grep 'FOP2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```

Output from stdout included:
```
	Processed 22484 loci.                        [*************************] 100%
	Map Properties:
	     Normalized Map Mass: 50507412.55
	     Raw Map Mass: 50507412.55
	     Fragment Length Distribution: Empirical (learned)
	                   Estimated Mean: 181.98
	                Estimated Std Dev: 78.39
	[13:02:48] Assembling transcripts and estimating abundances.
	Processed 22506 loci.                        [*************************] 100%
```

The Estimated Mean: 181.98 allowed calculation of of the mean insert gap to be
-20bp 182-(97*2) where 97 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2'); do
	# for Assembly in $(ls assembly/merged_canu_spades/*/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta); do
	# for Assembly in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -v -e '_rep'); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			InsertGap='-20'
			InsertStdDev='78'
			Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
			while [ $Jobs -gt 1 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
			done
			printf "\n"
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
		done
	done
```


#### Braker prediction

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
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
	mkdir -p alignment/$Organism/$Strain/concatenated
	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
	alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
	alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
	OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
	AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
	GeneModelName="$Organism"_"$Strain"_braker_new
	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
	qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
	done
```

Fasta and gff files were extracted from Braker1 output.

```bash
	for File in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff); do
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

```bash
	InBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated.bam
	ViewBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_view.bam
	SortBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_sorted
	samtools view -b $InBam > $ViewBam
	samtools sort $ViewBam $SortBam
	samtools index $SortBam.bam
```


Cufflinks was run to compare the predicted genes to assembled transcripts:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		done
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```

<!--
The number of Fo47 genes was determined for comparison to number predicted by Braker (16269):
```bash
	fo47_transcripts=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf
	cat $fo47_transcripts | grep 'gene_id' | cut -f2 -d '"' | sort | uniq | wc -l
	# 18191
```
 -->



## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		mkdir -p $OutDir
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -w -e 'Fus2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
	for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_pacbio/*/augustus.gff3 | grep -e 'Fus2'); do
		Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
		Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
		Assembly=$(ls repeat_masked/$Organism/$Strain/*/"$Strain"_contigs_softmasked.fa)
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
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
		$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

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
for DirPath in $(ls -d gene_pred/codingquary/F.*/*/final | grep -w -e'Fus2'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
<!--
## Suplimenting gene models with known genes

Additional gene models were consrtucted in braker / transfered from other isolates
and were exported to the following locations:

```bash
	Fus2Six9=assembly/spades/F.oxysporum_fsp_cepae/Fus2_edited_v2/filtered_contigs/Fus2_edited_v2_contig_1090_six9.gff
```

These gene models were then edited with nano to give names and IDs to these genes. -->

## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
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
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w 'Fus2'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```


## B) SwissProt
<!--
```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
    mkdir -p $OutDir
    blastp \
    -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
    -query $ProjDir/$Proteome \
    -out $OutDir/swissprot_v2015_10_hits.tbl \
    -evalue 1e-100 \
    -outfmt 6 \
    -num_threads 16 \
    -num_alignments 10
  done
``` -->


```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred/swissprot/$Organism/$Strain
		SwissDbDir=../../uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl | grep -w -e 'Fus2'); do
	# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
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

Convert top blast hits into gff annotations

```bash
	for BlastHitsCsv in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv | grep -w -e 'Fus2'); do
		Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
		Column2="$Strain"_gene_homolog
		NumHits=1
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
	done
```

#### Intersecting blast hits with genes from FoL

```bash
	for HitsGff in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.gff | grep -w -e 'Fus2'); do
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $HitsGff| rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsDir=$(dirname $HitsGff)
		FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
		FoLIntersect=$HitsDir/4287_chromosomal_final_genes_combined_intersect.bed
		bedtools intersect -wo -a $HitsGff -b $FoLGenes > $FoLIntersect
	done
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP
 * D) From ORF fragments - Signal peptide & RxLR motif
 * E) From ORF fragments - Hmm evidence of WY domains
 * F) From ORF fragments - Hmm evidence of RxLR effectors


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
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
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/Fus2); do
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
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
		qsub $ProgDir/submit_TMHMM.sh $Proteome
	done
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```

### C) Identification of MIMP-flanking genes

```bash
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
		Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
		Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
		BrakerGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_CodingQuary.gff3)
		QuaryGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_Braker.gff3)
		OutDir=analysis/mimps/$Organism/$Strain
		mkdir -p "$OutDir"
		echo "$Organism - $Strain"
		ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
		$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
		$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
		echo "The number of mimps identified:"
		cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
		bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		bedtools intersect -u -a $QuaryGff -b $OutDir/"$Strain"_mimps_exp.gff >> $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		echo "The following transcripts intersect mimps:"
		MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
		cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' > $MimpGenesTxt
		cat $MimpGenesTxt | wc -l
		echo ""
	done
```

```
	F.oxysporum_fsp_cepae - Fus2
	The number of mimps identified:
	151
	The following transcripts intersect mimps:
	155
```

```bash
	cat gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp.aa | grep '>' | sed 's/>//g' | cut -f1 -d '-' | sed -r 's/\s//g' > gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp_headers.txt
	cat analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP.txt | grep 'Effector' | cut -f1 | sed 's/|//g' > analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP_headers.txt
	cat gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos.txt | cut -f1 > gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos_headers.txt
```



# 4. Genomic analysis
<!--
## 4.1 Identifcation of protospacers

To facilitate CriprCas editing of the Fusarium oxysporum genome target sites
known as protospacers were identified.

This was done using a published program OPTIMuS as well as a parser script. The
commands to do this were as follows:

```bash
	for GeneSeq in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.codingseq | grep -v 'old'); do
	  Organism=$(echo $GeneSeq | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $GeneSeq | rev | cut -f2 -d '/' | rev)
		echo $Organism
		echo $Strain
	  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
	  OutDir=analysis/protospacers/$Organism/$Strain
	  mkdir -p $OutDir
	  $ProgDir/OPTIMuS_EMR.pl $GeneSeq "threshold" 1 > $OutDir/"$Strain"_protospacer_sites.txt
	  $ProgDir/Optimus2csv.py --inp $OutDir/"$Strain"_protospacer_sites.txt  --out $OutDir/"$Strain"_protospacer_by_gene.csv
	done
``` -->


## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |



<!--
### 4.2.a) Orthology between pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogenic_orthology.md


### 4.2.b) Orthology between non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_non-pathogenic_orthology.md
 -->

### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md


<!-- ### 4.2.d) Orthology between all isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_isolates.md -->


## 5.  Searches

## 5.1.A) Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemblies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2'); do
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv | grep -w 'Fus2'); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=SIX_homolog
		NumHits=1
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	done
```

	The blast hits were summarised in a single table for all the genomes. The top
	identity of the top blast hit in relation to the enquire query sequence was
	presented for each blast hit.

```bash
	OutFile=analysis/blast_homology/Fo_path_genes_CRX_summary.tab
	echo "Organism" > tmp2.tab
	cat analysis/blast_homology/F.proliferatum/A8/A8_Fo_path_genes_CRX.fa_homologs.csv | cut -f1 >> tmp2.tab
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		echo "$Organism" > tmp.tab
		echo "$Strain" >> tmp.tab
		cat $BlastHits | cut -f10  | tail -n +2 >> tmp.tab
		paste tmp2.tab tmp.tab > $OutFile
		cp $OutFile tmp2.tab
	done
	rm tmp.tab
	rm tmp2.tab
```

```bash
	for HitsGff in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.gff | grep -v 'trinity' | grep -w 'Fus2'); do
		Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		# GffBraker=gene_pred/codingquary/$Organism/$Strain/final/final_genes_Braker.gff3
		# GffQuary=gene_pred/codingquary/$Organism/$Strain/final/final_genes_CodingQuary.gff3
		GffAppended=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		OutDir=$(dirname $HitsGff)
		SixIntersect=$OutDir/"$Strain"_Fo_path_genes_CRX.fa_hit_genes.bed
		# bedtools intersect -wo -a $HitsGff -b $GffBraker > $SixIntersect
		# bedtools intersect -wo -a $HitsGff -b $GffQuary >> $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended > $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended | cut -f9,18 | grep -v 'Parent'
		echo ""
	done > analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab
```


## 5.1.B) Identifying FTF genes

Previously published FTF genes from Sanchez et al 2016 were blasted against
Fusarium genomes.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'cepae' | grep -v -e '_edited_v2' -e 'HB17'); do
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
		OutDir=analysis/FTF/$Organism/$Strain
		mkdir -p $OutDir
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
	done
```

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
	for BlastHits in $(ls analysis/FTF/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_hits.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/FTF/$Organism/$Strain
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=FTF_homolog
		NumHits=1
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

		GffAppended=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
		bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_FTF_hits_intersected.bed
	done
```



## 5.2 Identifying PHIbase homologs

The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.
<!--
```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	printf "header\n" > $PhibaseDir/PHI_headers.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/\r//g' >> $PhibaseDir/PHI_headers.csv
	printf "effect\n" > .$PhibaseDir/PHI_virulence.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '|' | rev  >> $PhibaseDir/PHI_virulence.csv
```


```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	PhibaseHeaders=$PhibaseDir/PHI_headers.csv
	PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
	for BlastCSV in $(ls analysis/blast_homology/F*/*/*_PHI_36_accessions.fa_homologs.csv); do
		Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
		echo "$Strain"
		OutDir=$(dirname $BlastCSV)
		paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
		cat $OutDir/"$Strain"_PHIbase_virulence.csv | grep 'NODE_' | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
	done
```
results were:

```
	125
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     66 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    131 Reduced virulence
	     87 Unaffected pathogenicity
	55
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     58 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    121 Reduced virulence
	     83 Unaffected pathogenicity
	A23
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     79 Unaffected pathogenicity
	A28
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     58 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    125 Reduced virulence
	     79 Unaffected pathogenicity
	D2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     62 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     81 Unaffected pathogenicity
	Fus2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     59 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    120 Reduced virulence
	     79 Unaffected pathogenicity
	HB17
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     84 Unaffected pathogenicity
	PG
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     61 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    126 Reduced virulence
	     78 Unaffected pathogenicity
	N139
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     72 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    149 Reduced virulence
	     93 Unaffected pathogenicity
	PG18
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     72 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    153 Reduced virulence
	     97 Unaffected pathogenicity
	PG3
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     65 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    139 Reduced virulence
	     93 Unaffected pathogenicity
	A8
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      4 Increased virulence (Hypervirulence)
	     15 Lethal
	     71 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    152 Reduced virulence
	     82 Unaffected pathogenicity

```
-->

<!-- The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Proteins in $(ls gene_pred/braker/F.oxysporum_fsp_cepae/*/*/*.aa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteins protein ../../phibase/v3.8/PHI_accessions.fa
	done
``` -->

### Expressed Genes

As a preliminary analysis of the RNAseq data, highly expressed genes at 72hrs
post infection were identified in Fus2.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	# samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam
	for AcceptedHits in $(ls alignment/*/*/concatenated/concatenated.bam | grep -v -e 'Fus2_edited_v2'); do
		Strain=$(echo $AcceptedHits | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $AcceptedHits | rev | cut -f4 -d '/' | rev)
		OutDir=$(dirname $AcceptedHits)
		echo "$Organism - $Strain"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

```bash
	for Transcripts in $(ls alignment/F.*/*/concatenated/transcripts.gtf); do
		Strain=$(echo $Transcripts | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Transcripts | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		GeneGff=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		ExpressedGenes=alignment/$Organism/$Strain/concatenated/expressed_genes.bed
		bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
	done
```


<!--
```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=Fus2
	samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam

	cufflinks -o timecourse/2016_genes/Fus2/72hrs/cufflinks -p 16 --max-intron-length 4000 alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam

	Transcripts=timecourse/2016_genes/Fus2/72hrs/cufflinks/transcripts.gtf
	GeneGff=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_appended.gff3
	ExpressedGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes.gff
	bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
```
-->

<!--
It was noted that not all the top expressed transcripts had gene models. The
transcripts without gene models were identified:

```bash
	ExpressedTranscriptsTxt=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | cut -f2 -d'"' | uniq | head -n20 > $ExpressedTranscriptsTxt
	ExpressedTranscriptsInter=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts_intersect.gff
	cat $ExpressedGenes | grep -v -P "\ttranscript_id" | sort -r -n -t'"' -k6 | grep -w -f $ExpressedTranscriptsTxt > $ExpressedTranscriptsInter
	ExpressedTranscriptsGff=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.gff
	cat $ExpressedTranscriptsInter | cut -f1-9 | uniq> $ExpressedTranscriptsGff
	cat $ExpressedTranscriptsInter | grep -v -w 'exon' | cut -f2,18 -d'"' --output-delimiter " - " | uniq | less
```

The top 20 expressed transcripts are shown below with genes they correspond to.
9 of the 20 top expressed transcripts have not been predicted as genes. This was
determined, in large, to be a result of transposon being expressed. It was also
noted that cufflinks doesn't always predict direction of a transcript and
therefore when features are intersected, they can not be restricted by strand.

```
	CUFF.4396 - 33 bp - A-rich
	CUFF.2444 - 69 bp - T-rich - repmasked
	CUFF.540 - 60 bp - C-Rich - repmasked
	CUFF.11337 - 92 bp - T-rich - repmasked
	CUFF.14001 - g11792 - Six5
	CUFF.5007 - repmasked
	CUFF.4291 - repmasked
	CUFF.12587 - g11716 - Six3
	CUFF.4090 - repmasked
	CUFF.12214 - g4762
	CUFF.8542 - 103 bp - A-rich - repmasked
	CUFF.13998 - g12630
	CUFF.13635 - g10978
	CUFF.4290 - no gene - Six9
	CUFF.2922 - g12344
	CUFF.9876 - g10974
	CUFF.13435 - g10975
	CUFF.9559 - g10973
	CUFF.12530 - g4720
	CUFF.6539 - g10346
```
-->


<!--
```bash
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | grep -w 'gene' | head -n 20 | cut -f18 > $TopGenes
	InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
	TopFunctions=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.tsv
	cat $InterPro | grep -w -f $TopGenes > $TopFunctions

	IsolateAbrv=FoC_path_vs_non_path
	WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
	Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
	PathDir=$WorkDir/path_orthogroups
	echo "The number of pathogen common genes in the top 20 expressed Fus2 genes at 72hr"
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $TopGenes | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt
	# secreted
	PathOrthogroupsFus2Secreted=$PathDir/Fus2_path_orthogroup_secreted.txt
	cat $PathOrthogroupsFus2Secreted | grep -w -f $TopGenes
	# near mimps
	MimpGenesTxt=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_genes_in_2kb_mimp.txt
  cat $MimpGenesTxt | grep -w -f $TopGenes
	# Intersecting SIX genes

``` -->


## 6. Summarising the Fusarium Proteome

```bash

```


```
