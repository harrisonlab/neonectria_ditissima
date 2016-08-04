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
