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
	for Genome in $(ls repeat_masked/N.*/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```
** Number of cegma genes present and complete: 95.56%
** Number of cegma genes present and partial: 96.77%

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/N.*/Hg199/*_dna_cegma.completeness_report); do
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
  for FilePath in $(ls -d raw_rna/paired/N.*/Hg199); do
    echo $FilePath;
    FileF=$(ls $FilePath/F/*.gz);
    FileR=$(ls $FilePath/R/*.gz);
    IlluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa; ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc;
    qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA;
  done
```

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls qc_rna/paired/N.ditissima/Hg199/*/*.fq.gz); do
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
	for Assembly in $(ls repeat_masked/*/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199); do
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
79.1% overall read mapping rate.
70.2% concordant pair alignment rate.

Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	# AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
  AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
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
  Processed 16528 loci.                        [*************************] 100%
  Map Properties:
 	Normalized Map Mass: 13483696.89
	Raw Map Mass: 13483696.89
	Fragment Length Distribution: Empirical (learned)
	              Estimated Mean: 217.61
	           Estimated Std Dev: 43.68
[09:20:01] Assembling transcripts and estimating abundances.
  Processed 16573 loci.                        [*************************] 100%
```

The Estimated Mean: 217.61 allowed calculation of of the mean insert gap to be
-123bp 217-(170*2) where 170? was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.

Then Rnaseq data was aligned to each genome assembly:

```bash
		for Assembly in $(ls repeat_masked/*/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
			Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
			Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
			echo "$Organism - $Strain"
			for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199 | grep -v -e '_rep'); do
				Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
				echo "$Timepoint"
				FileF=$(ls $RNADir/F/*_trim.fq.gz)
				FileR=$(ls $RNADir/R/*_trim.fq.gz)
				OutDir=alignment/$Organism/$Strain/$Timepoint
				InsertGap='-123'
				InsertStdDev='44'
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
79.1% overall read mapping rate.
70.2% concordant pair alignment rate.

#### Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:
```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash
#for Assembly in $(ls repeat_masked/N.ditissima/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
#Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
#while [ $Jobs -gt 1 ]; do
#sleep 10
#printf "."
#Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
#done
#printf "\n"
#Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
#Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
#echo "$Organism - $Strain"
#mkdir -p alignment/$Organism/$Strain/concatenated
#samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
#alignment/$Organism/Hg199/Hg199/accepted_hits.bam
#OutDir=gene_pred/braker/$Organism/"$Strain"_braker_first
#AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
#GeneModelName="$Organism"_"$Strain"_braker_first
#rm -r /home/gomeza/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_first
#ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
#qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
#done
```


```bash
for Assembly in $(ls repeat_masked/N.ditissima/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
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
OutDir=gene_pred/braker/$Organism/"$Strain"_braker_first
AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
GeneModelName="$Organism"_"$Strain"_braker_first
rm -r /home/gomeza/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_first
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/N.*/Hg199_braker_first/*/augustus.gff); do
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
    for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/*/accepted_hits.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
    done
```

Secondly, genes were predicted using CodingQuary:

```bash
		for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
		done
```
