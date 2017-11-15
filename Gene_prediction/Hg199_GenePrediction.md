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

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
for BrakerGff in $(ls gene_pred/braker/N.*/Hg199_braker_first/*/augustus.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_first//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
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

gene_pred/codingquary/N.ditissima/Hg199/final
12807
954
13761

gene_pred/codingquary/N.ditissima/R0905_pacbio_canu/final
12946
53
12999

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
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/N.*/*/*/*_contigs_unmasked.fa); do
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
	for Genes in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
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
	for Proteome in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
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
	# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
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
	SplitfileDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
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
	for SplitDir in $(ls -d gene_pred/final_genes_split/N.*/*); do
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
	for Proteome in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
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
	for Proteome in $(ls gene_pred/codingquary/N.*/Hg199/*/final_genes_combined.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```
