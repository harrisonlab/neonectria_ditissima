=================================
# Gene  Prediction
=================================

Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Busco
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

```bash
	for Assembly in $(ls repeat_masked/Ref_Genomes/*/*/*/*unmasked.fa); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
		BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
		OutDir=gene_pred/busco/Ref_Genomes/$Organism/$Strain/CSAR_assembly
		qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
	done
```

## Gene model training

```bash
for Strain in Hg199 R0905; do
	for Assembly in $(ls repeat_masked/Ref_Genomes/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
	  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	  echo "$Organism - $Strain"
	  for RNADir in $(ls -d qc_rna/paired/N.ditissima/$Strain); do
	    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
	    echo "$Timepoint"
	    FileF=$(ls $RNADir/F/*_trim.fq.gz)
	    FileR=$(ls $RNADir/R/*_trim.fq.gz)
	    OutDir=alignment/Ref_Genomes/$Organism/$Strain/$Timepoint
	    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
	    qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
done
```

Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

screen -a

```bash
for Strain in Hg199 R0905; do
	for Assembly in $(ls repeat_masked/Ref_Genomes/N.*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	AcceptedHits=alignment/Ref_Genomes/$Organism/$Strain/$Strain/accepted_hits.bam
	OutDir=gene_pred/cufflinks/Ref_Genomes/$Organism/$Strain/vs_"$Strain"reads/concatenated_prelim
	echo "$Organism - $Strain"
	mkdir -p $OutDir
	cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
done
```
I have aligned every isolate with their RNA reads.

Output from stdout included:

```
N.ditissima - Hg199
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[16:36:09] Inspecting reads and determining fragment length distribution.
> Processed 16735 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 13257861.80
>       Raw Map Mass: 13257861.80
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 217.50
>                  Estimated Std Dev: 43.70
[16:38:43] Assembling transcripts and estimating abundances.
> Processed 16782 loci.                        [*************************] 100%
N.ditissima - R0905
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[16:49:16] Inspecting reads and determining fragment length distribution.
> Processed 19118 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 12100829.47
>       Raw Map Mass: 12100829.47
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 219.53
>                  Estimated Std Dev: 39.59
[16:52:14] Assembling transcripts and estimating abundances.
> Processed 19196 loci.                        [*************************] 100%

The Estimated Mean: 219.68 allowed calculation of of the mean insert gap to be
-140bp 182-(180*2) where 180? was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.
```

Then Rnaseq data was aligned to each genome assembly:

```bash
	for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		cho "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199 | grep -v -e '_rep'); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			InsertGap='-140'
			InsertStdDev='43'
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

	for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	cho "$Organism - $Strain"
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

## Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash
  for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
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
    OutDir=gene_pred/braker/Ref_Genomes/$Organism/$Strain
    AcceptedHits=alignment/Ref_Genomes/N.ditissima/Hg199/Hg199/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```
```bash
  for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
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
    OutDir=gene_pred/braker/Ref_Genomes/$Organism/$Strain
    AcceptedHits=alignment/Ref_Genomes/N.ditissima/R0905/R0905/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

Fasta and gff files were extracted from Braker1 output.

```bash
  for File in $(ls gene_pred/braker/Ref_Genomes/N.*/*/*/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```

## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/*/*/*_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  OutDir=gene_pred/cufflinks/Ref_Genomes/$Organism/$Strain/concatenated_prelim
  mkdir -p $OutDir
  AcceptedHits=alignment/Ref_Genomes/N.ditissima/$Strain/$Strain/accepted_hits.bam
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
  qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/Ref_Genomes/N.ditissima/*/*/*_contigs_unmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	  Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/
		mkdir -p $OutDir
		CufflinksGTF=gene_pred/cufflinks/Ref_Genomes/$Organism/$Strain/concatenated_prelim/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note: This needs to run step by step.

```bash
  BrakerGff=$(ls gene_pred/braker/Ref_Genomes/*/R0905/*/augustus.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuaryGff=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/out/PredictedPass.gff3
	PGNGff=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/out/PGN_predictedPass.gff3
	AddDir=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/additional
	FinalDir=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/final
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
```

The final number of genes per isolate was observed using:
```bash
	for DirPath in $(ls -d gene_pred/codingquary/Ref_Genomes/N.*/*/final); do
	echo $DirPath;
	cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
	cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
	cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
	echo "";
	done
```

gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final
13665
1418
15083

gene_pred/codingquary/Ref_Genomes/N.ditissima/R0905/final
13300
1134
14434

Remove duplicate genes

```bash
for GffAppended in $(ls gene_pred/codingquary/Ref_Genomes/*/*/final/final_genes_appended.gff3);
do
  Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
  Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  FinalDir=gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/final
  GffFiltered=$FinalDir/filtered_duplicates.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary/
  $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
	GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
	LogFile=$FinalDir/final_genes_appended_renamed.log
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary
	$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
	rm $GffFiltered
	Assembly=$(ls repeat_masked/Ref_Genomes/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
	$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed
	# The proteins fasta file contains * instead of Xs for stop codons, these should
	# be changed
	sed -i 's/\*/X/g' gene_pred/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```

No duplicated genes were detected.


=================================
# ORF Finder
=================================

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
for Genome in $(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/*/*_contigs_unmasked.fa); do
	qsub $ProgDir/run_ORF_finder.sh $Genome
done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation
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
qalter -l h=blacklace09.blacklace

=================================
# Functional annotation
=================================

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interproscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/Ref_Genomes/N.*/*/*/final_genes_appended_renamed.pep.fasta); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/Ref_Genomes/N.*/*/*/final_genes_appended_renamed.pep.fasta); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/Ref_Genomes/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```

## B) SwissProt

```bash
	for Proteome in $(ls gene_pred/codingquary/Ref_Genomes/N.*/*/*/final_genes_appended_renamed.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred/swissprot/Ref_Genomes/$Organism/$Strain
		SwissDbDir=../../../../home/groups/harrisonlab/uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_vJul2016_10_hits.tbl); do
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```

=================================
# Effector genes
=================================

## A) From Augustus gene models - Identifying secreted proteins

 Required programs:
  * SignalP-4.1
  * TMHMM

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
  screen -a
	for Strain in Hg199 R0905; do
  SplitfileDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  CurPath=$PWD
  for Proteome in $(ls gene_pred/codingquary/Ref_Genomes/N.*/$Strain/*/final_genes_appended_renamed.pep.fasta); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  SplitDir=gene_pred/final_genes_split/Ref_Genomes/$Organism/$Strain
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
  qsub $ProgDir/pred_sigP.sh $File
  qsub $ProgDir/pred_sigP.sh $File signalp-3.0
  qsub $ProgDir/pred_sigP.sh $File signalp-4.1
  done
  done
  done
```
 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
for SplitDir in $(ls -d gene_pred/final_genes_split/N.*/Hg199_minion); do
Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
for SigpDir in $(ls -d gene_pred/final_genes_signalp-4.1 | cut -f2 -d'/')
do
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
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
done
```

 Some proteins that are incorporated into the cell membrane require secretion.
 Therefore proteins with a transmembrane domain are not likely to represent
 cytoplasmic or apoplastic effectors.

 Proteins containing a transmembrane domain were identified:

 ```bash
 	for Proteome in $(ls gene_pred/codingquary/N.*/Hg199_minion/*/final_genes_combined.pep.fasta); do
 		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
 		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
 		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
 		qsub $ProgDir/submit_TMHMM.sh $Proteome
 	done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
for File in $(ls gene_pred/trans_mem/*/Hg199_minion/*_TM_genes_neg.txt); do
 Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
 Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
 echo "$Organism - $Strain"
 TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
 cat $File | cut -f1 > $TmHeaders
 SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
 OutDir=$(dirname $SigP)
 ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
 $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
 cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
 ```

 N.ditissima - Hg199_minion
1033

## B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
cp gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.pep.fasta prog/EffectorP/EffectorP_1.0/Scripts/
cd prog/EffectorP/EffectorP_1.0/Scripts/
python EffectorP.py -i final_genes_combined.pep.fasta -o N.ditissima_Hg199_minion_EffectorP.txt
mv N.ditissima_Hg199_minion_EffectorP.txt ../../../../analysis/effectoP/N.ditissima/Hg199_minion
```

Number of proteins that were tested: 14919
Number of predicted effectors: 2777

-----------------
18.6 percent are predicted to be effectors.

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_appended.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
    cat $EffectorP_Gff | grep -w 'gene' | wc -l
  done > tmp.txt
```

## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
  for Proteome in $(ls gene_pred/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/Ref_Genomes/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

  ```bash
  for File in $(ls gene_pred/CAZY/N.*/*/*CAZY.out.dm); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $File)
  echo "$Organism - $Strain"
  ProgDir=dbCAN
  $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
  CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
  cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
  echo "number of CAZY genes identified:"
  cat $CazyHeaders | wc -l
  Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

  SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
  echo "number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
```
N.ditissima - Hg199_minion
number of CAZY genes identified:
753
number of Secreted CAZY genes identified:
290


=================================
# PhiBase genes:Identifying PHIbase homologs
=================================

```bash
cd /data/scratch/gomeza
dbFasta=$(ls /home/groups/harrisonlab/phibase/v4.4/phi_accessions.fa)
dbType="prot"
QueryFasta=$(ls gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta)
Prefix="Hg199_phi_accessions"
Eval="1e-30"
 # WorkDir=$TMPDIR
 OutDir=analysis/blast_homology/N.ditissima/Hg199
 mkdir -p $OutDir
 #-------------------------------------------------------
 # 		Step 1.		Make database
 #-------------------------------------------------------
 makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
 #-------------------------------------------------------
 # 		Step 2.		Blast Search
 #-------------------------------------------------------
 blastx -num_threads 4 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
 #-------------------------------------------------------
 # 		Step 3.		Summarise hits
 #-------------------------------------------------------
 cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
```

=================================
# Looking for Transcription Factors
=================================

```bash
for Interpro in $(ls /data/scratch/gomeza/gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv); do
  Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/transcription_factors/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
  $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
  echo "total number of transcription factors"
  cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
  cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
done
```
N.ditissima - Hg199_minion
total number of transcription factors
967

Interproscan results use transcript IDs rather gene IDs. I produce an additional file with gene names

```bash
cat analysis/transcription_factors/N.ditissima/Hg199_minion/Hg199_minion_TF_gene_headers.txt | sed -e "s/.t.*//g" > analysis/transcription_factors/N.ditissima/Hg199_minion/Hg199_minion_TF_geneid_headers.txt
```

=================================

Secondary metabolites (Antismash and SMURF)
=================================

Antismash was run to identify clusters of secondary metabolite genes within the genome. Antismash was run using the weserver at: http://antismash.secondarymetabolites.org

Results of web-annotation of gene clusters within the assembly were downloaded to the following directories:

```bash
for AntiSmash in $(ls analysis/secondary_metabolites/antismash/Hg199_minion/fungi-dfc734cf-18aa-414d-b034-0da05c627613/*.gbk); do
  OutDir=analysis/secondary_metabolites/antismash/Hg199_minion/
  Prefix=$OutDir/Hg199_antismash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix
done

# Identify secondary metabolites within predicted clusters
printf "Number of secondary metabolite detected:\t"
cat "$Prefix"_secmet_clusters.gff | wc -l
GeneGff=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3
bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_secmet_genes.txt
bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
printf "Number of predicted proteins in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.txt | wc -l
printf "Number of predicted genes in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

# Identify cluster finder additional non-secondary metabolite clusters
printf "Number of cluster finder non-SecMet clusters detected:\t"
cat "$Prefix"_clusterfinder_clusters.gff | wc -l
GeneGff=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3
bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt
printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
cat "$Prefix"_clusterfinder_genes.txt | wc -l
printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
done

These clusters represented the following genes. Note that these numbers just show the number of intersected genes with gff clusters and are not confirmed by function

Number of secondary metabolite detected:	45
Number of predicted proteins in secondary metabolite clusters:	392
Number of predicted genes in secondary metabolite clusters:	380
Number of cluster finder non-SecMet clusters detected:	110
Number of predicted proteins in cluster finder non-SecMet clusters:	1339
Number of predicted genes in cluster finder non-SecMet clusters:	1321
```
