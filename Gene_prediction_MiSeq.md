

Busco has replaced CEGMA and was run to check gene space in assemblies

Previous isolates

```bash
for Assembly in $(ls repeat_masked/N.ditissima/*/*unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls repeat_masked/N.ditissima/*/*/*unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
short_summary_Hg199_contigs_unmasked.txt
3673    Complete BUSCOs (C) 3658    Complete and single-copy BUSCOs (S)
15 Complete and duplicated BUSCOs (D) 24      Fragmented BUSCOs (F)
28 Missing BUSCOs (M) 3725    Total BUSCO groups searched
short_summary_AgN04_contigs_unmasked.txt	3663	17	23	39	3725
short_summary_R0905_contigs_unmasked.txt	3651	16	20	54	3725
short_summary_R45-15_contigs_unmasked.txt	3674	15	22	29	3725
short_summary_R0905_contigs_unmasked.txt	3668	16	21	36	3725 canu_2017

Gene prediction

```bash
for Assembly in $(ls repeat_masked/*/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
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
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
screen -a
for Assembly in $(ls repeat_masked/*/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
# AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
OutDir=gene_pred/cufflinks/$Organism/$Strain/vs_Hg199reads/concatenated_prelim
echo "$Organism - $Strain"
mkdir -p $OutDir
cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
```
I have aligned every isolate (Ag04, R45/15, Hg199 and R0905_canu_2017_v2) with the Hg199 RNA reads.

Output from stdout included:
```
N.ditissima - AgN04
tee: gene_pred/cufflinks/N.ditissima/AgN04/vs_Hg199reads/concatenated_prelim/cufflinks/cufflinks.log: No such file or directory
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[09:54:31] Inspecting reads and determining fragment length distribution.
> Processed 18766 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 12803452.83
>	Raw Map Mass: 12803452.83
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 217.01
>	           Estimated Std Dev: 43.68
[09:56:57] Assembling transcripts and estimating abundances.
> Processed 18835 loci.                        [*************************] 100%
N.ditissima - R45-15
tee: gene_pred/cufflinks/N.ditissima/R45-15/vs_Hg199reads/concatenated_prelim/cufflinks/cufflinks.log: No such file or directory
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[10:04:15] Inspecting reads and determining fragment length distribution.
> Processed 18553 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 12833084.96
>	Raw Map Mass: 12833084.96
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 217.19
>	           Estimated Std Dev: 43.63
[10:06:35] Assembling transcripts and estimating abundances.
> Processed 18638 loci.                        [*************************] 100%
N.ditissima - Hg199
tee: gene_pred/cufflinks/N.ditissima/Hg199/vs_Hg199reads/concatenated_prelim/cufflinks/cufflinks.log: No such file or directory
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[10:48:22] Inspecting reads and determining fragment length distribution.
> Processed 16497 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 13479906.07
>	Raw Map Mass: 13479906.07
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 217.62
>	           Estimated Std Dev: 43.68
[10:51:11] Assembling transcripts and estimating abundances.
> Processed 16537 loci.                        [*************************] 100%
N.ditissima - R0905_canu_2017_v2
Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
[14:23:01] Inspecting reads and determining fragment length distribution.
> Processed 18535 loci.                        [*************************] 100%
> Map Properties:
>	Normalized Map Mass: 12759384.33
>	Raw Map Mass: 12759384.33
>	Fragment Length Distribution: Empirical (learned)
>	              Estimated Mean: 217.05
>	           Estimated Std Dev: 43.58
[14:25:44] Assembling transcripts and estimating abundances.
> Processed 18606 loci.                        [*************************] 100%
```

The Estimated Mean: 219.68 allowed calculation of of the mean insert gap to be
-140bp 182-(180*2) where 180? was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/N*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
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
  for Assembly in $(ls repeat_masked/N.ditissima/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
    done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker_Nov2017
    AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker_Nov2017
    rm -r /home/gomeza/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_Nov2017
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

Fasta and gff files were extracted from Braker1 output.

```bash
  for File in $(ls gene_pred/braker/N.*/*Nov2017/*/augustus.gff); do
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
    for Assembly in $(ls repeat_masked/N*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
    done

    for Assembly in $(ls repeat_masked/N*/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/Hg199/accepted_hits.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
    done
```

Secondly, genes were predicted using CodingQuary:

```bash
    for Assembly in $(ls repeat_masked/N*/Hg199/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/NewHg199
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
    done

    for Assembly in $(ls repeat_masked/N*/*/R0905_canu*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
    done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/$Organism/Hg199/*/augustus.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
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


gene_pred/codingquary/N.ditissima/Ag04/final
12805
854
13659

gene_pred/codingquary/N.ditissima/Hg199/final
12818
905
13723

gene_pred/codingquary/N.ditissima/R0905_canu_2017_v2/final
12848
898
13746

gene_pred/codingquary/N.ditissima/R45-15/final
12635
989
13624



## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
for Genome in $(ls repeat_masked/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
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
	for SwissTable in $(ls gene_pred/swissprot/*/*/*_hits.tbl); do
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2017_tophit_parsed.tbl
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
  screen -a
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
 	for SplitDir in $(ls -d gene_pred/final_genes_split/N.*/*_201*); do
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
 	for Proteome in $(ls gene_pred/codingquary/N.*/*_201*/*/final_genes_combined.pep.fasta); do
 		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
 		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
 		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
 		qsub $ProgDir/submit_TMHMM.sh $Proteome
 	done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
   for File in $(ls gene_pred/trans_mem/*/*_201*/*_TM_genes_neg.txt); do
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
    N.ditissima - R0905_canu_2017_v2
    1017
    N.ditissima - R0905_merged_2017
    991

### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
	for Proteome in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```

```bash
  for File in $(ls analysis/effectorP/*/R0905_canu_2017_v2/*_EffectorP.txt); do
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

### C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
  for Proteome in $(ls gene_pred/codingquary/N.*/*/*/final_genes_combined.pep.fasta); do
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
  for File in $(ls gene_pred/CAZY/N.*/*_201*/*CAZY.out.dm); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $File)
  echo "$Organism - $Strain"
  ProgDir=/home/groups/harrisonlab/dbCAN
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

  N.ditissima - R0905_canu_2017_v2
  number of CAZY genes identified:
  721
  number of Secreted CAZY genes identified:
  287
  N.ditissima - R0905_merged_2017
  number of CAZY genes identified:
  720
  number of Secreted CAZY genes identified:
  283

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

### Alignment of raw reads vs the Fus2 genome

Sequence data for isolates with a data from a single sequencing run was aligned against the R0905 genome

```bash
  Reference=$(ls repeat_masked/N.ditissima/Ref_Genomes/R0905_canu_2017_v2/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/new_isolates_2017/*); do
    Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_canu_2017/
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
  ```
