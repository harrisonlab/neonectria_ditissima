# R0905 isolate. Neonectria Reference Genome Assembly

==================================================================
# Gene  Prediction
==================================================================

Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Busco
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
	for Genome in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*unmasked.fa); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
	for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*unmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
		BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
		OutDir=$(dirname $Assembly)/busco
		qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
	done
```

## Gene model training

Gene prediction was performed for Neonectria genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1 and prediction of all putative ORFs in the genome using the ORF finder (atg.pl) approach.

####Align mycelium reads to the reference assembly with STAR

```bash
for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*unmasked.fa)
do
    Strain=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/paired/N.ditissima/R0905/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
        echo $FileF
        echo $FileR
        Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
        echo "$Timepoint"
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment_vAG/star/$Organism/$Strain
        ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq/
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done

for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*unmasked.fa)
do
    Strain=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/rna_seq/Hg199/mycelium/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
        echo $FileF
        echo $FileR
        Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
        echo "$Timepoint"
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment_vAG/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq/
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done
```
### Alignment outputs were concatenated and Braker1 prediction was run

```bash
  Strain=R0905
  Organism=N.ditissima
  mkdir -p alignment_vAG/star/$Organism/$Strain/mycelium/concatenated
  samtools merge -f alignment_vAG/star/$Organism/$Strain/mycelium/concatenated/concatenated.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_1/star_aligmentAligned.sortedByCoord.out.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_2/star_aligmentAligned.sortedByCoord.out.bam \
  alignment_vAG/star/$Organism/$Strain/mycelium/Hg199_3/star_aligmentAligned.sortedByCoord.out.bam
```

## Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/2018/gm_key_64 ~/.gm_key
```

```bash
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred_vAG/braker/Ref_Genomes/$Organism/$Strain
    AcceptedHits=alignment_vAG/star/N.ditissima/Ref_Genomes/R0905/Hg199_reads/mycelium/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    #rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```

Fasta and gff files were extracted from Braker1 output.

```bash
  for File in $(ls gene_pred_vAG/braker/Ref_Genomes/N.*/R0905/N.ditissima_R0905_braker/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```

## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
    Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred_vAG/cufflinks/Ref_Genomes/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=alignment_vAG/star/N.ditissima/Ref_Genomes/R0905/Hg199_reads/mycelium/concatenated/concatenated.bam
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_unmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
	  Organism=$(echo $Assembly| rev | cut -d '/' -f6 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/
		mkdir -p $OutDir
		CufflinksGTF=gene_pred_vAG/cufflinks/Ref_Genomes/$Organism/$Strain/concatenated_prelim/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note: This needs to run step by step.

```bash
  BrakerGff=$(ls gene_pred_vAG/braker/Ref_Genomes/*/R0905/N.ditissima_R0905_braker/augustus.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuaryGff=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/out/PredictedPass.gff3
	PGNGff=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/out/PGN_predictedPass.gff3
	AddDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/additional
	FinalDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final
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
	for DirPath in $(ls -d gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/final); do
	echo $DirPath;
	cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
	cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
	cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
	echo "";
	done
```

```
gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/Hg199/final
13086
1752
14838

gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/R0905/final
12883
2272
15155
```

Remove duplicate genes

```bash
for GffAppended in $(ls gene_pred_vAG/codingquary/Ref_Genomes/*/R0905/final/final_genes_appended.gff3);
do
  Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
  Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  FinalDir=gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final
  GffFiltered=$FinalDir/filtered_duplicates.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary/
  $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
	GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
	LogFile=$FinalDir/final_genes_appended_renamed.log
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/codingquary
	$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
	rm $GffFiltered
	Assembly=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed
	# The proteins fasta file contains * instead of Xs for stop codons, these should
	# be changed
	sed -i 's/\*/X/g' gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```

No duplicated genes

==================================================================
# ORF Finder
==================================================================

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. Additional functionality was added to this script by
also printing ORFs in .gff format.

Script needs to be edited
```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  for Genome in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked/*_contigs_unmasked.fa); do
  	qsub $ProgDir/run_ORF_finder.sh $Genome
  done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for ORF_Gff in $(ls gene_pred_vAG/ORF_finder/N.ditissima/Hg199/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_corrected.gff3/g')
		echo ""
		echo "Correcting the following file:"
		echo $ORF_Gff
		echo "Redirecting to:"
		echo $ORF_Gff_mod
		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
	done
```

==================================================================
# Functional annotation
==================================================================

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
	for Genes in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	echo $Genes
  	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteins in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
  	Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
  	echo "$Organism - $Strain"
  	echo $Strain
  	InterProRaw=gene_pred_vAG/interproscan/Ref_Genomes/$Organism/$Strain/raw
  	$ProgDir/append_interpro.sh $Proteins $InterProRaw
  done
```
## B) SwissProt

```bash
	for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred_vAG/swissprot/Ref_Genomes/$Organism/$Strain
		SwissDbDir=../../../../home/groups/harrisonlab/uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
for SwissTable in $(ls gene_pred_vAG/swissprot/Ref_Genomes/*/R0905/swissprot_vMar2018_10_hits.tbl); do
	Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
	Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
	echo "$Organism - $Strain"
	OutTable=gene_pred_vAG/swissprot/Ref_Genomes/$Organism/$Strain/swissprot_vMar2018_tophit_parsed.tbl
	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
	$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta > $OutTable
done
```

==================================================================
# Effector genes
==================================================================

## A) From Augustus gene models - Identifying secreted proteins

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
  for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred_vAG/final_genes_split/Ref_Genomes/$Organism/$Strain
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
      #qsub $ProgDir/pred_sigP.sh $File
      #qsub $ProgDir/pred_sigP.sh $File signalp-3.0
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
  for SplitDir in $(ls -d gene_pred_vAG/final_genes_split/Ref_Genomes/N.*/R0905); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    for SigpDir in $(ls -d gene_pred_vAG/Ref_Genomes_signalp-4.1 | cut -f2 -d'/')
    do
      InStringAA=''
      InStringNeg=''
      InStringTab=''
      InStringTxt=''
      for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
        InStringAA="$InStringAA gene_pred_vAG/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
        InStringNeg="$InStringNeg gene_pred_vAG/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
        InStringTab="$InStringTab gene_pred_vAG/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
        InStringTxt="$InStringTxt gene_pred_vAG/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
      done
      cat $InStringAA > gene_pred_vAG/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
      cat $InStringNeg > gene_pred_vAG/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
      tail -n +2 -q $InStringTab > gene_pred_vAG/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
      cat $InStringTxt > gene_pred_vAG/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
    done
  done
```

SigP v2 did not predict any gene. Only v4 prediction will be used.

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

 ```bash
  for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
  	Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  	Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
  	qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
 ```
Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
  for File in $(ls gene_pred_vAG/trans_mem/Ref_Genomes/*/R0905/*_TM_genes_neg.txt); do
   Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
   cat $File | cut -f1 > $TmHeaders
   SigP=$(ls gene_pred_vAG/Ref_Genomes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
   OutDir=$(dirname $SigP)
   ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
   cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
  done
 ```
 ```
 N.ditissima - Hg199
1032
 ```

## B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
	for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis_vAG/effectorP/Ref_Genomes/$Organism/$Strain
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```

```bash
  for File in $(ls analysis_vAG/effectorP/Ref_Genomes/*/R0905/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    Secretome=$(ls gene_pred_vAG/Ref_Genomes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
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
  for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred_vAG/CAZY/Ref_Genomes/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done


New database
  for Proteome in $(ls gene_pred_vAG/codingquary/Ref_Genomes/N.*/R0905/*/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred_vAG/CAZY_NewDatabase/Ref_Genomes/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=dbCAN/dbCAN-HMMdb-V7.txt
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
  for File in $(ls gene_pred_vAG/CAZY/Ref_Genomes/N.*/R0905/*CAZY.out.dm); do
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
  Gff=$(ls gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

  SecretedProts=$(ls gene_pred_vAG/Ref_Genomes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
  echo "number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
```
N.ditissima - R0905
number of CAZY genes identified:
723
number of Secreted CAZY genes identified:
291

```bash
  for File in $(ls gene_pred_vAG/CAZY_NewDatabase/Ref_Genomes/N.*/R0905/*CAZY.out.dm); do
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
  Gff=$(ls gene_pred_vAG/codingquary/Ref_Genomes/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

  SecretedProts=$(ls gene_pred_vAG/Ref_Genomes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
  echo "number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
```
N.ditissima - R0905
number of CAZY genes identified:
737
number of Secreted CAZY genes identified:
296

==================================================================
# PhiBase genes:Identifying PHIbase homologs
==================================================================

```bash
cd /data/scratch/gomeza
dbFasta=$(ls /home/groups/harrisonlab/phibase/v4.5/phi_accessions.fa)
dbType="prot"
QueryFasta=$(ls gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/R0905/final/final_genes_appended_renamed.cdna.fasta)
Prefix="R0905_phi_accessions"
Eval="1e-30"
OutDir=analysis_vAG/blast_homology/Ref_Genomes/N.ditissima/R0905
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

==================================================================
# Looking for Transcription Factors
==================================================================

```bash
  for Interpro in $(ls gene_pred_vAG/interproscan/Ref_Genomes/N.ditissima/Hg199/*_interproscan.tsv); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis_vAG/transcription_factors/Ref_Genomes/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
  done
```

N.ditissima - Hg199
total number of transcription factors
934

Interproscan results use transcript IDs rather gene IDs. I produce an additional file with gene names

```bash
cat analysis_vAG/transcription_factors/Ref_Genomes/N.ditissima/Hg199/Hg199_TF_gene_headers.txt | sed -e "s/.t.*//g" > analysis_vAG/transcription_factors/Ref_Genomes/N.ditissima/Hg199/Hg199_TF_geneid_headers.txt
```

==================================================================
# Secondary metabolites (Antismash and SMURF)
==================================================================

Antismash was run to identify clusters of secondary metabolite genes within the genome. Antismash was run using the weserver at: http://antismash.secondarymetabolites.org

Results of web-annotation of gene clusters within the assembly were downloaded to the following directories:

```bash
for AntiSmash in $(ls analysis_vAG/secondary_metabolites/antismash/Ref_Genomes/N.ditissima/Hg199/fungi-2ab971b9-127d-4760-94fc-96f18c99b501/*.gbk); do
  OutDir=analysis_vAG/secondary_metabolites/antismash/Ref_Genomes/N.ditissima/Hg199
  Prefix=$OutDir/Hg199
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix
done

# Identify secondary metabolites within predicted clusters
printf "Number of secondary metabolite detected:\t"
cat "$Prefix"_secmet_clusters.gff | wc -l
GeneGff=gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3
bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_secmet_genes.txt
bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
printf "Number of predicted proteins in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.txt | wc -l
printf "Number of predicted genes in secondary metabolite clusters:\t"
cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

# Identify cluster finder additional non-secondary metabolite clusters. Not supported in the Antismash v5.
printf "Number of cluster finder non-SecMet clusters detected:\t"
cat "$Prefix"_clusterfinder_clusters.gff | wc -l
GeneGff=gene_pred_vAG/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3
bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt
printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
cat "$Prefix"_clusterfinder_genes.txt | wc -l
printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
done

These clusters represented the following genes. Note that these numbers just show the number of intersected genes with gff clusters and are not confirmed by function

Number of secondary metabolite detected:	40
Number of predicted proteins in secondary metabolite clusters:	579
Number of predicted genes in secondary metabolite clusters:	563
```
