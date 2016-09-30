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
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/*); do
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
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -w -e 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done
```

## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final_genes/N.*/*/*/final_genes_combined.pep.fasta); do
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
	for File in $(ls gene_pred/CAZY/F.*/*/*CAZY.out.dm); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		OutDir=$(dirname $File)
		echo "$Organism - $Strain"
		ProgDir=/home/groups/harrisonlab/dbCAN
		$ProgDir/hmmscan-parser.sh $OutDir/Fus2_canu_new_CAZY.out.dm > $OutDir/Fus2_canu_new_CAZY.out.dm.ps
		SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
		SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
		cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
		Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
		CazyGff=$OutDir/Fus2_canu_new_CAZY.gff
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
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
