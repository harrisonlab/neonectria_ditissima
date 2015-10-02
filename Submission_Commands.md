# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the N. ditissima genome for
submisson are shown below.


# Final Submission

These commands were used in the final submission of the N. ditissima genome:


## Output directory
An output and working directory was made for genome submission:

```bash
  	cd /home/groups/harrisonlab/project_files/neonectria_ditissima
  	OutDir="genome_submission/N.ditissima/R0905"
  	mkdir -p $OutDir
  	ProjDir=/home/groups/harrisonlab/project_files/neonectria_ditissima
```

## SbtFile
The genbank submission template tool was used at:
http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
This produce a template file detailing the submission.

## Setting varibales
Vairables containing locations of files and options for scripts were set:

```bash
  	# Program locations:
  	AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
  	ProgDir="/home/armita/git_repos/emr_repos/tools/genbank_submission"
  	# File locations:
  	SbtFile="/home/groups/harrisonlab/project_files/neonectria_ditissima/collaboration/genome_submission/tbl2asn_out/genome.sbt"
  	Assembly="/home/groups/harrisonlab/project_files/neonectria_ditissima/assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_filtered_renamed.fa"
  	InterProTab="gene_pred/interproscan/spades/N.ditissima/N.ditissima_interproscan.tsv"
  	SwissProtBlast="gene_pred/uniprot/N.ditissima/R0905/swissprot_v2015_09_hits.tbl"
  	SwissProtFasta="/home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta"
  	GffFile="gene_pred/augustus/N.ditissima/R0905_v2/R0905_v2_EMR_aug_preds.gff"
  	# tbl2asn options:
  	Organism="Neonectria ditissima"
  	Strain="R09/05"
  	# ncbi_tbl_corrector script options:
  	SubmissionID="AK830"
  	LabID="ArmitageEMR"
  	GeneSource='ab initio prediction:Augustus:3.1'
  	IDSource='similar to AA sequence:SwissProt:2015_09'
  	# Final submisison file name:
  	FinalName="Nd_Gomez_2015"
```

## Preparing Gff input file

Parse the Augustus Gff file.
Transcripts should be renamed as mRNA features. Exons should be added to the
Gff and unique IDs should be given to all features in the file.

```bash
  	cat $GffFile | sed 's/transcript/mRNA/g' > $OutDir/GffMRNA.gff
  	$ProgDir/generate_tbl_file/exon_generator.pl $OutDir/GffMRNA.gff > $OutDir/corrected_exons.gff
  	$ProgDir/generate_tbl_file/gff_add_id.py --inp_gff $OutDir/corrected_exons.gff --out_gff $OutDir/corrected_exons_id.gff
```

## Generating .tbl file (GAG)

The Genome Annotation Generator (GAG.py) can be used to convert gff files into
.tbl format, for use by tbl2asn.

It can also add annotations to features as provided by Annie the Annotation
extractor.

### Extracting annotations (Annie)

Interproscan and Swissprot annotations were extracted using annie, the
ANNotation Information Extractor. The output of Annie was filtered to
keep only annotations with references to ncbi approved databases.
Note - It is important that transcripts have been re-labelled as mRNA by this
point.

```bash
  	python3 $AnnieDir/annie.py -ipr $InterProTab -g $OutDir/corrected_exons_id.gff -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
  	$ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv
```

### Running GAG

Gag was run using the modified gff file as well as the annie annotation file.
Gag was noted to output database references incorrectly, so these were modified.

```bash
  	mkdir -p $OutDir/gag/round1
  	gag.py -f $Assembly -g $OutDir/corrected_exons_id.gff -a $OutDir/annie_corrected_output.csv -o $OutDir/gag/round1
  	sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl
```

## tbl2asn round 1

tbl2asn was run an initial time to collect error reports on the current
formatting of the .tbl file.
Note - all input files for tbl2asn need to be in the same directory and have the
same basename.

```bash
  	cp $Assembly $OutDir/gag/round1/genome.fsa  
  	cp $SbtFile $OutDir/gag/round1/genome.sbt
  	mkdir -p $OutDir/tbl2asn/round1
  	tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -Z discrep -j "[organism=$Organism] [strain=$Strain]"
```

## Editing .tbl file

The tbl2asn .val output files were observed and errors corrected. THis was done
with an in house script. The .val file indicated that some cds had premature
stops, so these were marked as pseudogenes ('pseudo' - SEQ_FEAT.InternalStop)
and that some genes had cds coordinates that did not match the end of the gene
if the protein was hanging off a contig ('stop' - SEQ_FEAT.NoStop).
Furthermore a number of other edits were made to bring the .tbl file in line
with ncbi guidelines. This included: Marking the source of gene
predictions and annotations ('add_inference'); Correcting locus_tags to use the
given ncbi_id ('locus_tag'); Correcting the protein and transcript_ids to
include the locus_tag and reference to submitter/lab id ('lab_id'), removal of
annotated names of genes if you don't have high confidence in their validity
(--gene_id 'remove'). If 5'-UTR and 3'-UTR were not predicted during gene
annotation then genes, mRNA and exon features need to reflect this by marking
them as incomplete ('unknown_UTR').

```bash
  	mkdir -p $OutDir/gag/edited
  	$ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --add_inference "$GeneSource" "$IDSource" --edits stop pseudo unknown_UTR --out_tbl $OutDir/gag/edited/genome.tbl
```

## Final run of tbl2asn

Following correction of the GAG .tbl file, tbl2asn was re-run to provide the
final genbank submission file.

```bash
  	cp $Assembly $OutDir/gag/edited/genome.fsa
  	cp $SbtFile $OutDir/gag/edited/genome.sbt
  	mkdir $OutDir/tbl2asn/final
  	tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -Z discrep -j "[organism=$Organism] [strain=$Strain]"
  	cp $OutDir/tbl2asn/final/genome.sqn $OutDir/tbl2asn/final/$FinalName.sqn
```

The final error report contained the following warnings. These were judged to be
legitimate concerns but biologically explainable.

67 WARNING: SEQ_FEAT.PartialProblem
 5 WARNING: SEQ_FEAT.ProteinNameEndsInBracket
211 WARNING: SEQ_FEAT.ShortExon
18 WARNING: SEQ_FEAT.SuspiciousFrame
 5 INFO:    SEQ_FEAT.PartialProblem

 Note -
 *SEQ_FEAT.partial problem. In this case, upon investigation these genes were hannging
 off the end of a contig but did not have an mRNA feature that went off of the
 end of the contig. This was occuring due to an intron being predicted hanging
 off the contig. An example on the ncbi guidelines here shows this to be
 acceptable:
 http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation#Partialcodingregionsinincompletegenomes
 *SEQ_FEAT.ProteinNameEndsInBracket. These gene names include brackets for good
 reason