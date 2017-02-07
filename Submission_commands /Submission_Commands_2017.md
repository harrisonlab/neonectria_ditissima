# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the F. oxysporum f. sp.
narcissi genome for submisson are shown below.

# Preliminary submission

A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov

A preliminary submission was made for the .fasta assembly to check if
any contigs needed to be split. This step was performed early in the annotation
process (prior to gene prediction) to ensure that annotation did not have to
be repeated at the end of the project.


The following note was provided in the WGS submission page on NCBI in the box
labeled "Private comments to NCBI staff":

I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'

# Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  ls genome_submission/SRA_metadata_acc.txt genome_submission/SRA_metadata_acc.xlsx
```

As these files included a file > 500 Mb, a presubmission folder was requested.
This aids submission of large data files. This file was created on the ftp server
at ftp-private.ncbi.nlm.nih.gov, with a private folder named
uploads/andrew.armitage@emr.ac.uk_6L2oakBI. Ncbi provided a username a password.
Files were uploaded into a folder created within my preload folder using ftp.

```bash
# Bioproject="PRJNA369735"
SubFolder="Ndit_PRJNA369735"
mkdir $SubFolder
for Read in $(ls raw_dna/paired/*/Hg199/*/*.fastq.gz); do
echo $Read;
cp $Read $SubFolder/.
done
cd $SubFolder


ftp ftp-private.ncbi.nlm.nih.gov

# Enter user name and password. Copy and paste allowed. It has to be quick.

cd uploads/antonio.gomez@emr.ac.uk_EpDamnXg
mkdir Ndit_PRJNA369735
cd Ndit_PRJNA369735
# put Nd_PRJNA369735
prompt
mput *
bye
cd ../
rm -r $SubFolder
```
# The process of submissions ends here. You might need to provide more data, for example coverage



For WGS
# Calculate coverage using the count_nucl.pl script

## For one library

```bash
for DataDir in $(ls -d raw_dna/paired/*/Hg199)
do
F_Read=$(ls $DataDir/F/*.gz)
R_Read=$(ls $DataDir/R/*.gz)
Strain=$(echo $DataDir | rev | cut -f1 -d '/' | rev)
Organism=$(echo $DataDir | rev | cut -f2 -d '/' | rev)
WorkDir=tmp_dir/$Strain
mkdir -p $WorkDir
cp -r $F_Read $WorkDir
cp -r $R_Read $WorkDir
cd $WorkDir
Read1=*F*
Read2=*R*
gunzip $Read1
gunzip $Read2
Sub1=*F*.fq
Sub2=*R*.fq
echo "$Organism - $Strain"
count_nucl.pl -i $Sub1 -i $Sub2 -g 35
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/clocks
done
```

WGS submission: SUB2105171 Vd12161
I would like to ask to NCBI to run the assembly thorough the contamination screen as this is not the final submission.
contigs_min_500bp.fasta


V.dahliae - 51
The estimated genome size is: 35000000 bp


The input file is: wilt_51_F_appended_trim.fq

Results for: wilt_51_F_appended_trim.fq
 Within this file of 1116178220 bp there were 6515637 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_51_R_appended_trim.fq

Results for: wilt_51_R_appended_trim.fq
 Within this file of 1111891504 bp there were 6515637 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 2228069724 nucleotides in this file.

 This equates to an estimated genome coverage of 63.66 .


V.dahliae - 53
The estimated genome size is: 35000000 bp


The input file is: wilt_53_F_appended_trim.fq

Results for: wilt_53_F_appended_trim.fq
 Within this file of 790141507 bp there were 3597693 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_53_R_appended_trim.fq

Results for: wilt_53_R_appended_trim.fq
 Within this file of 757410064 bp there were 3597693 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 1547551571 nucleotides in this file.

 This equates to an estimated genome coverage of 44.22 .


V.dahliae - 58
The estimated genome size is: 35000000 bp


The input file is: wilt_58_F_appended_trim.fq

Results for: wilt_58_F_appended_trim.fq
 Within this file of 1404169996 bp there were 6746472 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_58_R_appended_trim.fq

Results for: wilt_58_R_appended_trim.fq
 Within this file of 1385769457 bp there were 6746472 fastq sequences
 of these 0 lines were empty.

Total results:

 There are a total of 2789939453 nucleotides in this file.

 This equates to an estimated genome coverage of 79.71 .


V.dahliae - 61
The estimated genome size is: 35000000 bp


The input file is: wilt_61_F_appended_trim.fq

Results for: wilt_61_F_appended_trim.fq
 Within this file of 1087035735 bp there were 4990922 fastq sequences
 of these 0 lines were empty.


The input file is: wilt_61_R_appended_trim.fq

Results for: wilt_61_R_appended_trim.fq
 Within this file of 1072541654 bp there were 4990922 fastq sequences
 of these 0 lines were empty.
