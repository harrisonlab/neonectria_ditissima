

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
  for File in $(ls gene_pred/braker/N.*/R0905_pacbio_canu_braker_fourth/*/augustus.gff); do
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
#InBam=alignment/N.ditissima/R0905_pacbio_canu/concatenated/concatenated.bam
#ViewBam=alignment/N.ditissima/R0905_pacbio_canu/concatenated/concatenated_view.bam
#SortBam=alignment/N.ditissima/R0905_pacbio_canu/concatenated/concatenated_sorted
#samtools view -b $InBam > $ViewBam
#samtools sort $ViewBam $SortBam
#samtools index $SortBam.bam
```

[bam_header_read] EOF marker is absent. The input is probably truncated.
[main_samview] truncated file.

Cufflinks was run to compare the predicted genes to assembled transcripts:

```bash
	#for Assembly in $(ls repeat_masked/*/R0905_pacbio_canu/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
		#Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		#while [ $Jobs -gt 1 ]; do
		#sleep 10
		#printf "."
		#Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		#done
		#Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		#Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		#AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		#OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		#echo "$Organism - $Strain"
		#mkdir -p $OutDir
		#samtools merge -f $AcceptedHits \
		#alignment/$Organism/R0905_pacbio_canu/R0905/accepted_hits.bam
		#ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
		#qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	#done
```

<!--
The number of Fo47 genes was determined for comparison to number predicted by Braker (16269):
```bash
	fo47_transcripts=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf
	cat $fo47_transcripts | grep 'gene_id' | cut -f2 -d '"' | sort | uniq | wc -l
	# 18191
```
 -->
