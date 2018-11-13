# Hg199 minion assembly using miniasm

Thisis fast approach for assembling and correcting PacBio and MinION data using miniasm and racon.

## Assembly

### Assembly will with full trimming of reads:

The next 2 steps were done in the first assembly of the minion data, but repeated here.

#### Splitting reads and trimming adapters using porechop

```bash
    RawReads=raw_dna/minion/N.ditissima/Hg199/*allfiles.fastq.gz
    OutDir=qc_dna/minion/N.ditissima/Hg199_2ndround
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
```

#### Read correction using Canu

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_minion2/N.ditissima/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction.sh $TrimReads 46m $Strain $OutDir
  done
```

#### Assembly

```bash
# If reads have same names or same part splited by space, fix them using rename.sh from bbtools.
# This can only be done in blacklace11, since it needs a specific library for blasr.
ssh blacklace11.blacklace
/home/gomeza/prog/bbtools/bbmap/rename.sh in=Hg199.trimmedReads.fasta.gz out=trimmed_renamed.fasta prefix=Hg199

# Fast all-against-all overlap of raw reads
# Overlap for MinION reads (or use "-x ava-pb" for Pacbio read overlapping)
/home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > Hg199_fastq_allfiles.paf.gz

# Concatenate pieces of read sequences to generate the final sequences.
# Thus the per-base error rate is similar to the raw input reads. Make sure you correct your reads.
# Layout
/home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > Hg199_fastq_allfiles.paf.gz

# Convert gfa file to fasta file.
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > Hg199_miniasm.fa
```

## Error correction using racon:

```bash
  for Assembly in $(ls Hg199_miniasm2/Hg199_miniasm.fa); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls Hg199_miniasm2/racon_10/*10.fasta); do
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls Hg199_miniasm2/racon_10/*10.fasta); do
  Strain=Hg199
  Organism=N.ditissima
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=Hg199_miniasm2/racon_10/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Assembly correction using nanopolish

```bash
Assembly=$(ls Hg199_miniasm2/racon_10/*10.fasta)
Strain=Hg199
Organism=N.ditissima
echo "$Organism - $Strain"
ReadDir=raw_dna/nanopolish2/$Organism/$Strain
mkdir -p $ReadDir
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Event information would have been used from all of the runs, however MinKnow doesnt
# produce event-level information and therefore just the albacore data was used.
ReadsFq1=$(ls /home/groups/harrisonlab/project_files/neonectria_ditissima/raw_dna/minion/N.ditissima/Hg199/Hg199_fastq_allfiles.fastq.gz)
#ReadsFq2=$(ls raw_dna/minion/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz)
#ReadsFq1=$(ls raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/fastq_runid_298a8dbc00c3db453901232f1ad01b11fd094980_pass.fastq.gz)
#ReadsFq2=$(ls raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/fastq_runid_5832f037a56936787d17e66d1e3b8ac05572199f_pass.fastq.gz)
cat $ReadsFq1 | gunzip -cf > $ReadDir/"$Strain"_concatenated_reads.fastq
/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish/nanopolish_remove_dup_reads.py --fastq $ReadDir/"$Strain"_concatenated_reads.fastq --out $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
Fast5Dir1=$(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/raw_dna/minion/N.ditissima/Hg199/03-12-17/rebasecalled/pass/Nd_Hg199_20171203)
Fast5Dir2=$(ls -d /home/groups/harrisonlab/project_files/neonectria_ditissima/raw_dna/minion/N.ditissima/Hg199/25-10-17/rebasecalled/pass/Nd_Hg199_20171025)
nanopolish index -d $Fast5Dir1 -d $Fast5Dir2 $ReadDir/"$Strain"_concatenated_reads_filtered.fastq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $ReadDir/"$Strain"_concatenated_reads_filtered.fastq $OutDir/nanopolish
```

Split the assembly into 50Kb fragments and submit each to the cluster for nanopolish correction

```bash
screen -a

Assembly=$(ls Hg199_miniasm2/racon_10/*10.fasta)
Strain=Hg199
Organism=N.ditissima
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish2/$Organism/$Strain/"$Strain"_concatenated_reads_filtered.fastq)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)
NanoPolishDir=/home/gomeza/prog/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt
Ploidy=1
echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
```






```bash
Assembly=$(ls assembly/SMARTdenovo/N.ditissima/Hg199/racon_10/Hg199_smartdenovo_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py assembly/SMARTdenovo/$Organism/$Strain/racon_10/*/*.fa > $OutDir/"$Strain"_nanoplish.fa
# for File in $(ls assembly/SMARTdenovo/F.*/*/racon2/*/*.fa | grep ":*-"); do
#     # echo $File
#     cat $File >> $OutDir/"$Strain"_nanoplish.fa
# done
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/nanopolish/Hg199_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/N.ditissima/Hg199/nanopolish/Hg199_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
