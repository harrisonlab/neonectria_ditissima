#Salmon tool for quantifying the expression of transcripts using RNA-seq data

Note that it is designed to align to predicted transcripts and not to the whole genome

Nd transcripts

```bash
for Transcriptome in $(ls gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/rna_seq/*/*); do
  FileNum=$(ls $RNADir/F/*trim.fq.gz | wc -l)
  Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    for num in $(seq 1 $FileNum); do
    printf "\n"
    FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
    FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
    echo $FileF
    echo $FileR
    Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    echo "$Timepoint"
    OutDir=alignment/salmon/$Organism/$Strain/$Prefix/$Timepoint/$Sample_Name
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
  done
done
```
Apple transcripts

```bash
for Transcriptome in $(ls apple_genome/transcripts_modified3.fasta); do
#Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
#Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/rna_seq/*/* | grep -v "mycelium"); do
  FileNum=$(ls $RNADir/F/*trim.fq.gz | wc -l)
  Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    for num in $(seq 1 $FileNum); do
    printf "\n"
    FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
    FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
    echo $FileF
    echo $FileR
    Prefix=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    echo "$Timepoint"
    OutDir=alignment/salmon/N.ditissima/Cultivar/$Prefix/$Timepoint/$Sample_Name
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
  done
done
```


Convert Salmon quasi-quanitifcations to gene counts using an awk script:

https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

```bash
mkdir -p alignment/N.ditissima/Hg199_minion/salmon/DeSeq2
for File in $(ls alignment/salmon/*/Hg199_minion/*/*/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199_minion/DeSeq2/trans2gene.txt
done
# Put files in a convenient location for DeSeq. Analysis was not performed on
# Strawberry control samples.
for File in $(ls alignment/salmon/*/*/*/*/quant.sf | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_'); do
  # cat $File | grep 'g6.t1'
  Prefix=$(echo $File | cut -f6 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done

# cp /home/groups/harrisonlab/project_files/idris/alignment/star/P.cactorum/414_v2/DeSeq/P.cactorum_RNAseq_design_parsed.txt alignment/salmon/DeSeq2/.
nano alignment/salmon/DeSeq2/P.cactorum_RNAseq_design_parsed.txt
```

/home/deakig/R3.4/bin/R
