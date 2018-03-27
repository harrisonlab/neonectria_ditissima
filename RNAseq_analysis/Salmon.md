#Salmon tool for quantifying the expression of transcripts using RNA-seq data

Note that it is designed to align to predicted transcripts and not to the whole genom

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
