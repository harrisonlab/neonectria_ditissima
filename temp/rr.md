

Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
for File in $(ls repeat_masked/spades/N.*/R0905_v2/filtered_contigs_repmask/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done
```

repeat_masked/spades/N.ditissima/R0905_v2/filtered_contigs_repmask/R0905_v2_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5504315

q2. Alignment of raw reads vs the Fus2 genome

Sequence data for isolates with a data from a single sequencing run was aligned against the Fus2 genome

```bash
Reference=$(ls repeat_masked/spades/N.*/R0905_v2/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
for StrainPath in $(ls -d qc_dna/paired/N.ditissima/R0905); do
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*.fq.gz)
R_Read=$(ls $StrainPath/R/*.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_appended3
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
  ```

  ```bash
Reference=$(ls repeat_masked/spades/N.*/Hg199/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
for StrainPath in $(ls -d qc_dna/paired/N.ditissima/Hg199); do
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*.fq.gz)
R_Read=$(ls $StrainPath/R/*.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Hg199_appended3
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
    ```

    ```bash
Reference=$(ls repeat_masked/N.*/R0905_merged_assembly/*/*_contigs_softmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/N.ditissima/*); do
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*.fq.gz)
R_Read=$(ls $StrainPath/R/*.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_merged_assembly_softmasked
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
      ```
