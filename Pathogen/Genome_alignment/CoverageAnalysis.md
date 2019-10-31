# Identify read coverage over each bp

## R0905 reference genome aligments

```bash
Reference=assembly_vAG/canu/N.ditissima/R0905/polished/repeat_masked/*_contigs_unmasked.fa
for Strain in Ag02 Ag06 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 R0905_all Hg199 Ag04 Ag05 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 OPC304 P112 NMaj; do
  for Reads in $(ls -d qc_dna/paired/N.*/$Strain)
  do
    Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
      echo "$Organism - $Strain"
      FRead=$Reads/F/*.fq.gz
      RRead=$Reads/R/*.fq.gz
      OutDir=alignment_vAG/bwa/vs_R0905_Ref
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
      qsub $ProgDir/sub_bwa.sh $Strain $Reference $FRead $RRead $OutDir
  done
done
```
```bash
  for Bam in $(ls analysis_vAG/genome_alignment/vs_R0905_Ref/*/*/*sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_R0905_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_R0905_depth.tsv > $OutDir/${Organism}_${Strain}_vs_R0905_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_R0905_depth_10kb.tsv
  done
  OutDir=analysis_vAG/genome_alignment/bwa/vs_R0905_Ref/grouped
  mkdir -p $OutDir
  cat analysis_vAG/genome_alignment/bwa_vsR0905/*/*/*_*_vs_R0905_depth_10kb.tsv > analysis_vAG/genome_alignment/bwa_vsR0905/grouped/vs_R0905_grouped_depth.tsv
```

## Hg199 reference genome aligments

```bash
Reference=assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/filtered_contigs/repeat_masked/*_contigs_unmasked.fa
for Strain in Ag02 Ag06 R37-15 R39-15 R41-15 R42-15 R45-15 R6-17-2 R6-17-3 R68-17-C2 R68-17-C3 SVK1 SVK2 R0905_all Hg199 Ag04 Ag05 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 OPC304 P112 NMaj; do
  for Reads in $(ls -d qc_dna/paired/N.*/$Strain)
  do
    Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
      echo "$Organism - $Strain"
      FRead=$Reads/F/*.fq.gz
      RRead=$Reads/R/*.fq.gz
      OutDir=alignment_vAG/bwa/vs_Hg199_Ref
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
      qsub $ProgDir/sub_bwa.sh $Strain $Reference $FRead $RRead $OutDir
  done
done
```
```bash
  for Bam in $(ls analysis_vAG/genome_alignment/bowtie/vs_Hg199_Ref/*/*/*sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_Hg199_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_Hg199_depth.tsv > $OutDir/${Organism}_${Strain}_vs_Hg199_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_Hg199_depth_10kb.tsv
  done
  OutDir=analysis_vAG/genome_alignment/bowtie/vs_Hg199_Ref/grouped
  mkdir -p $OutDir
  cat analysis_vAG/genome_alignment/bowtie/vs_Hg199_Ref/*/*/*_*_vs_Hg199_depth_10kb.tsv > analysis_vAG/genome_alignment/bowtie/vs_Hg199_Ref/grouped/vs_Hg199_grouped_depth.tsv
```
