Identify read coverage over each bp

R0905 reference genome aligments

```bash
  for Bam in $(ls analysis_vAG/genome_alignment/bwa_vsR0905/*/*/*sorted.bam); do
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

Hg199 CSAR genome aligments

```bash
  for Bam in $(ls analysis_vAG/genome_alignment/bowtie/vs_Hg199_CSAR/*/*/*sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_Hg199_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_Hg199_depth.tsv > $OutDir/${Organism}_${Strain}_vs_Hg199_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_Hg199_depth_10kb.tsv
  done
  OutDir=analysis_vAG/genome_alignment/bowtie/vs_Hg199_CSAR/grouped
  mkdir -p $OutDir
  cat analysis_vAG/genome_alignment/bowtie/vs_Hg199_CSAR/*/*/*_*_vs_Hg199_depth_10kb.tsv > analysis_vAG/genome_alignment/bowtie/vs_Hg199_CSAR/grouped/vs_Hg199_grouped_depth.tsv
```
