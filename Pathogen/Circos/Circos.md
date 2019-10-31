# Representation of the genome of Nd using Circos plots

### MiSeq reads alignment to the Reference genome using bowtie

```bash
  Reference=$(ls analysis_vAG/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/N.ditissima/Ag0*); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis_vAG/genome_alignment/bowtie/$Organism/$Strain/vs_R0905_Ref
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
```

### Circos plot of R0905

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos
Nd_genome=assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
OutDir=analysis_vAG/circos/N.ditissima/R0905
mkdir -p $OutDir

# Convert the Nd genome into circos format
$ProgDir/fasta2circos.py --genome $Nd_genome --contig_prefix "" > $OutDir/R0905_genome.txt

# Make 20kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Nd_genome --size 20 > $OutDir/R0905_20kb_windows.gff

# Convert Nd isolates MiSeq reads aligning in 100kb windows into coverage stats
# NOTE - this step must not be run on the head node as it uses a lot of RAM.
# R6-17-2 (moderate pathogenic UK isolate)
for ReadsBam in $(ls analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/R6-17-2/R6-17-2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f2 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -a $OutDir/R0905_20kb_windows.gff -b $ReadsBam > $AlignDir/"$Strain"_coverage_vs_R0905.bed

# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

# R45-15 (high pathogenic UK isolate)
for ReadsBam in $(ls analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/R45-15/R45-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f2 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -a $OutDir/R0905_20kb_windows.gff -b $ReadsBam > $AlignDir/"$Strain"_coverage_vs_R0905.bed

# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

# R37-15 (low pathogenic Belgium isolate)
for ReadsBam in $(ls analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/R37-15/R37-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f2 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -a $OutDir/R0905_20kb_windows.gff -b $ReadsBam > $AlignDir/"$Strain"_coverage_vs_R0905.bed

# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

# Hg199 (high pathogenic UK isolate)
for ReadsBam in $(ls analysis_vAG/genome_alignment/bowtie/vs_R0905/N.ditissima/R0905/Hg199/Hg199_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f2 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -a $OutDir/R0905_20kb_windows.gff -b $ReadsBam > $AlignDir/"$Strain"_coverage_vs_R0905.bed

# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

# Plot location of secreted CAZymes and effectorP genes as a scatterplot
GffCAZY=gene_pred_vAG/CAZY_NewDatabase/Ref_Genomes/N.ditissima/R0905/R0905_CAZY_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/R0905_CAZY_plot.txt

GffEffP=analysis_vAG/effectorP/Ref_Genomes/N.ditissima/R0905/N.ditissima_R0905_EffectorP_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/R0905_effectorP_plot.txt

circos -conf /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/circos/circos.conf -outputdir ./$OutDir
```

### Circos plot of Hg199

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos
Nd_genome=repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
OutDir=analysis_vAG/circos/N.ditissima/Hg199
mkdir -p $OutDir

# Convert the Nd genome into circos format
$ProgDir/fasta2circos.py --genome $Nd_genome --contig_prefix "" > $OutDir/Hg199_genome.txt

# Make 20kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Nd_genome --size 20 > $OutDir/Hg199_20kb_windows.gff

# Convert Nd isolates MiSeq reads aligning in 100kb windows into coverage stats
# NOTE - this step must not be run on the head node as it uses a lot of RAM.

qlogin -pe smp 16 -l virtual_free=1G

#for Strain in Ag04 Ag05 BGV344 OPC304 R68-17-C2 SVK1 SVK2 NMaj; do
for Strain in R6-17-2 R37-15 R6-17-3; do
  for ReadsBam in $(ls ../../../home/groups/harrisonlab/project_files/neonectria_ditissima/analysis/genome_alignment/bowtie/N.*/$Strain/"$Strain"_unmasked.fa_aligned.bam); do
  Organism=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $ReadsBam | rev | cut -f2 -d '/' | rev)
  AlignDir=$(dirname $ReadsBam)
  echo "$Organism - $Strain"
  bedtools coverage -a $OutDir/Hg199_20kb_windows.gff -b $ReadsBam > $AlignDir/"$Strain"_coverage_vs_Hg199.bed

  # Convert coverage bed files into circos format
  $ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_Hg199.bed > $OutDir/"$Strain"_coverage_vs_Hg199_scatterplot.txt
  done
done

# Plot location of secreted CAZymes and effectorP genes as a scatterplot
GffCAZY=gene_pred/CAZY/Ref_Genomes/N.ditissima/Hg199/Hg199_CAZY_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/Hg199_CAZY_plot.txt

GffEffP=analysis/effectorP/Ref_Genomes/N.ditissima/Hg199/N.ditissima_Hg199_EffectorP_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/Hg199_effectorP_plot.txt

circos -conf /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos/Hg199_Ref/circos.conf -outputdir ./$OutDir
```
