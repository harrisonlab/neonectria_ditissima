# Circos plots

```bash

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Nd_genome=repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
OutDir=analysis/circos/N.ditissima/Hg199
mkdir -p $OutDir

# Convert the Nd genome into circos format
$ProgDir/fasta2circos.py --genome $Nd_genome --contig_prefix "" > $OutDir/Hg199_genome.txt

# Make 20kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Nd_genome --size 20 > $OutDir/Hg199_20kb_windows.gff

# Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
# NOTE - this step must not be  run on the head node as it uses a lot of RAM.

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

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
GffCAZY=gene_pred/CAZY/Ref_Genomes/N.ditissima/Hg199/Hg199_CAZY_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/Hg199_CAZY_plot.txt
GffEffP=analysis/effectorP/Ref_Genomes/N.ditissima/Hg199/N.ditissima_Hg199_EffectorP_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/Hg199_effectorP_plot.txt

circos -conf /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos/Hg199_Ref/circos.conf -outputdir ./$OutDir
```
