Circos plots were generated for John's APS meeting:

```bash

ProgDir=/home/gomeza/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Nd_genome=repeat_masked/N.ditissima/R0905_pacbio_canu/filtered_contigs_repmask/R0905_contigs_unmasked.fa
OutDir=analysis/circos/N.ditissima/R0905
mkdir -p $OutDir

# Convert the Fus2 genome into circos format
$ProgDir/fasta2circos.py --genome $Nd_genome --contig_prefix "" > $OutDir/R0905_pacbio_canu_genome.txt

# Make 100kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Nd_genome --size 20 > $OutDir/R0905_pacbio_canu_20kb_windows.gff

# Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
for ReadsBam in $(ls analysis/genome_alignment/bowtie/N.*/*/vs_R0905_pacbio_new/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa_aligned.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/R0905_pacbio_canu_20kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_R0905.bed

# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
GffCAZY=gene_pred/CAZY/N.ditissima/R0905_pacbio_canu/R0905_pacbio_canu_CAZY_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/R0905_CAZY_plot.txt
GffEffP=analysis/effectorP/N.ditissima/R0905_pacbio_canu/N.ditissima_R0905_pacbio_canu_EffectorP_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/R0905_effectorP_plot.txt


circos -conf /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/circos/circos.conf -outputdir ./$OutDir
```