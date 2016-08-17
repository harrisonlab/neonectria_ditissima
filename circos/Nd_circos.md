Circos plots were generated for John's APS meeting:

```bash

  ProgDir=home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fus2_genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_unmasked.fa
  OutDir=analysis/circos/N.ditissima/R0905
  mkdir -p $OutDir

  # Convert the Fus2 genome into circos format
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > $OutDir/Fus2_genome.txt

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Fus2_genome > $OutDir/Fus2_100kb_windows.gff

  # Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
for ReadsBam in $(ls analysis/genome_alignment/bowtie/N.*/*/*/R0905_contigs_unmasked.fa_aligned_sorted.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
#bedtools coverage -abam $ReadsBam -b $OutDir/R0905_pacbio_canu_100kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_R0905.bed

    # Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_R0905.bed > $OutDir/"$Strain"_coverage_vs_R0905_scatterplot.txt
done

  # Plot location of FoL gene Blast hits as a scatterplot
  for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_test_merged/*_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > $OutDir/FoL_chr"$Chr"_genes.txt
  done

  # Plot location of Fus2 genes in pathogen-shared orthogroups as scatterplot
  GffFile=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/Fus2_genes/Fus2_path_shared_genes.gff
  $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature gene --value '1' > $OutDir/Fus2_path_shared_genes_plot.txt

  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/APS_plot/Fus2_APS_circos.conf -outputdir ./$OutDir
```
circos -conf /home/groups/harrisonlab/project_files/neonectria_ditissima/analysis/circos/N.ditissima/R0905/Fus2_APS_circos.conf -outputdir ./$OutDir

circos -conf /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/circos/circos.conf -outputdir ./$OutDir
```
