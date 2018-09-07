# Hg199 vs R0905 circos plot

This program is used to convert fasta files into input format for circos

```bash
  OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos

  Hg199_genome=$(ls repeat_masked/*/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'Hg199')
  $ProgDir/fasta2circos.py --genome $Hg199_genome --contig_prefix "Hg_" > $OutDir/Hg199_genome.txt

  R0905_genome=$(ls repeat_masked/*/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'R0905')
  $ProgDir/fasta2circos.py --genome $R0905_genome --contig_prefix "R9_" > $OutDir/R0905_genome.txt

  cat $OutDir/Hg199_genome.txt > $OutDir/Hg199_R0905_genome.txt
  tac $OutDir/R0905_genome.txt >> $OutDir/Hg199_R0905_genome.txt

  # cat $OutDir/At_Ag_genome.txt | grep -v 'DS231' | grep -v -e 'chr23' -e 'chr24' -e 'chr25' -e 'chr26' -e 'chr27' -e 'chr28' -e 'chr29' -e 'chr30' -e 'chr31' -e 'chr32' -e 'chr33' -e 'chr34' > $OutDir/At_Ag_genome_edited.txt
```

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos
Coords=$(ls analysis/genome_alignment/mummer/N.ditissima/Hg199/Hg199_vs_R0905/Hg199_vs_R0905_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id Hg --ref_id R9 > $OutDir/Hg_vs_R9_links.txt
cat $OutDir/Hg_vs_R9_links.txt > $OutDir/Hg_vs_R9_links_edited.txt
```

A file showing contig orientations was made:
```bash
  cat $OutDir/Hg199_R0905_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/Hg_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Hg_vs_R9_links_edited.txt > $OutDir/Hg_vs_R9_contig_orientation.txt
```
<!--
The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/At_vs_Ag_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/At_vs_Ag_syntenous_contigs.txt
  cat $OutDir/At_Ag_genome.txt | grep -f $OutDir/At_vs_Ag_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
``` -->

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/Hg199_R0905_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
cat $OutDir/Hg199_R0905_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
```


Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

cat $OutDir/At_vs_As_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/At_As_genome.txt | grep 'As' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/As/, As/g'
# cat $OutDir/At_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen At contigs and remaining As contigs"
cat $OutDir/At_As_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/At/, At/g' | sed 's/As/, As/g'
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/At_vs_As/At_vs_As_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/At_vs_As_circos.png
mv $OutDir/circos.svg $OutDir/At_vs_As_circos.svg
ls $PWD/$OutDir/At_vs_As_circos.png



echo "Order of unseen At contigs and remaining As contigs"
cat $OutDir/Hg199_R0905_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g' | sed 's/R9/, R9/g'
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/Hg_vs_R9/Hg_vs_R9_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/At_vs_As_circos.png
mv $OutDir/circos.svg $OutDir/At_vs_As_circos.svg
ls $PWD/$OutDir/At_vs_As_circos.png





cat $OutDir/At_Ag_genome.txt | grep 'Ag' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Ag/, Ag/g'
cat $OutDir/At_Ag_genome.txt | grep 'At' | cut -f3 -d ' ' | tr -d '\n' | sed 's/At/, At/g' >> tmp.txt





```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/Ag_vs_At_genome_alignment/At_vs_Ag_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/At_vs_Ag_circos.png
mv $OutDir/circos.svg $OutDir/At_vs_Ag_circos.svg
```
