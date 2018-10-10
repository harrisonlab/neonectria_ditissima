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
```
Telomere locations on contigs:

```bash
cat analysis/telomere/N.ditissima/Hg199/telomere_hits_circos.txt | sed 's/contig/Hg_contig/g' | sort -k3 -n -t'_' > $OutDir/Hg_vs_R9_telomere_hits.txt
cat analysis/telomere/N.ditissima/R0905/telomere_hits_circos.txt  | sed 's/contig/R9_contig/g' | sort -k3 -n -t'_' >> $OutDir/Hg_vs_R9_telomere_hits.txt
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

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/Hg199_R0905_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
cat $OutDir/Hg199_R0905_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
```
```bash
echo "Order of unseen Hg contigs and remaining R9 contigs"
cat $OutDir/Hg199_R0905_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R9/, R9/g' | sed 's/Hg/, Hg/g'
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos
circos -conf $ProgDir/Hg_vs_R9_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Hg_vs_R9_circos.png
mv $OutDir/circos.svg $OutDir/Hg_vs_R9_circos.svg
ls $PWD/$OutDir/Hg_vs_R9_circos.png
```
#2nd Method

orthologs were used to create a file with coordonates
```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/OrthoFinder/formatted/Results_Sep24/Orthogroups.txt  --name1 199R --gff1 gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3  > $OutDir/Hg_vs_R9_links_gene.txt --name2 R09R --gff2 gene_pred/codingquary/Ref_Genomes/N.ditissima/R0905/final/final_genes_appended_renamed.gff3
cat $OutDir/Hg_vs_R9_links_gene.txt > $OutDir/Hg_vs_R9_links_gene_edited.txt
```
A file showing contig orientations was made:
```bash
  cat $OutDir/Hg199_R0905_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/Hg_contig_order2.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Hg_vs_R9_links_gene_edited.txt > $OutDir/Hg_vs_R9_contig_orientation2.txt
```

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/Hg_vs_R9_contig_orientation2.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Hg_vs_R9_contig_orientation2.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/Hg199_R0905_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
cat $OutDir/Hg199_R0905_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
```
```bash
echo "Order of unseen Hg contigs and remaining R9 contigs"
cat $OutDir/Hg199_R0905_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R9/, R9/g' | sed 's/Hg/, Hg/g'
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos
circos -conf $ProgDir/Hg_vs_R9_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Hg_vs_R9_circos2.png
mv $OutDir/circos.svg $OutDir/Hg_vs_R9_circos2.svg
ls $PWD/$OutDir/Hg_vs_R9_circos2.png
```
#3rd Method

This program is used to convert fasta files into input format for circos

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos_NoCSAR
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos

Hg_genome=$(ls assembly/NoCSAR/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_renamed.fasta)
$ProgDir/fasta2circos.py --genome $Hg_genome --contig_prefix "Hg_" > $OutDir/Hg_genome.txt

R9_genome=$(ls assembly/NoCSAR/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
$ProgDir/fasta2circos.py --genome $R9_genome --contig_prefix "R9_" > $OutDir/R9_genome.txt

cat $OutDir/Hg_genome.txt > $OutDir/Hg_R9_genome.txt
tac $OutDir/R9_genome.txt >> $OutDir/Hg_R9_genome.txt
```

Identify Telomere repeats:
Telomeric repeats were identified in assemblies

```bash
for Assembly in $(ls assembly/NoCSAR/N.ditissima/*/filtered_contigs/*_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/telomere/NoCSAR/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
$ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
done
```

Telomere locations on contigs:

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos_NoCSAR
cat analysis/telomere/NoCSAR/N.ditissima/Hg199/telomere_hits_circos.txt | sed 's/contig/Hg_contig/g' | sort -k3 -n -t'_' > $OutDir/Hg_vs_R9_telomere_hits.txt
cat analysis/telomere/NoCSAR/N.ditissima/R0905/telomere_hits_circos.txt  | sed 's/contig/R9_contig/g' | sort -k3 -n -t'_' >> $OutDir/Hg_vs_R9_telomere_hits.txt
```

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_circos_NoCSAR
Coords=$(ls analysis/genome_alignment/mummer/N.ditissima/Hg199/Hg199_vs_R0905_50contigs/Hg199_vs_R0905_50contigs_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id Hg --ref_id R9 > $OutDir/Hg_vs_R9_links.txt
cat $OutDir/Hg_vs_R9_links.txt > $OutDir/Hg_vs_R9_links_edited.txt
```

A file showing contig orientations was made:
```bash
  cat $OutDir/Hg_R9_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/Hg_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Hg_vs_R9_links_edited.txt > $OutDir/Hg_vs_R9_contig_orientation.txt
```

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/Hg_R9_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
cat $OutDir/Hg_R9_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
```
```bash
echo "Order of unseen Hg contigs and remaining R9 contigs"
cat $OutDir/Hg_R9_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R9/, R9/g' | sed 's/Hg/, Hg/g'
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos/Hg_vs_R9_NoCSAR
circos -conf $ProgDir/Hg_vs_R9_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Hg_vs_R9_circos.png
mv $OutDir/circos.svg $OutDir/Hg_vs_R9_circos.svg
ls $PWD/$OutDir/Hg_vs_R9_circos.png
```

#4rd Method

```bash
  OutDir=analysis/circos/Hg_vs_R9_genome_alignment_hard
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos

  Hg199_genome=$(ls repeat_masked/Ref_Genomes/*/*/filtered_contigs/*_contigs_hardmasked.fa | grep 'Hg199')
  $ProgDir/fasta2circos.py --genome $Hg199_genome --contig_prefix "Hg_" > $OutDir/Hg199_genome.txt

  R0905_genome=$(ls repeat_masked/Ref_Genomes/*/*/filtered_contigs/*_contigs_hardmasked.fa | grep 'R0905')
  $ProgDir/fasta2circos.py --genome $R0905_genome --contig_prefix "R9_" > $OutDir/R0905_genome.txt

  cat $OutDir/Hg199_genome.txt > $OutDir/Hg199_R0905_genome.txt
  tac $OutDir/R0905_genome.txt >> $OutDir/Hg199_R0905_genome.txt
```

Identify Telomere repeats:
Telomeric repeats were identified in assemblies

```bash
for Assembly in $(ls repeat_masked/Ref_Genomes/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/telomere/CSAR_Hard/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
$ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
done
```

Telomere locations on contigs:

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_hard
cat analysis/telomere/CSAR_Hard/N.ditissima/Hg199/telomere_hits_circos.txt | sed 's/contig/Hg_contig/g' | sort -k3 -n -t'_' > $OutDir/Hg_vs_R9_telomere_hits.txt
cat analysis/telomere/CSAR_Hard/N.ditissima/R0905/telomere_hits_circos.txt  | sed 's/contig/R9_contig/g' | sort -k3 -n -t'_' >> $OutDir/Hg_vs_R9_telomere_hits.txt
```

```bash
OutDir=analysis/circos/Hg_vs_R9_genome_alignment_hard
Coords=$(ls analysis/genome_alignment/mummer/N.ditissima/Hg199/Hg199_vs_R0905_hard/Hg199_vs_R0905_hard_coords.tsv)
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

Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Hg_vs_R9_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/Hg199_R0905_genome.txt | grep 'Hg' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Hg/, Hg/g'
cat $OutDir/Hg199_R0905_genome.txt | grep 'R9' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R9/, R9/g' >> tmp.txt
echo "Order of unseen Hg contigs and remaining R9 contigs"
cat $OutDir/Hg199_R0905_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R9/, R9/g' | sed 's/Hg/, Hg/g'
```
```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Pathogen/Circos/Hg_vs_R9_Hard
circos -conf $ProgDir/Hg_vs_R9_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Hg_vs_R9_circos.png
mv $OutDir/circos.svg $OutDir/Hg_vs_R9_circos.svg
ls $OutDir/Hg_vs_R9_circos.png
```
