# MUMmer was run to align assemblies against the reference genome.

MUMmer was run to align assemblies against the reference genome.

### Hg199 against R0905 genome.

```bash
Reference=$(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/filtered_contigs/R0905_contigs_unmasked.fa)
for Query in $(ls repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_R0905
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

### R0905 against Hg199 genome.

```bash
Reference=$(ls repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa)
for Query in $(ls repeat_masked/Ref_Genomes/N.ditissima/R0905/filtered_contigs/R0905_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Hg199
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

### No CSAR genomes. Hg199 against R0905 genome (50 contigs).

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/canu_pacbio/N.ditissima/R0905/Original_v3/polished/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/NoCSAR/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv

  touch tmp.csv
  for Assembly in $(ls assembly/merged_SMART_spades/Hg199_minion_5k/pilon/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/NoCSAR/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

```bash
  Reference=$(ls assembly/NoCSAR/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
    for Query in $(ls assembly/NoCSAR/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_renamed.fasta); do
    Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    Prefix="$Strain"_vs_R0905_50contigs
    OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
    qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
  done
```
