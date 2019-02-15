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

### Hg199 miniasm genome against R0905 canu genome

```bash
Reference=$(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/R0905_pilon10_renamed.fasta)
for Query in $(ls assembly_vAG/miniasm/N.ditissima/Hg199/racon_10/pilon/Hg199_pilon10_renamed.fasta); do
Strain=$(echo $Query | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_R0905_vAG
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
