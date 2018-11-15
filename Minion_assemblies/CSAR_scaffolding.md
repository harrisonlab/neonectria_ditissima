#CSAR scaffolding

CSAR is a web server of contig scaffolding using algebraic rearrangements
https://lu168.cs.nthu.edu.tw/CSAR-web/

I did a few tests with this tool. As reference genome, I used the R0905 pilon_5, since it was the most contiguous and complete genome assembly that I had so far. As a target I used the Hg199_minion_20k merged assembly for the same reasons.

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/Scaffold_test/*/*/*.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/Scaffold_test/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/Scaffold_test/*/*/*.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=assembly/Scaffold_test/Busco/$Organism/$Strain
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

Next will the same genomes, changing the reference to target and target to references. Hg199_minion_5k was uses as a reference.

```bash
    mv 20180806_SN104_ref.fna_scaffolding /data/scratch/gomeza/assembly/Scaffold_2/N.ditisima/Hg199/Hg199.fasta
    mv 20180806_SN104_target.fna_scaffolding /data/scratch/gomeza/assembly/Scaffold_2/N.ditisima/R0905/R0905.fasta

    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    touch tmp.csv
    for Assembly in $(ls assembly/Scaffold_2/N.ditissima/*/*.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/Scaffold_2/$Organism/$Strain
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
    done
    rm tmp.csv
```
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/Scaffold_2/*/*/*renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/Scaffold_2/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/Scaffold_2/*/*/*renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=assembly/Scaffold_2/Busco/$Organism/$Strain
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```

This merged assembly was polished using Pilon

```bash
    for Assembly in $(ls assembly/CSAR/Scaffold_v2/*/*/*renamed.fasta); do
      Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
      IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
      TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz);
      TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz);
      Iterations=5
      OutDir=assembly/CSAR/Scaffold_v2/$Organism/$Strain/polished
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
      qsub -R y $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done
```
```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/CSAR/Scaffold_v2/*/*/polished/*5.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/CSAR/Scaffold_v2/$Organism/$Strain/polished
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```
```bash
    for Assembly in $(ls assembly/CSAR/Scaffold_v2/*/*/polished/*5.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=assembly/CSAR/Scaffold_v2/Busco/$Organism/$Strain
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
```
After pilon, only two additional gene was predicted in the R0905 genome.
Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/CSAR/Scaffold_v2/N.ditissima/*/polished/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/CSAR/Scaffold_v2/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

## Repeat Masking

```bash
for Assembly in $(ls assembly/CSAR/Scaffold_v2/N.ditissima/*/filtered_contigs/*_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/Ref_Genomes/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```
The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/Ref_Genomes/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
```
```
repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4060878
repeat_masked/Ref_Genomes/N.ditissima/R0905/filtered_contigs/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5516834
```

Identify Telomere repeats:
Telomeric repeats were identified in assemblies

```bash
for Assembly in $(ls repeat_masked/Ref_Genomes/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/telomere/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
$ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
done
cat $OutDir/telomere_hits.txt | sort -nr -k5 | less
```
