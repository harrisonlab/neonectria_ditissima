

Busco has replaced CEGMA and was run to check gene space in assemblies

Previous isolates

```bash
for Assembly in $(ls repeat_masked/N.ditissima/*/*unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls repeat_masked/N.ditissima/*/*/*unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
Gene prediction

```bash
for Assembly in $(ls repeat_masked/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  for RNADir in $(ls -d qc_rna/paired/N.ditissima/Hg199); do
    Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
    echo "$Timepoint"
    FileF=$(ls $RNADir/F/*_trim.fq.gz)
    FileR=$(ls $RNADir/R/*_trim.fq.gz)
    OutDir=alignment/$Organism/$Strain/$Timepoint
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
  done
done
```
