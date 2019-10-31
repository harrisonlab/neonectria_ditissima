# KAT: a K-mer analysis toolkit

Kat tool was used to assess level of error and duplications in the genome assemblies generated. Mapleson et al., 2016.

### Hg199 genomes

```bash
for Assembly in $(ls assembly_vAG/canu_1step/*/Hg199/pilon/*10.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev )
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  IlluminaDir=$(ls -d qc_dna/paired/*/Hg199)
  cat $IlluminaDir/F/*_trim.fq.gz > $IlluminaDir/F/F_trim_appended.fq.gz
  cat $IlluminaDir/R/*_trim.fq.gz > $IlluminaDir/R/R_trim_appended.fq.gz
  ReadsF=$(ls $IlluminaDir/F/F_trim_appended.fq.gz)
  ReadsR=$(ls $IlluminaDir/R/R_trim_appended.fq.gz)
  OutDir=$(dirname $Assembly)/kat
  Prefix="${Strain}"
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/kat
  qsub $ProgDir/sub_kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 200
done
```

after KAT jobs have finished running, then remove appended trimmed reads
```bash
rm qc_dna/paired/*/*/*/F_trim_appended.fq.gz
rm qc_dna/paired/*/*/*/R_trim_appended.fq.gz
```

### R0905 genomes

```bash
for Assembly in $(ls assembly_vAG/SMARTdenovo/N.ditissima/R0905/racon_10/pilon/R0905_pilon10_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev )
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
  cat $IlluminaDir/F/*_trim.fq.gz > $IlluminaDir/F/F_trim_appended.fq.gz
  cat $IlluminaDir/R/*_trim.fq.gz > $IlluminaDir/R/R_trim_appended.fq.gz
  ReadsF=$(ls $IlluminaDir/F/F_trim_appended.fq.gz)
  ReadsR=$(ls $IlluminaDir/R/R_trim_appended.fq.gz)
  OutDir=$(dirname $Assembly)/kat
  Prefix="${Strain}"
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/kat
  qsub $ProgDir/sub_kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 200
done
```

after KAT jobs have finished running, then remove appended trimmed reads
```bash
rm qc_dna/paired/*/*/*/F_trim_appended.fq.gz
rm qc_dna/paired/*/*/*/R_trim_appended.fq.gz
```
