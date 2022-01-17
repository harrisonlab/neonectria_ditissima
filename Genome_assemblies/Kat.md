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

### Spades assemblies

```bash
for Strain in 118923 226-31 Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C3 SVK2 118924 227-31 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C2 SVK1; do
  for Assembly in $(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev )
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/$Strain)
    ReadsF=$(ls $IlluminaDir/F/*_trim.fq.gz)
    ReadsR=$(ls $IlluminaDir/R/*_trim.fq.gz)
    OutDir=$(dirname $Assembly)/kat
    Prefix="${Strain}"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    sbatch -p long $ProgDir/kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 200
  done
done
```
