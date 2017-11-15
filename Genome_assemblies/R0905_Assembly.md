## Data extraction


for N. ditissima data:
```bash
  	cd /home/groups/harrisonlab/project_files/idris
  	RawDatDir=/home/harrir/projects/pacbio_test/n_dit/
  	mkdir -p raw_dna/pacbio/N.ditissima/R0905
  	cp -r $RawDatDir/D08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/E08_1 raw_dna/pacbio/N.ditissima/R0905/.
  	cp -r $RawDatDir/F08_1 raw_dna/pacbio/N.ditissima/R0905/.
    OutDir=raw_dna/pacbio/N.ditissima/R0905/extracted
  	mkdir -p $OutDir
  	cat raw_dna/pacbio/N.ditissima/R0905/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
```


## Assembly


### Canu assembly

```bash
  	Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq)
  	GenomeSz="46m"
  	Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  	Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  	Prefix="$Strain"_canu
  	OutDir="assembly/canu/$Organism/$Strain"
  	ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  	qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

Assemblies were polished using Pilon

```bash
  	for Assembly in $(ls assembly/canu/*/R0905/*_canu.contigs.fasta); do
  	Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/NG-R0905)
    TrimF1_Read=$(ls $IlluminaDir/F/NG-R0905_qc_F.fastq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/NG-R0905_qc_R.fastq.gz);
    OutDir=assembly/canu/$Organism/$Strain/polished
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  	done
```
cat pilon.fasta | grep 'tig' | wc -l
125

### Spades Assembly


For N. ditissima

```bash
  	for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $StrainPath
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/NG-R0905)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/NG-R0905_qc_F.fastq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/NG-R0905_qc_R.fastq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    OutDir=assembly/spades_pacbio/$Organism/$Strain
    qsub $ProgDir/sub_spades_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $OutDir 15
  	done
```
cat contigs.fasta | grep 'NODE' | wc -l
641

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/spades/*/*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

Contigs shorter than 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs_min_500bp
    FilterDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs_min_500bp/contigs_min_500bp.fasta
  done
```
cat contigs_min_500bp.fasta | grep 'NODE' | wc -l
364

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs_min_500bp
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

## Merging pacbio and hybrid assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir
  done
```

```bash
  	ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  	for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  	done
```

This merged assembly was polished using Pilon

```bash
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```


#Contigs were renamed in accordance with ncbi recomendations.


```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/R0905/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Assembly stats were collected using quast

```bash
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/merged_canu_spades/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

#Canu assembly contigs were renamed in accordance with ncbi recomendations.


```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/canu/N.ditissima/R0905/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

## Renaming assemblies - temporarily
Fus2 was temporarily renamed for preliminary analysis

```bash
  cp -r assembly/canu/F.oxysporum_fsp_cepae/Fus2 assembly/canu/F.oxysporum_fsp_cepae/Fus2_pacbio_test_canu
  cp -r assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2 assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_pacbio_test_merged
  cp -r assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged_richards
```

# Repeatmasking assemblies

```bash
  R0905_pacbio_merged=$(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
  # for Assembly in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  for Assembly in $(ls $R0905_pacbio_merged); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

```bash
  R0905_pacbio_canu=$(ls assembly/canu/*/R0905/filtered_contigs/R0905_contigs_renamed.fasta)
  # for Assembly in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  for Assembly in $(ls $R0905_pacbio_canu); do
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

** % bases masked by repeatmasker: 11.57% (bases masked:5290972 bp)

** % bases masked by transposon psi: 10.29% (bases masked:4704506 bp)


Up till now we have been using just the repeatmasker/repeatmodeller fasta file when we have used softmasked fasta files. You can merge in transposonPSI masked sites using the following command:

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "$OutFile"
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done
```
repeat_masked/N.ditissima/R0905_pacbio_canu/filtered_contigs_repmask/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:5398957


# Preliminary analysis

## Checking PacBio coverage against Fus2 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.


Merged canu spades assembly

```bash
    Assembly=assembly/merged_canu_spades/N.ditissima/R0905/filtered_contigs/R0905_contigs_renamed.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```

Canu assembly

```bash
    Assembly=assembly/canu/N.ditissima/R0905/polished/pilon.fasta
    Reads=raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq
    OutDir=analysis/genome_alignment/bwa/N.ditissima/R0905/
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```
