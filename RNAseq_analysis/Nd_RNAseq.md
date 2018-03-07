
##RNA-Seq data was downloaded from novogenes servers with the following commands

```bash
mkdir -p /home/scratch/gomeza/rna_seq/
tar -xvf /home/scratch/gomeza/rna_seq/20171211/C101HW17030405_20180102_5_Yvad6z.tar
```

##Reorganise raw data

```bash
mkdir -p N.ditissima/Hg199/mycellium/F/
mkdir -p N.ditissima/Hg199/mycellium/R/
mkdir -p N.ditissima/GD/t0/F/
mkdir -p N.ditissima/GD/t1/F/
mkdir -p N.ditissima/GD/t2/F/
mkdir -p N.ditissima/GD/t0/R/
mkdir -p N.ditissima/GD/t1/R/
mkdir -p N.ditissima/GD/t2/R/
mkdir -p N.ditissima/M9/t0/F/
mkdir -p N.ditissima/M9/t1/F/
mkdir -p N.ditissima/M9/t2/F/
mkdir -p N.ditissima/M9/t0/R/
mkdir -p N.ditissima/M9/t1/R/
mkdir -p N.ditissima/M9/t2/R/

mv 20171211/C101HW17030405/raw_data/GD_C1_3_1.fq.gz N.ditissima/GD/t0/F/
mv 20171211/C101HW17030405/raw_data/GD_C1_3_2.fq.gz N.ditissima/GD/t0/R/
mv 20171211/C101HW17030405/raw_data/GD_6A3_1.fq.gz N.ditissima/GD/t1/F/
mv 20171211/C101HW17030405/raw_data/GD_6A3_2.fq.gz N.ditissima/GD/t1/R/
mv 20171211/C101HW17030405/raw_data/M9_C1_3_2.fq.gz N.ditissima/M9/t0/R/
mv 20171211/C101HW17030405/raw_data/M9_C1_3_1.fq.gz N.ditissima/M9/t0/F/
mv 20171211/C101HW17030405/raw_data/M9_6A3_1.fq.gz N.ditissima/M9/t1/F/
mv 20171211/C101HW17030405/raw_data/M9_6A3_2.fq.gz N.ditissima/M9/t1/R/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A4_1.fq.gz N.ditissima/GD/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A4_2.fq.gz N.ditissima/GD/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A3_1.fq.gz N.ditissima/GD/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A3_2.fq.gz N.ditissima/GD/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A5_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_4A5_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A2_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_6A2_1.fq.gz N.ditissima/GD/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_6A2_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_5A2_2.fq.gz N.ditissima/GD/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C2_3_1.fq.gz N.ditissima/GD/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C2_3_2.fq.gz N.ditissima/GD/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C3_3_2.fq.gz N.ditissima/GD/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/GD_C3_3_1.fq.gz N.ditissima/GD/t0/F/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C2_3_1.fq.gz N.ditissima/M9/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C3_3_1.fq.gz N.ditissima/M9/t0/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C3_3_2.fq.gz N.ditissima/M9/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_C2_3_2.fq.gz N.ditissima/M9/t0/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A3_1.fq.gz N.ditissima/M9/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A3_2.fq.gz N.ditissima/M9/t1/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_2A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_6A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_6A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A2_2.fq.gz N.ditissima/M9/t2/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A2_1.fq.gz N.ditissima/M9/t2/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A4_1.fq.gz N.ditissima/M9/t1/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/M9_5A4_2.fq.gz N.ditissima/M9/t1/R/

cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_1_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_2_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_3_1.fq.gz N.ditissima/Hg199/mycelium/F/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_3_2.fq.gz N.ditissima/Hg199/mycelium/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_2_2.fq.gz N.ditissima/Hg199/mycelium/R/
cp /data/seq_data/external/20180118_novogene_gomeza/C101HW17030405/data_release/C101HW17030405/raw_data/Hg199_1_2.fq.gz N.ditissima/Hg199/mycelium/R/
```

##Perform qc on RNA-Seq timecourse and mycelium data

```bash
for FilePath in $(ls -d rna_seq/N.ditissima/*/*)
do
    echo $FilePath
    FileNum=$(ls $FilePath/F/*.gz | wc -l)
    for num in $(seq 1 $FileNum)
    do
        FileF=$(ls $FilePath/F/*.gz | head -n $num | tail -n1)
        FileR=$(ls $FilePath/R/*.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Jobs=$(qstat -u "*" | grep 'rna_qc' | grep 'qw' | wc -l)
        while [ $Jobs -gt 16 ]
        do
            sleep 5m
            printf "."
            Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
        done        
        printf "\n"
        IlluminaAdapters=/home/gomeza/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
        ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/rna_qc
        qsub -h $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
        JobID=$(qstat | grep 'rna' | tail -n 1 | cut -d ' ' -f1)
        Queue_Status=$(qstat | grep 'rna' | grep 'hqw' | wc -l)
        while (($Queue_Status > 0))
        do
            Queue_Status=$(qstat | grep 'rna' | grep 'hqw' | wc -l)
            load02=$(qstat -u "*" | grep 'blacklace02'| grep 'rna' | wc -l)
            load05=$(qstat -u "*" | grep 'blacklace05'| grep 'rna' | wc -l)
            load06=$(qstat -u "*" | grep 'blacklace06'| grep 'rna' | wc -l)
            load10=$(qstat -u "*" | grep 'blacklace10'| grep 'rna' | wc -l)
            if (($load02 < 3))
            then
                qalter $JobID -l h=blacklace02.blacklace
                sleep 5s
                qalter $JobID -h U
                sleep 5s
                echo "Submitted to node 2"
            elif (($load05 < 3))
            then
                qalter $JobID -l h=blacklace05.blacklace
                sleep 5s
                qalter $JobID -h U
                sleep 5s
                echo "Submitted to node 5"
            elif (($load06 < 3))
            then
                qalter $JobID -l h=blacklace06.blacklace
                sleep 5s
                qalter $JobID -h U
                sleep 5s
                echo "Submitted to node 6"
            elif (($load10 < 3))
            then
                qalter $JobID -l h=blacklace10.blacklace
                sleep 5s
                qalter $JobID -h U
                sleep 5s
                echo "Submitted to node 10"
            else
                echo "all nodes full, waiting ten minutes"
                sleep 10m
            fi
        done    
    done
done
```

###Visualise data quality using fastqc

Only submit three jobs at a time, copying 30 files is too much!

```bash
for RawData in $(ls qc_rna/N.ditissima/Hg199/*/*/*_trim.fq.gz | grep 'Hg199_1')
do
    echo $RawData
    Jobs=$(qstat -u "*" | grep 'run_fastqc' | wc -l)
    while [ $Jobs -gt 3 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat -u "*" | grep 'run_fastqc' | wc -l)
    done
    ProgDir=/home/adamst/git_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/run_fastqc.sh $RawData
done
```

```
All seem okay to me
```

#Align mycelium reads to Hg199 minion assembly with STAR

```bash
for Assembly in $(ls repeat_masked/N.ditissima/Ref_Genomes/Hg199_minion/*/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/RNAseq/N.ditissima/Hg199/mycelium/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
        echo $FileF
        echo $FileR
        Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
        echo "$Timepoint"
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq/
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done
```
2nd run
```bash
for Assembly in $(ls Hg199_genome/*_contigs_unmasked.fa)
do
    Strain=Hg199_minion
    Organism=N.ditissima
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/N.ditissima/Hg199/mycelium/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
        echo $FileF
        echo $FileR
        Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
        echo "$Timepoint"
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq/
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done
```

#Align all timepoints to apple genome.

```bash
for FileF in $(ls /data/scratch/gomeza/qc_rna/N.ditissima/M9/t1/F/M9_6*_trim.fq.gz | grep -v "mycelium")
do
    Strain=$(echo $FileF | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $FileF | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    done
    printf "\n"
    FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1/_2/g')
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=/data/scratch/gomeza/alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    Assembly=/data/scratch/gomeza/apple_genome/GDDH13_1-1_formatted.fasta
    GFF=/home/groups/harrisonlab/project_files/neonectria_ditissima/apple_genome/gene_models_20170612.gff3
    qsub $ProgDir/sub_star_sensitive.sh $Assembly $FileF $FileR $OutDir $GFF
    done
done
```
2nd run
```bash
for FileF in $(ls /data/scratch/gomeza/qc_rna/N.ditissima/M9/t*/F/*_trim.fq.gz | grep -v "mycelium")
do
    Strain=$(echo $FileF | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $FileF | rev | cut -f5 -d '/' | rev)
    echo "$Organism - $Strain"
    Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
    done
    printf "\n"
    FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1/_2/g')
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=/data/scratch/gomeza/alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    mkdir -p $OutDir
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    Assembly=/data/scratch/gomeza/apple_genome/GDDH13_1-1_formatted.fasta
    GFF=/data/scratch/gomeza/apple_genome/gene_models_20170612.gff3
    qsub $ProgDir/sub_star_sensitive.sh $Assembly $FileF $FileR $OutDir $GFF
    done
done
```

##Gzip output files to save space on the disk and allow star to run correctly downstream. ONLY RUN THIS ONCE

```bash
for AlignDir in $(ls -d alignment/star/N.ditissima/*/*/*)
do
    cat $AlignDir/star_aligmentUnmapped.out.mate1 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate1.fq.gz
    cat $AlignDir/star_aligmentUnmapped.out.mate2 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate2.fq.gz
done
```

##Align compressed files of unmapped reads from aligning to Golden Delicious genome

This star script had the following options added to the sub_star.sh script in the ProgDir specified in the below commands:
--winAnchorMultimapNmax 200
--seedSearchStartLmax 30

```bash
for Assembly in $(ls /data/scratch/gomeza/Hg199_genome/*_contigs_unmasked.fa)
do
  #Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  #Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  #echo "$Organism - $Strain"
  for AlignDir in $(ls -d /data/scratch/gomeza/alignment/star/N.ditissima/Hg199_minion/M9/t1/M9_2*)
  do
    Strain=$(echo $AlignDir | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $AlignDir | rev | cut -f5 -d '/' | rev)
      echo "$Organism - $Strain"
      printf "\n"
      File1=$AlignDir/star_aligmentUnmapped.out.mate1.fq.gz
      File2=$AlignDir/star_aligmentUnmapped.out.mate2.fq.gz
      echo $File1
      echo $File2
      Timepoint=$(echo $AlignDir | rev | cut -d '/' -f2 | rev)
      echo "$Timepoint"
      Sample_Name=$(echo $AlignDir | rev | cut -d '/' -f1 | rev)
      Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
      while [ $Jobs -gt 1 ]
      do
          sleep 1m
          printf "."
          Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
      done
      OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/sub_star_TA.sh $Assembly $File1 $File2 $OutDir
  done
done
```

```bash
for Assembly in $(ls Hg199_genome/*_contigs_unmasked.fa)
do
  #Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  #Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
  #echo "$Organism - $Strain"
  for AlignDir in $(ls -d alignment/star/N.ditissima/Hg199_minion/M9/t1/*)
  do
    Strain=$(echo $AlignDir | rev | cut -f4 -d '/' | rev)
    Organism=$(echo $AlignDir | rev | cut -f5 -d '/' | rev)
      echo "$Organism - $Strain"
      printf "\n"
      File1=$AlignDir/star_aligmentUnmapped.out.mate1.fq.gz
      File2=$AlignDir/star_aligmentUnmapped.out.mate2.fq.gz
      echo $File1
      echo $File2
      Timepoint=$(echo $AlignDir | rev | cut -d '/' -f2 | rev)
      echo "$Timepoint"
      Sample_Name=$(echo $AlignDir | rev | cut -d '/' -f1 | rev)
      Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
      while [ $Jobs -gt 1 ]
      do
          sleep 1m
          printf "."
          Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
      done
      OutDir=alignment/star/star2/$Organism/$Strain/$Timepoint/$Sample_Name
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/sub_star.sh $Assembly $File1 $File2 $OutDir
  done
done
```

#Quantification of gene models

```bash
    for BamFile in $(ls -d alignment/star/N.ditissima/Hg199_minion/*/GD*/star_aligmentAligned.sortedByCoord.out.bam)
    do
        Gff=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3
        OutDir=$(dirname $BamFile)
        Prefix=$(echo $BamFile | rev | cut -f2 -d '/' | rev)
        Jobs=$(qstat | grep 'sub_fea' | wc -l)
        while [ $Jobs -gt 5 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_fea' | wc -l)
        done
        printf "\n"
        echo $Prefix
        ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_featureCounts.sh $BamFile $Gff $OutDir $Prefix
    done
```
```bash
      for BamFile in $(ls -d alignment/star/N.ditissima/Hg199_minion/GD/t2/*/star_aligmentAligned.sortedByCoord.out.bam)
      do
          Gff=apple_genome/gene_models_20170612.gff3
          Time=$(echo $BamFile | rev | cut -f3 -d '/' | rev)
          Prefix=$(echo $BamFile | rev | cut -f2 -d '/' | rev)
          OutDir=alignment/star/N.ditissima/Hg199_minion/cultivars/$Time/$Prefix
          Jobs=$(qstat | grep 'sub_fea' | wc -l)
          while [ $Jobs -gt 5 ]
          do
              sleep 1m
              printf "."
              Jobs=$(qstat | grep 'sub_fea' | wc -l)
          done
          printf "\n"
          echo $Prefix
          ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
          qsub $ProgDir/sub_featureCounts.sh $BamFile $Gff $OutDir $Prefix
      done
```

##A file was created with columns referring to experimental treatments

```bash
OutDir=alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq
mkdir -p $OutDir
printf "Sample.name\tCultivar\tTimepoint\n" > $OutDir/N.dit_Hg199_RNAseq_design.txt
# for File in $(ls alignment/star/P.cactorum/10300/Sample_*/Sample_*_featurecounts.txt); do
# Sample=$(echo $File | rev | cut -f2 -d '/' | rev)
# i=$(echo $Sample | sed 's/Sample_//g')
for i in $(seq 1 24)
do
  if [ $i == '10' ] || [ $i == '11' ] || [ $i == '12' ]
  then
      Timepoint='Control'
      Cultivar='mycelium'
  elif [ $i == '7' ] || [ $i == '8' ] || [ $i == '9' ]
  then
      Timepoint='t0'
      Cultivar='GD'
    elif [ $i == '19' ] || [ $i == '20' ] || [ $i == '21' ]
    then
        Timepoint='t0'
        Cultivar='M9'
  elif [ $i == '1' ] || [ $i == '4' ] || [ $i == '6' ]
  then
      Timepoint='t1'
      Cultivar='GD'
  elif [ $i == '14' ] || [ $i == '16' ] || [ $i == '18' ]
  then
      Timepoint='t1'
      Cultivar='M9'
  elif [ $i == '2' ] || [ $i == '3' ] || [ $i == '5' ]
  then
        Timepoint='t2'
        Cultivar='GD'
  elif [ $i == '13' ] || [ $i == '15' ] || [ $i == '17' ]
  then
        Timepoint='t2'
        Cultivar='M9'   
  fi
  if [ $i == '10' ]
  then
      printf "Hg199_1\t$Cultivar\t$Timepoint\n"
  elif [ $i == '11' ]
  then
      printf "Hg199_2\t$Cultivar\t$Timepoint\n"
  elif [ $i == '12' ]
  then
      printf "Hg199_3\t$Cultivar\t$Timepoint\n"
  elif [ $i == '7' ]
  then
      printf "GD_C1_3\t$Cultivar\t$Timepoint\n"
  elif [ $i == '8' ]
  then
      printf "GD_C2_3\t$Cultivar\t$Timepoint\n"
  elif [ $i == '9' ]
  then
      printf "GD_C3_3\t$Cultivar\t$Timepoint\n"
    elif [ $i == '19' ]
    then
        printf "M9_C1_3\t$Cultivar\t$Timepoint\n"
    elif [ $i == '20' ]
    then
        printf "M9_C2_3\t$Cultivar\t$Timepoint\n"
    elif [ $i == '21' ]
    then
        printf "M9_C3_3\t$Cultivar\t$Timepoint\n"

      elif [ $i == '1' ]
      then
          printf "GD_4A4\t$Cultivar\t$Timepoint\n"
      elif [ $i == '4' ]
      then
          printf "GD_5A3\t$Cultivar\t$Timepoint\n"
      elif [ $i == '6' ]
      then
          printf "GD_6A3\t$Cultivar\t$Timepoint\n"
        elif [ $i == '14' ]
        then
            printf "M9_2A3\t$Cultivar\t$Timepoint\n"
        elif [ $i == '16' ]
        then
            printf "M9_5A4\t$Cultivar\t$Timepoint\n"
        elif [ $i == '18' ]
        then
            printf "M9_6A3\t$Cultivar\t$Timepoint\n"
          elif [ $i == '2' ]
          then
              printf "GD_4A5\t$Cultivar\t$Timepoint\n"
          elif [ $i == '3' ]
          then
              printf "GD_5A2\t$Cultivar\t$Timepoint\n"
          elif [ $i == '5' ]
          then
              printf "GD_6A2\t$Cultivar\t$Timepoint\n"
            elif [ $i == '13' ]
            then
                printf "M9_2A2\t$Cultivar\t$Timepoint\n"
            elif [ $i == '15' ]
            then
                printf "M9_5A2\t$Cultivar\t$Timepoint\n"
            elif [ $i == '17' ]
            then
                printf "M9_6A2\t$Cultivar\t$Timepoint\n"
  fi
done >> $OutDir/N.dit_Hg199_RNAseq_design.txt

# Edit header lines of feature counts files to ensure they have the treatment name rather than file name
OutDir=alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq
mkdir -p $OutDir
for File in $(ls alignment/star/N.ditissima/Hg199_minion/cultivars/*/*/*_featurecounts.txt)
do
    echo $File
    cp $File $OutDir/.
done
for File in $(ls $OutDir/*_featurecounts.txt)
do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g')
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done


#Create a gene sequence fastafile
gffread gene_models_20170612.gff -g GDDH13_1-1_formatted.fasta -w your_transcripts.fasta
cat your_transcripts.fasta | sed 's/gene=.*//g' > transcripts_modified.fasta

#DeSeq commands

```R
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("data.table", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treatment you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
indexes <- c("GD_4A4","GD_4A5","GD_5A2","GD_5A3","GD_6A2","GD_6A3","GD_C1_3","GD_C2_3","GD_C3_3","Hg199_1","Hg199_2", "Hg199_3","M9_2A2","M9_2A3","M9_5A2","M9_5A4","M9_6A2","M9_6A3","M9_C1_3","M9_C2_3","M9_C3_3")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/Hg199_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/Hg199_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/Hg199_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])

colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)
countData <- round(countData,0)

#Eliminate samples
colData <- colData[!(colData$Sample=="Hg199_1"),]      
countData <- subset(countData, select=-Hg199_1)
colData <- colData[!(colData$Sample=="Hg199_2"),]      
countData <- subset(countData, select=-Hg199_2)
colData <- colData[!(colData$Sample=="Hg199_3"),]      
countData <- subset(countData, select=-Hg199_3)

colData <- colData[!(colData$Sample=="M9_C1_3"),]      
countData <- subset(countData, select=-M9_C1_3)
colData <- colData[!(colData$Sample=="M9_C2_3"),]      
countData <- subset(countData, select=-M9_C2_3)
colData <- colData[!(colData$Sample=="M9_C3_3"),]      
countData <- subset(countData, select=-M9_C3_3)

colData <- colData[!(colData$Sample=="M9_2A2"),]      
countData <- subset(countData, select=-M9_2A2)
colData <- colData[!(colData$Sample=="M9_6A2"),]      
countData <- subset(countData, select=-M9_6A2)
colData <- colData[!(colData$Sample=="M9_5A2"),]      
countData <- subset(countData, select=-M9_5A2)

colData <- colData[!(colData$Sample=="M9_2A3"),]      
countData <- subset(countData, select=-M9_2A3)
colData <- colData[!(colData$Sample=="M9_5A4"),]      
countData <- subset(countData, select=-M9_5A4)
colData <- colData[!(colData$Sample=="M9_6A3"),]      
countData <- subset(countData, select=-M9_6A3)

colData <- colData[!(colData$Sample=="GD_C1_3"),]      
countData <- subset(countData, select=-GD_C1_3)
colData <- colData[!(colData$Sample=="GD_C2_3"),]      
countData <- subset(countData, select=-GD_C2_3)
colData <- colData[!(colData$Sample=="GD_C3_3"),]      
countData <- subset(countData, select=-GD_C3_3)

colData <- colData[!(colData$Sample=="GD_4A4"),]      
countData <- subset(countData, select=-GD_4A4)
colData <- colData[!(colData$Sample=="GD_6A3"),]      
countData <- subset(countData, select=-GD_6A3)
colData <- colData[!(colData$Sample=="GD_5A3"),]      
countData <- subset(countData, select=-GD_5A3)

colData <- colData[!(colData$Sample=="GD_5A2"),]      
countData <- subset(countData, select=-GD_5A2)
colData <- colData[!(colData$Sample=="GD_6A2"),]      
countData <- subset(countData, select=-GD_6A2)
colData <- colData[!(colData$Sample=="GD_4A5"),]      
countData <- subset(countData, select=-GD_4A5)

design <- ~Group
#design <- colData$Group

dds <- DESeqDataSetFromMatrix(countData,colData,design)
#sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("iterate")))
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")


#  #Run DESeq2 removing an outlier
#
#  library(DESeq2)
#  colData <- read.table("colData",header=T,sep="\t")
#  countData <- read.table("countData2",header=T,sep="\t")
#
#  colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)
#  #Eliminate Frq08_DD24_rep3 sample from colData and countData
#  colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]
#  countData <- subset(countData, select=-Frq08_DD24_rep3)
#
#  design <- ~Group
#  dds <-  DESeqDataSetFromMatrix(countData,colData,design)
#  sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
#  dds <- DESeq(dds, fitType="local")
#
#Sample Distances

library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
  trace="none",  # turns off trace lines inside the heat map
  col=colours, # use on color palette defined earlier
  margins=c(12,12), # widens margins around plot
  srtCol=45,
  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
dev.off()

pdf("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#GD

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 412 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 382, 93%
LFC < 0 (down)   : 30, 7.3%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t2_down.txt",sep="\t",na="",quote=F)

out of 5337 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 4058, 76%
LFC < 0 (down)   : 1279, 24%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t1_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 18903 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 9706, 51%
LFC < 0 (down)   : 9197, 49%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t1_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 18677 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 8620, 46%
LFC < 0 (down)   : 10057, 54%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","GD_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/GD_t2_vs_GD_t0_down.txt",sep="\t",na="",quote=F)

out of 15916 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 7590, 48%
LFC < 0 (down)   : 8326, 52%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","M9_t0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/M9_t2_vs_M9_t0_down.txt",sep="\t",na="",quote=F)

out of 11146 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 6156, 55%
LFC < 0 (down)   : 4990, 45%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)
tar xvjf filename.tar.bz2

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library("naturalsort",lib.loc="/home/adamst/R/x86_64-pc-linux-gnu-library/3.2/")
library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("/data/scratch/gomeza/apple_genome/transcripts_modified.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/cultivars/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)
```


#Produce a detailed table of analyses
'''python
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

Run with commands:

```bash
for GeneGff in $(ls apple_genome/gene_models_20170612.gff3); do
    #Strain=JR2
    #Organism=V.dahliae
    Assembly=$(ls apple_genome/GDDH13_1-1_formatted.fasta)
    InterPro=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/interproscan/V.dahliae/JR2/JR2_interproscan.tsv)
    SwissProt=$(ls /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/gene_pred/swissprot/V.dahliae/12008/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls public_genomes/JR2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa)
    Dir1=$(ls -d RNA_alignment/featureCounts/experiment_53)
    Dir2=$(ls -d RNA_alignment/featureCounts/experiment_53/WT53)
    DEG_Files=$(ls \
        $Dir1/Frq53_LD_06h.txt \
        $Dir1/Wc153_LD_06h.txt \
        $Dir1/Wt53_Frq53_bl06h.txt \
        $Dir1/Wt53_Frq53_d06h.txt \
        $Dir1/Wt53_bl06h_vs_d06h.txt \
        $Dir1/Wt53_Wc153_bl06h.txt \
        $Dir1/Wt53_Wc153_d06h.txt \
        $Dir2/Wt53_d06h_d12h.txt \
        $Dir2/Wt53_d06h_d18h.txt \
        $Dir2/Wt53_d06h_d24h.txt \
        $Dir2/Wt53_d12h_d18h.txt \
        $Dir2/Wt53_d12h_d24h.txt \
        $Dir2/Wt53_d24h_d18h.txt \
        $Dir2/Wt53_LD_06h.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts_53.txt)
    FPKM=$(ls $Dir1/countData_53.fpkm)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_incl_exp.tsv
done
```

#Inital analysis of tables of DEGs

```bash
  effector_names=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_secreted_headers.txt
  CAZY_names=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
  for File in $(ls alignment/star/N.ditissima/Hg199_minion/DeSeq/GD*.txt | grep -v "genes" | grep -v "countData")
  do
      Assessment=$(basename $File | sed "s/.txt//g")
      echo $Assessment
      echo "Total number of genes in dataset:"
      cat $File | grep -v 'baseMean' | wc -l
      echo "Total number of secretedeffector in dataset:"
      Effector_File=$(echo $File | sed "s/.txt/_Effector.txt/g")
      cat $File | head -n 1 > $Effector_File
      cat $File | grep -w -f $effector_names >> $Effector_File
      cat $Effector_File | tail -n +2 | wc -l
        echo "Total number of CAZY in dataset:"
        CAZY_File=$(echo $File | sed "s/.txt/_CAZY.txt/g")
        cat $File | head -n 1 > $CAZY_File
        cat $File | grep -w -f $CAZY_names >> $CAZY_File
        cat $CAZY_File | tail -n +2 | wc -l

    done
done
```


##Produce a more detailed table of analyses

```bash
for Strain in Bc1 Nov9
do
    for GeneGff in $(ls gene_pred/annotation/P.fragariae/$Strain/"$Strain"_genes_incl_ORFeffectors.gff3)
    do
        Organism=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
        Assembly=$(ls repeat_masked/P.fragariae/$Strain/ncbi_edits_repmask/*_contigs_unmasked.fa)
        InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/greedy/*_interproscan.tsv)
        SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/greedy/swissprot_vJul2016_tophit_parsed.tbl)
        OutDir=gene_pred/annotation/$Organism/$Strain
        mkdir -p $OutDir
        # GeneFasta=$(ls gene_pred/annotation/P.cactorum/414_v2/414_v2_genes_incl_ORFeffectors.pep.fasta)
        GeneFasta=$(ls gene_pred/annotation/P.fragariae/$Strain/"$Strain"_genes_incl_ORFeffectors.cds.fasta)
        SigP2=$(ls gene_pred/final_sigP/$Organism/$Strain/*_aug_sp.aa)
        SigP2_ORF=$(ls gene_pred/ORF_sigP/$Organism/$Strain/*_aug_sp.aa)
        SigP3=$(ls gene_pred/final_signalp-3.0/$Organism/$Strain/*_aug_sp.aa)
        SigP3_ORF=$(ls gene_pred/ORF_signalp-3.0/$Organism/$Strain/*_aug_sp.aa)
        SigP4=$(ls gene_pred/final_signalp-4.1/$Organism/$Strain/*_aug_sp.aa)
        SigP4_ORF=$(ls gene_pred/ORF_signalp-4.1/$Organism/$Strain/*_aug_sp.aa)
        TMHMM_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/greedy/*_TM_genes_pos_headers.txt)
        GPI_headers=$(ls gene_pred/GPIsom/$Organism/$Strain/greedy/GPI_pos.txt)
        PhobiusFa=$(ls analysis/phobius_CQ/$Organism/$Strain/*_phobius.fa)
        PhobiusFa_ORF=$(ls analysis/phobius_ORF/$Organism/$Strain/*_phobius.fa)
        #RxLR_Motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/*_RxLR_EER_regex.fa | grep -v 'ORF')
        #RxLR_Hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer.fa | grep -v 'ORF')
        #RxLR_WY=$(ls analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/*_WY_hmmer_headers.txt | grep -v 'ORF')
        RxLR_total=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_Total_RxLR_motif_hmm.txt)
        RxLR_ORF_total=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_total_ORF_RxLR_headers.txt)
        RxLR_EER_total=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_Total_RxLR_EER_motif_hmm.txt)
        RxLR_EER_ORF_total=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_total_ORF_RxLR_EER_headers.txt)
        #CRN_LFLAK=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_pub_CRN_LFLAK_hmm.fa | grep -v 'ORF')
        #CRN_DWL=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_pub_CRN_DWL_hmm.fa | grep -v 'ORF')
        CRN_total=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_final_CRN.txt)
        ApoP_total=$(ls analysis/ApoplastP/$Organism/$Strain/*_Total_ApoplastP.txt)
        #	OrthoName=Pcac
        #	OrthoFile=$(ls analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt)
        ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
        DEG_Files=$(ls alignment/star/P.fragariae/$Strain/DeSeq/*_vs_*.txt  | grep -v -e 'up' -e 'down' -e "CRN" -e "RxLR" | sed -e "s/$/ /g" | tr -d "\n")
        # $ProgDir/pacbio_anntoation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP4 $SigP4 --phobius $PhobiusTxt --RxLR_motif $RxLR_Motif --RxLR_Hmm $RxLR_Hmm --RxLR_WY $RxLR_WY --RxLR_total $RxLR_total --CRN_LFLAK $CRN_LFLAK --CRN_DWL $CRN_DWL --CRN_total $CRN_total --DEG_files $DEG_Files  > $OutDir/414_v2_gene_table_incl_exp.tsv
        # NormCount=$(ls alignment/star/P.cactorum/414_v2/DeSeq/normalised_counts.txt)
        RawCount=$(ls alignment/star/P.fragariae/$Strain/DeSeq/raw_counts.txt)
        FPKM=$(ls alignment/star/P.fragariae/$Strain/DeSeq/fpkm_counts.txt)
        OrthoName=$Strain
        OrthoFile=analysis/orthology/orthomcl/All_Strains_plus_rubi_no_removal/All_Strains_plus_rubi_no_removal_orthogroups.txt
        $ProgDir/pacbio_anntoation_tables_modified.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP2_ORF $SigP2_ORF --SigP3 $SigP3 --SigP3_ORF $SigP3_ORF --SigP4 $SigP4 --SigP4_ORF $SigP4_ORF --phobius $PhobiusFa --phobius_ORF $PhobiusFa_ORF --trans_mem $TMHMM_headers --GPI_anchor $GPI_headers --RxLR_total $RxLR_total --RxLR_total_ORF $RxLR_ORF_total --RxLR_EER_total $RxLR_EER_total --RxLR_EER_total_ORF $RxLR_EER_ORF_total --CRN_total $CRN_total --ApoP_total $ApoP_total --ortho_name $OrthoName --ortho_file $OrthoFile --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp.tsv
    done
done
```

#Extract fasta file of all DEGs for BLAST analysis

```bash
for Strain in Bc1 Nov9
do
    input=alignment/star/P.fragariae/$Strain/DeSeq/*_vs_*mycelium.txt
    DEGFile=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs.tsv
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
    $ProgDir/parse_RNA-Seq_1_timepoint.py --input $input --out_dir $DEGFile
    DEGNames=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs_names.txt
    Genes=gene_pred/annotation/P.fragariae/$Strain/"$Strain"_genes_incl_ORFeffectors.cds.fasta
    DEGFasta=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs.fa
    $ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
done
```

#Investigate enriched functional annotations in DEGs vs all genes

##Analysis of DEGs vs all genes

```bash
for Strain in Bc1 Nov9
do
    OutDir=analysis/enrichment/P.fragariae/$Strain/Whole_Genome
    mkdir -p $OutDir
    InterProTSV=gene_pred/interproscan/P.fragariae/$Strain/greedy/"$Strain"_interproscan.tsv
    ProgDir=/home/adamst/git_repos/scripts/fusarium/analysis/gene_enrichment
    $ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/"$Strain"_gene_GO_annots.tsv

    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
    AnnotTable=gene_pred/annotation/P.fragariae/$Strain/"$Strain"_gene_table_incl_exp.tsv
    DEGs=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs_names.txt
    AllGenes=$OutDir/"$Strain"_all_genes.txt
    cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
    Set1Genes=$OutDir/"$Strain"_DEGs.txt
    Set2Genes=$OutDir/"$Strain"_all_genes2.txt
    AllGenes=$OutDir/"$Strain"_all_genes.txt
    cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
    cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
    cat $Set1Genes $Set2Genes > $AllGenes

    $ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/"$Strain"_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
done
```

Method 2

```
#DeSeq commands

```R
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("data.table", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treatment you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
indexes <- c("GD_4A4","GD_4A5","GD_5A2","GD_5A3","GD_6A2","GD_6A3","GD_C1_3","GD_C2_3","GD_C3_3","Hg199_1","Hg199_2", "Hg199_3","M9_2A2","M9_2A3","M9_5A2","M9_5A4","M9_6A2","M9_6A3","M9_C1_3","M9_C2_3","M9_C3_3")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/Hg199_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/Hg199_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/N.dit_Hg199_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/Hg199_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])

colData$Group <- paste0(colData$Cultivar,'_', colData$Timepoint)
countData <- round(countData,0)

#Eliminate samples
colData <- colData[!(colData$Sample=="Hg199_1"),]      
countData <- subset(countData, select=-Hg199_1)
colData <- colData[!(colData$Sample=="Hg199_2"),]      
countData <- subset(countData, select=-Hg199_2)
colData <- colData[!(colData$Sample=="Hg199_3"),]      
countData <- subset(countData, select=-Hg199_3)

colData <- colData[!(colData$Sample=="M9_C1_3"),]      
countData <- subset(countData, select=-M9_C1_3)
colData <- colData[!(colData$Sample=="M9_C2_3"),]      
countData <- subset(countData, select=-M9_C2_3)
colData <- colData[!(colData$Sample=="M9_C3_3"),]      
countData <- subset(countData, select=-M9_C3_3)

colData <- colData[!(colData$Sample=="M9_2A2"),]      
countData <- subset(countData, select=-M9_2A2)
colData <- colData[!(colData$Sample=="M9_6A2"),]      
countData <- subset(countData, select=-M9_6A2)
colData <- colData[!(colData$Sample=="M9_5A2"),]      
countData <- subset(countData, select=-M9_5A2)

colData <- colData[!(colData$Sample=="M9_2A3"),]      
countData <- subset(countData, select=-M9_2A3)
colData <- colData[!(colData$Sample=="M9_5A4"),]      
countData <- subset(countData, select=-M9_5A4)
colData <- colData[!(colData$Sample=="M9_6A3"),]      
countData <- subset(countData, select=-M9_6A3)

colData <- colData[!(colData$Sample=="GD_C1_3"),]      
countData <- subset(countData, select=-GD_C1_3)
colData <- colData[!(colData$Sample=="GD_C2_3"),]      
countData <- subset(countData, select=-GD_C2_3)
colData <- colData[!(colData$Sample=="GD_C3_3"),]      
countData <- subset(countData, select=-GD_C3_3)

colData <- colData[!(colData$Sample=="GD_4A4"),]      
countData <- subset(countData, select=-GD_4A4)
colData <- colData[!(colData$Sample=="GD_6A3"),]      
countData <- subset(countData, select=-GD_6A3)
colData <- colData[!(colData$Sample=="GD_5A3"),]      
countData <- subset(countData, select=-GD_5A3)

colData <- colData[!(colData$Sample=="GD_5A2"),]      
countData <- subset(countData, select=-GD_5A2)
colData <- colData[!(colData$Sample=="GD_6A2"),]      
countData <- subset(countData, select=-GD_6A2)
colData <- colData[!(colData$Sample=="GD_4A5"),]      
countData <- subset(countData, select=-GD_4A5)

design <- ~Group
#design <- colData$Group

dds <- DESeqDataSetFromMatrix(countData,colData,design)
#sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("iterate")))
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")


#  #Run DESeq2 removing an outlier
#
#  library(DESeq2)
#  colData <- read.table("colData",header=T,sep="\t")
#  countData <- read.table("countData2",header=T,sep="\t")
#
#  colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)
#  #Eliminate Frq08_DD24_rep3 sample from colData and countData
#  colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]
#  countData <- subset(countData, select=-Frq08_DD24_rep3)
#
#  design <- ~Group
#  dds <-  DESeqDataSetFromMatrix(countData,colData,design)
#  sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
#  dds <- DESeq(dds, fitType="local")
#
#Sample Distances

library("RColorBrewer")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
  trace="none",  # turns off trace lines inside the heat map
  col=colours, # use on color palette defined earlier
  margins=c(12,12), # widens margins around plot
  srtCol=45,
  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plotsPl

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Cultivar", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Cultivar", "Timepoint"))
dev.off()

pdf("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#GD

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 4762 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2143, 45%
LFC < 0 (down)   : 2619, 55%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 5232 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1619, 31%
LFC < 0 (down)   : 3613, 69%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 792 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 558, 70%
LFC < 0 (down)   : 234, 30%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","GD_t2","mycelium_Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control_down.txt",sep="\t",na="",quote=F)

out of 1987 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 775, 39%
LFC < 0 (down)   : 1212, 61%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t1","GD_t1"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_GD_t1.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_GD_t1_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_GD_t1_down.txt",sep="\t",na="",quote=F)

out of 53 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 5, 9.4%
LFC < 0 (down)   : 48, 91%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","M9_t2","GD_t2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_GD_t2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_GD_t2_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_GD_t2_down.txt",sep="\t",na="",quote=F)

out of 1 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 1, 100%
LFC < 0 (down)   : 0, 0%
outliers [1]     : 0, 0%
low counts [2]   : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library("naturalsort",lib.loc="/home/adamst/R/x86_64-pc-linux-gnu-library/3.2/")
library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("/data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)
```


#Inital analysis of tables of DEGs

```bash
  effector_names=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
  CAZY_names=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
  for File in $(ls alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/*vs*.txt | grep -v "genes" | grep -v "countData")
  do
      Assessment=$(basename $File | sed "s/.txt//g")
      echo $Assessment
      echo "Total number of genes in dataset:"
      cat $File | grep -v 'baseMean' | wc -l
      echo "Total number of effectors in dataset:"
      Effector_File=$(echo $File | sed "s/.txt/_Effector.txt/g")
      cat $File | head -n 1 > $Effector_File
      cat $File | grep -w -f $effector_names >> $Effector_File
      cat $Effector_File | tail -n +2 | wc -l
      echo "Total number of CAZY in dataset:"
      CAZY_File=$(echo $File | sed "s/.txt/_CAZY.txt/g")
      cat $File | head -n 1 > $CAZY_File
      cat $File | grep -w -f $CAZY_names >> $CAZY_File
      cat $CAZY_File | tail -n +2 | wc -l
    done
done
```

GD_t1_vs_mycelium_Control_down
Total number of genes in dataset:
2382
Total number of secretedeffector in dataset:
479
Total number of CAZY in dataset:
91
GD_t1_vs_mycelium_Control
Total number of genes in dataset:
4762
Total number of secretedeffector in dataset:
777
Total number of CAZY in dataset:
335
GD_t1_vs_mycelium_Control_up
Total number of genes in dataset:
2070
Total number of secretedeffector in dataset:
233
Total number of CAZY in dataset:
242

GD_t2_vs_mycelium_Control_down
Total number of genes in dataset:
1207
Total number of secretedeffector in dataset:
189
Total number of CAZY in dataset:
56
GD_t2_vs_mycelium_Control
Total number of genes in dataset:
1987
Total number of secretedeffector in dataset:
283
Total number of CAZY in dataset:
191
GD_t2_vs_mycelium_Control_up
Total number of genes in dataset:
775
Total number of secretedeffector in dataset:
91
Total number of CAZY in dataset:
135

M9_t1_vs_GD_t1_down
Total number of genes in dataset:
47
Total number of secretedeffector in dataset:
6
Total number of CAZY in dataset:
3
M9_t1_vs_GD_t1
Total number of genes in dataset:
53
Total number of secretedeffector in dataset:
6
Total number of CAZY in dataset:
6
M9_t1_vs_GD_t1_up
Total number of genes in dataset:
5
Total number of secretedeffector in dataset:
0
Total number of CAZY in dataset:
3

M9_t1_vs_mycelium_Control_down
Total number of genes in dataset:
3343
Total number of secretedeffector in dataset:
587
Total number of CAZY in dataset:
122
M9_t1_vs_mycelium_Control
Total number of genes in dataset:
5232
Total number of secretedeffector in dataset:
821
Total number of CAZY in dataset:
344
M9_t1_vs_mycelium_Control_up
Total number of genes in dataset:
1594
Total number of secretedeffector in dataset:
175
Total number of CAZY in dataset:
219

M9_t2_vs_mycelium_Control_down
Total number of genes in dataset:
234
Total number of secretedeffector in dataset:
44
Total number of CAZY in dataset:
19
M9_t2_vs_mycelium_Control
Total number of genes in dataset:
792
Total number of secretedeffector in dataset:
120
Total number of CAZY in dataset:
114
M9_t2_vs_mycelium_Control_up
Total number of genes in dataset:
558
Total number of secretedeffector in dataset:
76
Total number of CAZY in dataset:
95
```
```
#Draw venn diagrams of differentially expressed genes

##All genes

###All DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/t1_all_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```


###Upregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_up.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control_up.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_up.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control_up.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_up.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_up.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/t1_up_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Downregulated DEGs

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_down.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t2_vs_mycelium_Control_down.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_down.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t2_vs_mycelium_Control_down.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
inp1=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/GD_t1_vs_mycelium_Control_down.txt
inp2=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/M9_t1_vs_mycelium_Control_down.txt
OutDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/t1_down_DEGs.tsv
$ProgDir/parse_RNA-Seq.py --input_1 $inp1 --input_2 $inp2 --out_dir $OutDir
```

###Venn diagrams

```bash
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis
WorkDir=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/GD_all_DEGs.tsv --out $WorkDir/GD_all_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_all_DEGs.tsv --out $WorkDir/t1_all_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_up_DEGs.tsv --out $WorkDir/t1_up_DEGs.pdf
$ProgDir/All_DEGs_venn_diag_2.r --inp $WorkDir/t1_down_DEGs.tsv --out $WorkDir/t1_down_DEGs.pdf
```



##Produce a more detailed table of analyses

```bash
    for GeneGff in $(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3)
    do
      Strain=Hg199_minion
      Organism=N.ditissima
        Assembly=$(ls /data/scratch/gomeza/Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/*_contigs_unmasked.fa)
        InterPro=$(ls /data/scratch/gomeza/gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv)
        SwissProt=$(ls /data/scratch/gomeza/gene_pred/swissprot/N.ditissima/Hg199_minion/swissprot_vJul2016_tophit_parsed.tbl)
        OutDir=gene_pred/annotation/$Organism/$Strain
        mkdir -p $OutDir
        GeneFasta=$(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.pep.fasta)

        SigP4=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/Hg199_minion_final_sp_no_trans_mem.aa)
        effector_total=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
        CAZY_total=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt

        TMHMM_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_pos_headers.txt)

        ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis

        Dir1=$(ls -d /data/scratch/gomeza/alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq)
        DEG_Files=$(ls \
            $Dir1/GD_t1_vs_mycelium_Control.txt \
            $Dir1/M9_t1_vs_mycelium_Control.txt \
            $Dir1/GD_t2_vs_mycelium_Control.txt \
            $Dir1/M9_t2_vs_mycelium_Control.txt \
            $Dir1/M9_t1_vs_GD_t1.txt \
            | sed -e "s/$/ /g" | tr -d "\n")

        #DEG_Files=$(ls alignment/star/N.ditisima/Hg199_minion/unmapped_to_Nd/DeSeq/*_vs_*.txt  | grep -v -e 'up' -e 'down' -e "CAZY" -e "Effector" | sed -e "s/$/ /g" | tr -d "\n")
        # $ProgDir/pacbio_anntoation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP4 $SigP4 --phobius $PhobiusTxt --RxLR_motif $RxLR_Motif --RxLR_Hmm $RxLR_Hmm --RxLR_WY $RxLR_WY --RxLR_total $RxLR_total --CRN_LFLAK $CRN_LFLAK --CRN_DWL $CRN_DWL --CRN_total $CRN_total --DEG_files $DEG_Files  > $OutDir/414_v2_gene_table_incl_exp.tsv
        # NormCount=$(ls alignment/star/P.cactorum/414_v2/DeSeq/normalised_counts.txt)
        RawCount=$(ls alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/raw_counts.txt )
        FPKM=$(ls alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/fpkm_counts.txt)
        #OrthoName=$Strain
      #OrthoFile=analysis/orthology/orthomcl/All_Strains_plus_rubi_no_removal/All_Strains_plus_rubi_no_removal_orthogroups.txt
        $ProgDir/pacbio_anntoation_tables_modified.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP4 $SigP4 --trans_mem $TMHMM_headers --effector_total $effector_total --CAZY_total $CAZY_total --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp.tsv
    done
done
```



#Produce a detailed table of analyses
'''python
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

Run with commands:

```bash
for GeneGff in $(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_appended.gff3); do
    Strain=Hg199_minion
    Organism=N.ditissima
    Assembly=$(ls /data/scratch/gomeza/Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    InterPro=$(ls /data/scratch/gomeza/gene_pred/interproscan/N.ditissima/Hg199_minion/Hg199_minion_interproscan.tsv)
    SwissProt=$(ls /data/scratch/gomeza/gene_pred/swissprot/N.ditissima/Hg199_minion/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    GeneFasta=$(ls /data/scratch/gomeza/gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.pep.fasta)


    Dir1=$(ls -d /data/scratch/gomeza/alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq)
    DEG_Files=$(ls \
        $Dir1/GD_t1_vs_mycelium_Control.txt \
        $Dir1/M9_t1_vs_mycelium_Control.txt \
        $Dir1/GD_t2_vs_mycelium_Control.txt \
        $Dir1/M9_t2_vs_mycelium_Control.txt \
        $Dir1/M9_t1_vs_GD_t1.txt \
        | sed -e "s/$/ /g" | tr -d "\n")

    RawCount=$(ls $Dir1/raw_counts.txt)
    FPKM=$(ls $Dir1/fpkm_counts.txt)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_clocks/annotation
    $ProgDir/Vd_annotation_tables.py --gene_gff $GeneGff --gene_fasta $GeneFasta --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_gene_table_incl_exp.tsv
done
```
```





























#Extract fasta file of all DEGs for BLAST analysis

```bash
for Cultivar in GD M9
do
    input=alignment/star/N.ditisima/Hg199_minion/unmapped_to_Nd/DeSeq/*_vs_mycelium_Control.txt
    DEGFile=alignment/star/N.ditissima/Hg199_minion/unmapped_to_Nd/DeSeq/"$Cultivar"_all_DEGs.tsv
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
    $ProgDir/parse_RNA-Seq_1_timepoint.py --input $input --out_dir $DEGFile

    DEGNames=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs_names.txt
    Genes=gene_pred/annotation/P.fragariae/$Strain/"$Strain"_genes_incl_ORFeffectors.cds.fasta
    DEGFasta=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs.fa
    $ProgDir/extract_DEG_Names_1_timepoint.py --input $DEGFile --output $DEGNames
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Genes --headers $DEGNames > $DEGFasta
done
```

#Investigate enriched functional annotations in DEGs vs all genes

##Analysis of DEGs vs all genes

```bash
for Strain in Bc1 Nov9
do
    OutDir=analysis/enrichment/P.fragariae/$Strain/Whole_Genome
    mkdir -p $OutDir
    InterProTSV=gene_pred/interproscan/P.fragariae/$Strain/greedy/"$Strain"_interproscan.tsv
    ProgDir=/home/adamst/git_repos/scripts/fusarium/analysis/gene_enrichment
    $ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/"$Strain"_gene_GO_annots.tsv

    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
    AnnotTable=gene_pred/annotation/P.fragariae/$Strain/"$Strain"_gene_table_incl_exp.tsv
    DEGs=alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"_all_DEGs_names.txt
    AllGenes=$OutDir/"$Strain"_all_genes.txt
    cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
    Set1Genes=$OutDir/"$Strain"_DEGs.txt
    Set2Genes=$OutDir/"$Strain"_all_genes2.txt
    AllGenes=$OutDir/"$Strain"_all_genes.txt
    cat $DEGs | sed -e 's/$/\t0.001/g' > $Set1Genes
    cat $AnnotTable | tail -n+2 | cut -f1 | grep -v $Set1Genes | sed -e 's/$/\t1.00/g' > $Set2Genes
    cat $Set1Genes $Set2Genes > $AllGenes

    $ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/"$Strain"_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
done
```
