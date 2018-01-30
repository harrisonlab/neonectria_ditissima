
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

mv qc_rna/Bc1/ qc_rna/novogene/P.fragariae/.
mv qc_rna/Nov9/ qc_rna/novogene/P.fragariae/.
```

###Visualise data quality using fastqc

Only submit three jobs at a time, copying 30 files is too much!

```bash
for RawData in $(ls qc_rna/N.ditissima/M9/t1/*/*)
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

#Align mycelium reads to BC-1 & NOV-9 assemblies with STAR

```bash
#BC-1
for Assembly in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/P.fragariae/Bc1/ncbi_edits_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_rna/novogene/P.fragariae/Bc1/mycelium/F/*_trim.fq.gz)
    do
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
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
        qsub $ProgDir/sub_star_TA.sh $Assembly $FileF $FileR $OutDir
    done
done

#NOV-9
for Assembly in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/P.fragariae/Nov9/ncbi_edits_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_rna/novogene/P.fragariae/Nov9/mycelium/F/*_trim.fq.gz)
    do
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
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
        qsub $ProgDir/sub_star_TA.sh $Assembly $FileF $FileR $OutDir
    done
done
```

##Align all timepoints to *Fragaria vesca* genome v1.1

```bash
for FileF in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_rna/novogene/P.fragariae/*/*/F/*_trim.fq.gz | grep -e "Bc1" -e "Nov9" | grep -v "Bc16" | grep -v "mycelium")
do
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
    OutDir=alignment/star/vesca_alignment/set2/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    Assembly=/home/sobczm/popgen/rnaseq/fvesca_v1.1_all.fa
    GFF=/home/sobczm/popgen/rnaseq/Fragaria_vesca_v1.1.a2.gff3
    qsub $ProgDir/sub_star_sensitive.sh $Assembly $FileF $FileR $OutDir $GFF
done
```

##Gzip output files to save space on the disk and allow star to run correctly downstream. ONLY RUN THIS ONCE

```bash
for AlignDir in $(ls -d /home/groups/harrisonlab/project_files/phytophthora_fragariae/alignment/star/vesca_alignment/set2/*/*)
do
    cat $AlignDir/star_aligmentUnmapped.out.mate1 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate1.fq.gz
    cat $AlignDir/star_aligmentUnmapped.out.mate2 | gzip -cf > $AlignDir/star_aligmentUnmapped.out.mate2.fq.gz
done
```

##Align compressed files of unmapped reads from aligning to vesca

This star script had the following options added to the sub_star.sh script in the ProgDir specified in the below commands:
--winAnchorMultimapNmax 200
--seedSearchStartLmax 30

```bash
#BC-1
for AlignDir in $(ls -d /home/groups/harrisonlab/project_files/phytophthora_fragariae/alignment/star/vesca_alignment/set2/48hr/*)
do
    Organism=P.fragariae
    Strain=Bc1
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
    if [ -f /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    else
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    fi
    OutDir=alignment/star/P.fragariae/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    qsub $ProgDir/sub_star_TA.sh $Assembly $File1 $File2 $OutDir
done

#NOV-9
for AlignDir in $(ls -d /home/groups/harrisonlab/project_files/phytophthora_fragariae/alignment/star/vesca_alignment/set2/72hr/*)
do
    Organism=P.fragariae
    Strain=Nov9
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
    if [ -f /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    else
        Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    fi
    OutDir=alignment/star/P.fragariae/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    qsub $ProgDir/sub_star_TA.sh $Assembly $File1 $File2 $OutDir
done
```

#Quantification of gene models

```bash
for Strain in Bc1 Nov9
do
    for BamFile in $(ls alignment/star/P.fragariae/$Strain/*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v "TA-")
    do
        Gff=gene_pred/annotation/P.fragariae/$Strain/*_genes_incl_ORFeffectors.gff3
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
done
```

##A file was created with columns referring to experimental treatments

```bash
#BC-1
OutDir=alignment/star/P.fragariae/Bc1/DeSeq
mkdir -p $OutDir
printf "Sample.name\tTimepoint\tIsolate\n" > $OutDir/P.frag_Bc1_RNAseq_design.txt
# for File in $(ls alignment/star/P.cactorum/10300/Sample_*/Sample_*_featurecounts.txt); do
# Sample=$(echo $File | rev | cut -f2 -d '/' | rev)
# i=$(echo $Sample | sed 's/Sample_//g')
for i in $(seq 1 6)
do
    if [ $i == '1' ] || [ $i == '2' ] || [ $i == '3' ]
    then
        Timepoint='48hr'
    elif [ $i == '4' ] || [ $i == '5' ] || [ $i == '6' ]
    then
        Timepoint='mycelium'
    fi
    Infection=Bc1
    if [ $i == '1' ]
    then
        printf "TA_B_P1\t$Timepoint\t$Infection\n"
    elif [ $i == '2' ]
    then
        printf "TA_B_P2\t$Timepoint\t$Infection\n"
    elif [ $i == '3' ]
    then
        printf "TA_B_P3\t$Timepoint\t$Infection\n"
    elif [ $i == '4' ]
    then
        printf "TA_B_M1\t$Timepoint\t$Infection\n"
    elif [ $i == '5' ]
    then
        printf "TA_B_M2\t$Timepoint\t$Infection\n"
    elif [ $i == '6' ]
    then
        printf "TA_B_M3\t$Timepoint\t$Infection\n"
    fi
done >> $OutDir/P.frag_Bc1_RNAseq_design.txt

# Edit header lines of feature counts files to ensure they have the treatment name rather than file name
OutDir=alignment/star/P.fragariae/Bc1/DeSeq
mkdir -p $OutDir
for File in $(ls alignment/star/P.fragariae/Bc1/*/*/*_featurecounts.txt)
do
    echo $File
    cp $File $OutDir/.
done
for File in $(ls $OutDir/*_featurecounts.txt)
do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g')
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done

#NOV-9
OutDir=alignment/star/P.fragariae/Nov9/DeSeq
mkdir -p $OutDir
printf "Sample.name\tTimepoint\tIsolate\n" > $OutDir/P.frag_Nov9_RNAseq_design.txt
# for File in $(ls alignment/star/P.cactorum/10300/Sample_*/Sample_*_featurecounts.txt); do
# Sample=$(echo $File | rev | cut -f2 -d '/' | rev)
# i=$(echo $Sample | sed 's/Sample_//g')
for i in $(seq 1 6)
do
    if [ $i == '1' ] || [ $i == '2' ] || [ $i == '3' ]
    then
        Timepoint='72hr'
    elif [ $i == '4' ] || [ $i == '5' ] || [ $i == '6' ]
    then
        Timepoint='mycelium'
    fi
        Infection='Nov9'
    if [ $i == '1' ]
    then
        printf "TA_NO_P1\t$Timepoint\t$Infection\n"
    elif [ $i == '2' ]
    then
        printf "TA_NO_P2\t$Timepoint\t$Infection\n"
    elif [ $i == '3' ]
    then
        printf "TA_NO_P3\t$Timepoint\t$Infection\n"
    elif [ $i == '4' ]
    then
        printf "TA_NO_M1\t$Timepoint\t$Infection\n"
    elif [ $i == '5' ]
    then
        printf "TA_NO_M2\t$Timepoint\t$Infection\n"
    elif [ $i == '6' ]
    then
        printf "TA_NO_M5\t$Timepoint\t$Infection\n"
    fi
done >> $OutDir/P.frag_Nov9_RNAseq_design.txt

# Edit header lines of feature counts files to ensure they have the treatment name rather than file name
OutDir=alignment/star/P.fragariae/Nov9/DeSeq
mkdir -p $OutDir
for File in $(ls alignment/star/P.fragariae/Nov9/*/*/*_featurecounts.txt)
do
    echo $File
    cp $File $OutDir/.
done
for File in $(ls $OutDir/*_featurecounts.txt)
do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g')
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done
```

#DeSeq commands

BC-1

```R
#install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.fragariae/Bc1/DeSeq","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

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
indexes <- c("TA_B_P1", "TA_B_P2", "TA_B_P3", "TA_B_M1", "TA_B_M2", "TA_B_M3")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/P.fragariae/Bc1/DeSeq/Bc1_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.fragariae/Bc1/DeSeq/Bc1_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/P.fragariae/Bc1/DeSeq/P.frag_Bc1_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/P.fragariae/Bc1/DeSeq/Bc1_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
countData <- round(countData,0)

design <- ~Group
#design <- colData$Group

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
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

pdf("alignment/star/P.fragariae/Bc1/DeSeq/heatmap_vst.pdf", width=12,height=12)
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

pdf("alignment/star/P.fragariae/Bc1/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plots

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/P.fragariae/Bc1/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/P.fragariae/Bc1/DeSeq/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.fragariae/Bc1/DeSeq/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.fragariae/Bc1/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#48hr vs mycelium

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Bc1_48hr","Bc1_mycelium"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.fragariae/Bc1/DeSeq/Bc1_48hr_vs_Bc1_mycelium.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.fragariae/Bc1/DeSeq/Bc1_48hr_vs_Bc1_mycelium_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.fragariae/Bc1/DeSeq/Bc1_48hr_vs_Bc1_mycelium_down.txt",sep="\t",na="",quote=F)

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.fragariae/Bc1/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.fragariae/Bc1/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("gene_pred/annotation/P.fragariae/Bc1/Bc1_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc1/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc1/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)
```

NOV-9

```R
#install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.fragariae/Nov9/DeSeq","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treament you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
indexes <- c("TA_NO_P1", "TA_NO_P2", "TA_NO_P3", "TA_NO_M1", "TA_NO_M2", "TA_NO_M5")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/P.fragariae/Nov9/DeSeq/Nov9_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.fragariae/Nov9/DeSeq/Nov9_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/P.fragariae/Nov9/DeSeq/P.frag_Nov9_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/P.fragariae/Nov9/DeSeq/Nov9_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
countData <- round(countData,0)

design <- ~Group
#design <- colData$Group

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
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

pdf("alignment/star/P.fragariae/Nov9/DeSeq/heatmap_vst.pdf", width=12,height=12)
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

pdf("alignment/star/P.fragariae/Nov9/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plots

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/P.fragariae/Nov9/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/P.fragariae/Nov9/DeSeq/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.fragariae/Nov9/DeSeq/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.fragariae/Nov9/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#72hr vs mycelium

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Nov9_72hr","Nov9_mycelium"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.fragariae/Nov9/DeSeq/Nov9_72hr_vs_Nov9_mycelium.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.fragariae/Nov9/DeSeq/Nov9_72hr_vs_Nov9_mycelium_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.fragariae/Nov9/DeSeq/Nov9_72hr_vs_Nov9_mycelium_down.txt",sep="\t",na="",quote=F)

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.fragariae/Nov9/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.fragariae/Nov9/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("gene_pred/annotation/P.fragariae/Nov9/Nov9_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Nov9/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Nov9/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)
```

#Inital analysis of tables of DEGs

```bash
for Strain in Bc1 Nov9
do
    RxLR_Names=analysis/RxLR_effectors/combined_evidence/P.fragariae/$Strain/"$Strain"_Total_RxLR_motif_hmm.txt
    CRN_Names=analysis/CRN_effectors/hmmer_CRN/P.fragariae/$Strain/"$Strain"_final_CRN.txt
    ApoP_Names=analysis/ApoplastP/P.fragariae/$Strain/"$Strain"_Total_ApoplastP.txt
    for File in $(ls alignment/star/P.fragariae/$Strain/DeSeq/"$Strain"*.txt | grep -v "genes" | grep -v "countData")
    do
        Assessment=$(basename $File | sed "s/.txt//g")
        echo $Assessment
        echo "Total number of genes in dataset:"
        cat $File | grep -v 'baseMean' | wc -l
        echo "Total number of RxLRs in dataset:"
        RxLR_File=$(echo $File | sed "s/.txt/_RxLRs.txt/g")
        cat $File | head -n 1 > $RxLR_File
        cat $File | grep -w -f $RxLR_Names >> $RxLR_File
        cat $RxLR_File | tail -n +2 | wc -l
        echo "Total number of CRNs in dataset:"
        CRN_File=$(echo $File | sed "s/.txt/_CRNs.txt/g")
        cat $File | head -n 1 > $CRN_File
        cat $File | grep -w -f $CRN_Names >> $CRN_File
        cat $CRN_File | tail -n +2 | wc -l
        echo "Total number of Apoplastic effectors in dataset:"
        ApoP_File=$(echo $File | sed "s/.txt/_ApoP.txt/g")
        cat $File | head -n 1 > $ApoP_File
        cat $File | grep -w -f $ApoP_Names >> $ApoP_File
        cat $ApoP_File | tail -n +2 | wc -l
    done
done
```

```
Bc1_48hr_vs_Bc1_mycelium_down
Total number of genes in dataset:
2,769
Total number of RxLRs in dataset:
100
Total number of CRNs in dataset:
22
Total number of Apoplastic effectors in dataset:
364
Bc1_48hr_vs_Bc1_mycelium
Total number of genes in dataset:
7,837
Total number of RxLRs in dataset:
258
Total number of CRNs in dataset:
31
Total number of Apoplastic effectors in dataset:
877
Bc1_48hr_vs_Bc1_mycelium_up
Total number of genes in dataset:
2,630
Total number of RxLRs in dataset:
121
Total number of CRNs in dataset:
1
Total number of Apoplastic effectors in dataset:
433

Nov9_72hr_vs_Nov9_mycelium_down
Total number of genes in dataset:
1,700
Total number of RxLRs in dataset:
74
Total number of CRNs in dataset:
11
Total number of Apoplastic effectors in dataset:
271
Nov9_72hr_vs_Nov9_mycelium
Total number of genes in dataset:
5,471
Total number of RxLRs in dataset:
211
Total number of CRNs in dataset:
18
Total number of Apoplastic effectors in dataset:
729
Nov9_72hr_vs_Nov9_mycelium_up
Total number of genes in dataset:
2,825
Total number of RxLRs in dataset:
125
Total number of CRNs in dataset:
5
Total number of Apoplastic effectors in dataset:
425
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

#Align mycelium reads to BC-16 assemblies with STAR

```bash
#BC-1
for Assembly in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa)
do
    Strain=Bc16
    Organism=P.fragariae
    echo "$Organism - $Strain"
    for FileF in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_rna/novogene/P.fragariae/Bc1/mycelium/F/*_trim.fq.gz)
    do
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
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
        qsub $ProgDir/sub_star_TA.sh $Assembly $FileF $FileR $OutDir
    done
done

#NOV-9
for Assembly in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/polished_contigs_unmasked.fa)
do
    Strain=Bc16
    Organism=P.fragariae
    echo "$Organism - $Strain"
    for FileF in $(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/qc_rna/novogene/P.fragariae/Nov9/mycelium/F/*_trim.fq.gz)
    do
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
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
        qsub $ProgDir/sub_star_TA.sh $Assembly $FileF $FileR $OutDir
    done
done
```

##Align compressed files of unmapped reads from aligning to vesca

This star script had the following options added to the sub_star.sh script in the ProgDir specified in the below commands:
--winAnchorMultimapNmax 200
--seedSearchStartLmax 30

```bash
#BC-1
for AlignDir in $(ls -d /home/groups/harrisonlab/project_files/phytophthora_fragariae/alignment/star/vesca_alignment/set2/48hr/*)
do
    Organism=P.fragariae
    Strain=Bc16
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
    Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
    echo $Assembly
    OutDir=alignment/star/P.fragariae/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    qsub $ProgDir/sub_star_TA.sh $Assembly $File1 $File2 $OutDir
done

#NOV-9
for AlignDir in $(ls -d /home/groups/harrisonlab/project_files/phytophthora_fragariae/alignment/star/vesca_alignment/set2/72hr/*)
do
    Organism=P.fragariae
    Strain=Bc16
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
    Assembly=$(ls /home/groups/harrisonlab/project_files/phytophthora_fragariae/repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
    echo $Assembly
    OutDir=alignment/star/P.fragariae/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/adamst/git_repos/scripts/popgen/rnaseq
    qsub $ProgDir/sub_star_TA.sh $Assembly $File1 $File2 $OutDir
done
```

#Quantification of gene models

```bash
for Strain in Bc16
do
    for BamFile in $(ls alignment/star/P.fragariae/$Strain/*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v "TA-" | grep -v "countData" | grep -v "genes")
    do
        Gff=gene_pred/annotation/P.fragariae/$Strain/*_genes_incl_ORFeffectors.gff3
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
done
```

##A file was created with columns referring to experimental treatments

```bash
#BC-1
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_Bc1
mkdir -p $OutDir
printf "Sample.name\tTimepoint\tIsolate\n" > $OutDir/P.frag_Bc1_RNAseq_design.txt
# for File in $(ls alignment/star/P.cactorum/10300/Sample_*/Sample_*_featurecounts.txt); do
# Sample=$(echo $File | rev | cut -f2 -d '/' | rev)
# i=$(echo $Sample | sed 's/Sample_//g')
for i in $(seq 1 6)
do
    if [ $i == '1' ] || [ $i == '2' ] || [ $i == '3' ]
    then
        Timepoint='48hr'
    elif [ $i == '4' ] || [ $i == '5' ] || [ $i == '6' ]
    then
        Timepoint='mycelium'
    fi
    Infection=Bc1
    if [ $i == '1' ]
    then
        printf "TA_B_P1\t$Timepoint\t$Infection\n"
    elif [ $i == '2' ]
    then
        printf "TA_B_P2\t$Timepoint\t$Infection\n"
    elif [ $i == '3' ]
    then
        printf "TA_B_P3\t$Timepoint\t$Infection\n"
    elif [ $i == '4' ]
    then
        printf "TA_B_M1\t$Timepoint\t$Infection\n"
    elif [ $i == '5' ]
    then
        printf "TA_B_M2\t$Timepoint\t$Infection\n"
    elif [ $i == '6' ]
    then
        printf "TA_B_M3\t$Timepoint\t$Infection\n"
    fi
done >> $OutDir/P.frag_Bc1_RNAseq_design.txt

# Edit header lines of feature counts files to ensure they have the treatment name rather than file name
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_Bc1
mkdir -p $OutDir
for File in $(ls alignment/star/P.fragariae/Bc16/DeSeq_Bc1/*/*_featurecounts.txt | grep -v 'TA-')
do
    echo $File
    cp $File $OutDir/.
done
for File in $(ls $OutDir/*_featurecounts.txt)
do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g')
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done

#NOV-9
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_Nov9
mkdir -p $OutDir
printf "Sample.name\tTimepoint\tIsolate\n" > $OutDir/P.frag_Nov9_RNAseq_design.txt
# for File in $(ls alignment/star/P.cactorum/10300/Sample_*/Sample_*_featurecounts.txt); do
# Sample=$(echo $File | rev | cut -f2 -d '/' | rev)
# i=$(echo $Sample | sed 's/Sample_//g')
for i in $(seq 1 6)
do
    if [ $i == '1' ] || [ $i == '2' ] || [ $i == '3' ]
    then
        Timepoint='72hr'
    elif [ $i == '4' ] || [ $i == '5' ] || [ $i == '6' ]
    then
        Timepoint='mycelium'
    fi
        Infection='Nov9'
    if [ $i == '1' ]
    then
        printf "TA_NO_P1\t$Timepoint\t$Infection\n"
    elif [ $i == '2' ]
    then
        printf "TA_NO_P2\t$Timepoint\t$Infection\n"
    elif [ $i == '3' ]
    then
        printf "TA_NO_P3\t$Timepoint\t$Infection\n"
    elif [ $i == '4' ]
    then
        printf "TA_NO_M1\t$Timepoint\t$Infection\n"
    elif [ $i == '5' ]
    then
        printf "TA_NO_M2\t$Timepoint\t$Infection\n"
    elif [ $i == '6' ]
    then
        printf "TA_NO_M5\t$Timepoint\t$Infection\n"
    fi
done >> $OutDir/P.frag_Nov9_RNAseq_design.txt

# Edit header lines of feature counts files to ensure they have the treatment name rather than file name
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_Nov9
mkdir -p $OutDir
for File in $(ls alignment/star/P.fragariae/Bc16/DeSeq_Nov9/*/*_featurecounts.txt | grep -v "TA-")
do
    echo $File
    cp $File $OutDir/.
done
for File in $(ls $OutDir/*_featurecounts.txt)
do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g')
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
done
```

#DeSeq commands

BC-1

```R
#install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.fragariae/Bc16/DeSeq_Bc1","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

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
indexes <- c("TA_B_P1", "TA_B_P2", "TA_B_P3", "TA_B_M1", "TA_B_M2", "TA_B_M3")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/P.frag_Bc1_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
countData <- round(countData,0)

design <- ~Group
#design <- colData$Group

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
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

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/heatmap_vst.pdf", width=12,height=12)
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

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plots

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.fragariae/Bc16/DeSeq_Bc1/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#48hr vs mycelium

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Bc1_48hr","Bc1_mycelium"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_48hr_vs_Bc1_mycelium.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_48hr_vs_Bc1_mycelium_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/Bc1_48hr_vs_Bc1_mycelium_down.txt",sep="\t",na="",quote=F)

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Bc1/fpkm_counts.txt",sep="\t",na="",quote=F)
```

NOV-9

```R
#install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

#install and load libraries
require("pheatmap")
require("data.table")

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.fragariae/Bc16/DeSeq_Nov9","*featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treament you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
indexes <- c("TA_NO_P1", "TA_NO_P2", "TA_NO_P3", "TA_NO_M1", "TA_NO_M2", "TA_NO_M5")

countData <- round(countData,0)

#output countData
write.table(countData,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]

#Running DeSeq2

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require("DESeq2")

unorderedColData <- read.table("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/P.frag_Nov9_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
countData <- round(countData,0)

design <- ~Group
#design <- colData$Group

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
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

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/heatmap_vst.pdf", width=12,height=12)
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

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)

#PCA plots

#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.fragariae/Bc16/DeSeq_Nov9/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)

#Analysis of gene expression

#72hr vs mycelium

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Nov9_72hr","Nov9_mycelium"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_72hr_vs_Nov9_mycelium.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_72hr_vs_Nov9_mycelium_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/Nov9_72hr_vs_Nov9_mycelium_down.txt",sep="\t",na="",quote=F)

#Make a table of raw counts, normalised counts and fpkm values:

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.fragariae/Bc16/DeSeq_Nov9/fpkm_counts.txt",sep="\t",na="",quote=F)
```

#Inital analysis of tables of DEGs

```bash
for Strain in Bc1 Nov9
do
    RxLR_Names=analysis/RxLR_effectors/combined_evidence/P.fragariae/Bc16/Bc16_Total_RxLR_motif_hmm.txt
    CRN_Names=analysis/CRN_effectors/hmmer_CRN/P.fragariae/Bc16/Bc16_final_CRN.txt
    ApoP_Names=analysis/ApoplastP/P.fragariae/Bc16/Bc16_Total_ApoplastP.txt
    for File in $(ls alignment/star/P.fragariae/Bc16/DeSeq_"$Strain"/"$Strain"*.txt | grep -v "genes" | grep -v "countData")
    do
        Assessment=$(basename $File | sed "s/.txt//g")
        echo $Assessment
        echo "Total number of genes in dataset:"
        cat $File | grep -v 'baseMean' | wc -l
        echo "Total number of RxLRs in dataset:"
        RxLR_File=$(echo $File | sed "s/.txt/_RxLRs.txt/g")
        cat $File | head -n 1 > $RxLR_File
        cat $File | grep -w -f $RxLR_Names >> $RxLR_File
        cat $RxLR_File | tail -n +2 | wc -l
        echo "Total number of CRNs in dataset:"
        CRN_File=$(echo $File | sed "s/.txt/_CRNs.txt/g")
        cat $File | head -n 1 > $CRN_File
        cat $File | grep -w -f $CRN_Names >> $CRN_File
        cat $CRN_File | tail -n +2 | wc -l
        echo "Total number of Apoplastic effectors in dataset:"
        ApoP_File=$(echo $File | sed "s/.txt/_ApoP.txt/g")
        cat $File | head -n 1 > $ApoP_File
        cat $File | grep -w -f $ApoP_Names >> $ApoP_File
        cat $ApoP_File | tail -n +2 | wc -l
    done
done
```

```
```

##Produce a more detailed table of analyses

```bash
for Strain in Bc16
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
        DEG_Files=$(ls alignment/star/P.fragariae/$Strain/DeSeq_*/*_vs_*.txt  | grep -v -e 'up' -e 'down' -e "CRN" -e "RxLR" | sed -e "s/$/ /g" | tr -d "\n")
        # $ProgDir/pacbio_anntoation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP4 $SigP4 --phobius $PhobiusTxt --RxLR_motif $RxLR_Motif --RxLR_Hmm $RxLR_Hmm --RxLR_WY $RxLR_WY --RxLR_total $RxLR_total --CRN_LFLAK $CRN_LFLAK --CRN_DWL $CRN_DWL --CRN_total $CRN_total --DEG_files $DEG_Files  > $OutDir/414_v2_gene_table_incl_exp.tsv
        # NormCount=$(ls alignment/star/P.cactorum/414_v2/DeSeq/normalised_counts.txt)
        RawCount=$(ls alignment/star/P.fragariae/Bc16/DeSeq_*/raw_counts.txt)
        FPKM=$(ls alignment/star/P.fragariae/Bc16/DeSeq_*/fpkm_counts.txt)
        OrthoName=$Strain
        Assessment=$(echo $RawCount | rev | cut -f2 -d "/" | rev | sed 's/.*_//g')
        OrthoFile=analysis/orthology/orthomcl/All_Strains_plus_rubi_no_removal/All_Strains_plus_rubi_no_removal_orthogroups.txt
        $ProgDir/pacbio_anntoation_tables_modified.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP2_ORF $SigP2_ORF --SigP3 $SigP3 --SigP3_ORF $SigP3_ORF --SigP4 $SigP4 --SigP4_ORF $SigP4_ORF --phobius $PhobiusFa --phobius_ORF $PhobiusFa_ORF --trans_mem $TMHMM_headers --GPI_anchor $GPI_headers --RxLR_total $RxLR_total --RxLR_total_ORF $RxLR_ORF_total --RxLR_EER_total $RxLR_EER_total --RxLR_EER_total_ORF $RxLR_EER_ORF_total --CRN_total $CRN_total --ApoP_total $ApoP_total --ortho_name $OrthoName --ortho_file $OrthoFile --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --Swissprot $SwissProt --InterPro $InterPro > $OutDir/"$Strain"_gene_table_incl_exp_"$Assessment".tsv
    done
done
```

#Draw venn diagrams of differenitally expressed genes

##All genes

###All DEGs

```bash
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
inp1=alignment/star/P.fragariae/Bc16/DeSeq/Bc16_all_DEGs_names.txt
inp2=alignment/star/P.fragariae/Bc16/DeSeq_Bc1/*_vs_*mycelium.txt
inp3=alignment/star/P.fragarie/Bc16/DeSeq_Nov9/*_vs_*mycelium.txt
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2/method_2_all_DEGs.tsv
$ProgDir/parse_RNA-Seq_2.py --input_1 $inp1 --input_2 $inp2 --input_3 $inp3 --out_dir $OutDir
```

###Upregulated DEGs

```bash
cat alignment/star/P.fragariae/Bc16/DeSeq/Bc16_up_DEGs.tsv | cut -f1 > alignment/star/P.fragariae/Bc16/DeSeq/Bc16_up_DEGs_names.txt
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
inp1=alignment/star/P.fragariae/Bc16/DeSeq/Bc16_up_DEGs_names.txt
inp2=alignment/star/P.fragariae/Bc16/DeSeq_Bc1/*_vs_*mycelium_up.txt
inp3=alignment/star/P.fragariae/Bc16/DeSeq_Nov9/*_vs_*mycelium_up.txt
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2/method_2_up_DEGs.tsv
$ProgDir/parse_RNA-Seq_2.py --input_1 $inp1 --input_2 $inp2 --input_3 $inp3 --out_dir $OutDir
```

###Downregulated DEGs

```bash
cat alignment/star/P.fragariae/Bc16/DeSeq/Bc16_down_DEGs.tsv | cut -f1 > alignment/star/P.fragariae/Bc16/DeSeq/Bc16_down_DEGs_names.txt
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
inp1=alignment/star/P.fragariae/Bc16/DeSeq/Bc16_down_DEGs_names.txt
inp2=alignment/star/P.fragariae/Bc16/DeSeq_Bc1/*_vs_*mycelium_down.txt
inp3=alignment/star/P.fragariae/Bc16/DeSeq_Nov9/*_vs_*mycelium_down.txt
OutDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2/method_2_down_DEGs.tsv
$ProgDir/parse_RNA-Seq_2.py --input_1 $inp1 --input_2 $inp2 --input_3 $inp3 --out_dir $OutDir
```

###Venn diagrams

```bash
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
WorkDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_all_DEGs.tsv --out $WorkDir/method_2_all_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_up_DEGs.tsv --out $WorkDir/method_2_up_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_down_DEGs.tsv --out $WorkDir/method_2_down_DEGs.pdf
```

##RxLRs

###Create cut down tables of DEGs just containing RxLRs

```bash
WorkDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2
RxLRs=analysis/RxLR_effectors/combined_evidence/P.fragarie/Bc16/Bc16_Total_RxLR_motif_hmm.txt
cat $WorkDir/method_2_all_DEGs.tsv | grep -w -f $RxLRs > $WorkDir/method_2_all_DEGs_RxLRs.tsv
cat $WorkDir/method_2_up_DEGs.tsv | grep -w -f $RxLRs > $WorkDir/method_2_up_DEGs_RxLRs.tsv
cat $WorkDir/method_2_down_DEGs.tsv | grep -w -f $RxLRs > $WorkDir/method_2_down_DEGs_RxLRs.tsv
```

###Venn diagrams

```bash
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
WorkDir=alignment/star/P.fragariae/$Strain/DeSeq_method_2
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_all_RxLRs_DEGs.tsv --out $WorkDir/method_2_all_RxLRs_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_up_RxLRs_DEGs.tsv --out $WorkDir/method_2_up_RxLRs_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_down_RxLRs_DEGs.tsv --out $WorkDir/method_2_down_RxLRs_DEGs.pdf
```

##CRNs

###Create cut down tables of DEGs just containing CRNs

```bash
WorkDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2
CNRs=analysis/CRN_effectors/hmmer_CRN/P.fragarie/Bc16/Bc16_final_CRN.txt
cat $WorkDir/method_2_all_DEGs.tsv | grep -w -f $CRNs > $WorkDir/method_2_all_DEGs_CRNs.tsv
cat $WorkDir/method_2_up_DEGs.tsv | grep -w -f $CRNs > $WorkDir/method_2_up_DEGs_CRNs.tsv
cat $WorkDir/method_2_down_DEGs.tsv | grep -w -f $CRNs > $WorkDir/method_2_down_DEGs_CRNs.tsv
```

###Venn diagrams

```bash
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
WorkDir=alignment/star/P.fragariae/$Strain/DeSeq_method_2
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_all_CRNs_DEGs.tsv --out $WorkDir/method_2_all_CRNs_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_up_CRNs_DEGs.tsv --out $WorkDir/method_2_up_CRNs_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_down_CRNs_DEGs.tsv --out $WorkDir/method_2_down_CRNs_DEGs.pdf
```

##Apoplastic effectors

###Create cut down tables of DEGs just containing CRNs

```bash
WorkDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2
ApoP=analysis/ApoplastP/P.fragariae/Bc16/Bc16_Total_ApoplastP.txt
cat $WorkDir/method_2_all_DEGs.tsv | grep -w -f $ApoP > $WorkDir/method_2_all_DEGs_ApoP.tsv
cat $WorkDir/method_2_up_DEGs.tsv | grep -w -f $ApoP > $WorkDir/method_2_up_DEGs_ApoP.tsv
cat $WorkDir/method_2_down_DEGs.tsv | grep -w -f $ApoP > $WorkDir/method_2_down_DEGs_ApoP.tsv
```

###Venn diagrams

```bash
ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
WorkDir=alignment/star/P.fragariae/Bc16/DeSeq_method_2
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_all_ApoP_DEGs.tsv --out $WorkDir/method_2_all_ApoP_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_up_ApoP_DEGs.tsv --out $WorkDir/method_2_up_ApoP_DEGs.pdf
$ProgDir/BC-16_method_2_All_DEGs_venn_diag.r --inp $WorkDir/method_2_down_ApoP_DEGs.tsv --out $WorkDir/method_2_down_ApoP_DEGs.pdf
```
