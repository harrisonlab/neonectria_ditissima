Sets up correct formatting for SNP calling analysis

```bash
input=/home/groups/harrisonlab/project_files/neonectria_ditissima/analysis/genome_alignment/bowtie
scripts=/home/gomeza/git_repos/scripts/popgen/snp
```

## Rename input mapping files in each folder by prefixing with the strain ID

```bash

for filename in $(ls -d analysis/genome_alignment/bowtie/*/$Strain/vs_R0905_canu_2017); do
Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
Strain=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
    cp "$filename/R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam" "$filename/$Strain_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam"
done

cd $input/*/Hg199/vs_R0905_canu_2017
for filename in R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam
do
    mv "$filename" "Hg199_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam"
done

cd $input/*/AgN04/vs_R0905_canu_2017
for filename in R0905_contigs_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam
do
    mv "$filename" "AgN04_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam"
done

cd $input/*/R0905/vs_NG-R0905_pacbio
for filename in polished_contigs_unmasked.fa_aligned.sam
do
    mv "$filename" "R0905_softmasked_repeatmasker_TPSI_appended.fa_aligned.sam"
done


## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)
Convention used:
qsub $scripts/sub_pre_snp_calling.sh <INPUT SAM FILE> <SAMPLE_ID>

```bash
for Strain in A4 Bc1 Bc16 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2 SCRP249 SCRP324 SCRP333
do
    Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'sub_pre_sn' | wc -l)
    done
    qsub $scripts/sub_pre_snp_calling.sh $input/*/$Strain/vs_Bc16_FALCON/"$Strain"_polished_contigs_unmasked.fa_aligned.sam $Strain
done
```

##Copy outputs from cleanup to alignment folder

```bash
for Strain in A4 Bc1 Bc16 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2 SCRP249 SCRP324 SCRP333
do
    Bam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
    rgBam="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
    Bai="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
    Txt="$Strain"_polished_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Bc16_FALCON/
    mv $Bam $Directory
    mv $rgBam $Directory
    mv $Bai $Directory
    mv $Txt $Directory
done
```