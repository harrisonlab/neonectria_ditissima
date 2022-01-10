# Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash
    for Strain in Hg199 R0905_all; do
        F_Read=$(ls qc_dna/paired/N.*/$Strain/F/*trim.fq.gz)
        R_Read=$(ls qc_dna/paired/N.*/$Strain/R/*trim.fq.gz)
        echo $F_Read
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        Outdir=assembly_VP/SPAdes/N.ditissima/$Strain
        sbatch $ProgDir/SPAdes.sh $F_Read $R_Read $Outdir correct
    done
```

```bash
    for Longreads in qc_dna/minion/N.ditissima/Hg199_2ndround/Hg199_fastq_allfiles_trim.fastq.gz; do
        F_Read=$(ls qc_dna/paired/N.*/Hg199/F/*trim.fq.gz)
        R_Read=$(ls qc_dna/paired/N.*/Hg199/R/*trim.fq.gz)
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        SeqType=nanopore
        Outdir=assembly_VP/hybridSPAdes/N.ditissima/Hg199
        sbatch $ProgDir/hybridSPAdes.sh $Longreads $F_Read $R_Read $SeqType $Outdir
    done

    for Longreads in raw_dna/pacbio/N.ditissima/R0905/extracted/concatenated_pacbio.fastq; do
        F_Read=$(ls qc_dna/paired/N.*/R0905_all/F/*trim.fq.gz)
        R_Read=$(ls qc_dna/paired/N.*/R0905_all/R/*trim.fq.gz)
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        SeqType=pacbio
        Outdir=assembly_VP/hybridSPAdes/N.ditissima/R0905
        sbatch $ProgDir/hybridSPAdes.sh $Longreads $F_Read $R_Read $SeqType $Outdir
    done
```

```bash
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/contigs.fasta); do
OutDir=$(dirname $Assembly)/contigs
sbatch $ProgDir/quast.sh $Assembly $OutDir
done

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/scaffolds.fasta); do
OutDir=$(dirname $Assembly)/scaffolds
sbatch $ProgDir/quast.sh $Assembly $OutDir
done
```

## Rename contigs

```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    # If split or remove contigs is needed, provide FCSreport file by NCBI.
    touch tmp.txt
    for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/scaffolds.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/Hg199_hybrid_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
```

    ```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/Hg199_hybrid_renamed.fasta); do
    OutDir=$(dirname $Assembly)/renamed
    sbatch $ProgDir/quast.sh $Assembly $OutDir
    done
    ```

```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    # If split or remove contigs is needed, provide FCSreport file by NCBI.
    touch tmp.txt
    for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/scaffolds.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_hybrid_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
```

    ```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/R0905_hybrid_renamed.fasta); do
    OutDir=$(dirname $Assembly)/renamed
    sbatch $ProgDir/quast.sh $Assembly $OutDir
    done
    ```

# Last seq

## Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
    cd /data/scratch/gomeza/raw_dna
    Species=N.ditissima
    mkdir -p paired/$Species/118923/F
    mkdir -p paired/$Species/118923/R
    mkdir -p paired/$Species/118924/F
    mkdir -p paired/$Species/118924/R
    mkdir -p paired/$Species/226-31/F
    mkdir -p paired/$Species/226-31/R
    mkdir -p paired/$Species/227-31/F
    mkdir -p paired/$Species/227-31/R

    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/118923_S2_L001_R1_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/118923_S2_L001_R1_001.fastq.gz > paired/$Species/118923/F/118923_S2_L001_R1_001.allfiles.fastq.gz
    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/118923_S2_L001_R2_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/118923_S2_L001_R2_001.fastq.gz > paired/$Species/118923/R/118923_S2_L001_R2_001.allfiles.fastq.gz

    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/118924_S1_L001_R1_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/118924_S1_L001_R1_001.fastq.gz > paired/$Species/118924/F/118924_S1_L001_R1_001.allfiles.fastq.gz
    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/118924_S1_L001_R2_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/118924_S1_L001_R2_001.fastq.gz > paired/$Species/118924/R/118924_S1_L001_R2_001.allfiles.fastq.gz

    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/226-31_S3_L001_R1_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/226-31_S3_L001_R1_001.fastq.gz > paired/$Species/226-31/F/226-31_S3_L001_R1_001.allfiles.fastq.gz
    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/226-31_S3_L001_R2_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/226-31_S3_L001_R2_001.fastq.gz > paired/$Species/226-31/R/226-31_S3_L001_R2_001.allfiles.fastq.gz

    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/227-31_S4_L001_R1_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/227-31_S4_L001_R1_001.fastq.gz > paired/$Species/227-31/F/227-31_S4_L001_R1_001.allfiles.fastq.gz
    cat /archives/2021_camb_miseq/RAW/210409_M04543_0001_000000000-DBMRD/Data/Intensities/BaseCalls/227-31_S4_L001_R2_001.fastq.gz /archives/2021_camb_miseq/RAW/210412_M04543_0002_000000000-JDCGC/Data/Intensities/BaseCalls/227-31_S4_L001_R2_001.fastq.gz > paired/$Species/227-31/R/227-31_S4_L001_R2.allfiles.fastq.gz
```


```bash
# Estimate coverage
    for DataDir in $(ls -d paired/N.ditissima/*); do
        F_Read=$(ls $DataDir/F/*.allfiles.fastq.gz)
        R_Read=$(ls $DataDir/R/*.allfiles.fastq.gz)
        echo $F_Read
        echo $R_Read
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch $ProgDir/count_nucl.sh $F_Read $R_Read 45 $DataDir
    done
```

```bash
# Run fastqc
    for Strain in 118923 118924 226-31 227-31; do
    RawData=$(ls raw_dna/paired/*/$Strain/F/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
    done

    for Strain in 118923 118924 226-31 227-31; do
    RawData=$(ls raw_dna/paired/*/$Strain/R/*.fastq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
    done
```

```bash
# Run fastq-mcf
for Strain in 118923 118924 226-31 227-31; do
    Read_F=raw_dna/paired/*/$Strain/F/*.fastq.gz
    Read_R=raw_dna/paired/*/$Strain/R/*.fastq.gz
    echo $Read_F;
    echo $Read_R;
    IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

```bash
# Run fastqc
    for Strain in 118923 118924 226-31 227-31; do
    RawData=$(ls qc_dna/paired/*/$Strain/F/*.fq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
    done

    for Strain in 118923 118924 226-31 227-31; do
    RawData=$(ls qc_dna/paired/*/$Strain/R/*.fq.gz)
    echo $RawData;
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
    sbatch $ProgDir/fastqc.sh $RawData
    done
```

## SPAdes

```bash
    for Strain in 118923 118924 226-31 227-31; do
        F_Read=$(ls ../../data/scratch/gomeza/qc_dna/paired/N.*/$Strain/F/*trim.fq.gz)
        R_Read=$(ls ../../data/scratch/gomeza/qc_dna/paired/N.*/$Strain/R/*trim.fq.gz)
        echo $F_Read
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        Outdir=assembly_VP/SPAdes/N.ditissima/$Strain
        sbatch $ProgDir/SPAdes.sh $F_Read $R_Read $Outdir correct
    done
```

## Rename contigs

Rename the sequences in assembly fasta file to have simple names.

```
Note: Old isolates were assembled before with an older version of SPAdes and masked. These assembies were not renamed before running repeatmasker, so genomes and gene models keep the contig names from SPAdes. If I need to update the annotation (probably not) I should use this renamed file.
```

```bash
    for Strain in 118923 118924 226-31 227-31 Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C2 SVK2 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C3 SVK1 NMaj; do
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    touch tmp.txt
        for Assembly in $(ls assembly_VP/SPAdes/N.*/$Strain/filtered_contigs/contigs_min_500bp.fasta); do
        OutDir=$(dirname $Assembly)
        echo "$Assembly"
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
        done
    rm tmp.txt
    done
```
## Quast

```bash
    for Strain in 118923 118924 226-31 227-31; do
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
        for Assembly in $(ls assembly/SPAdes/N.ditissima/$Strain/filtered_contigs/*renamed.fasta); do
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
        OutDir=$(dirname $Assembly)
        sbatch $ProgDir/quast.sh $Assembly $OutDir
        done
    done
```

## Repeatmask

```bash
    for Strain in 118923 118924 226-31 227-31; do
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
        BestAssembly=assembly_VP/SPAdes/N.ditissima/$Strain/filtered_contigs/*renamed.fasta
        OutDir=repeat_masked_VP/SPAdes_assembly/N.ditissima/$Strain/
        sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
        sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir
    done
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.


```bash
    for Strain in 118923 118924 226-31 227-31; do
        for File in $(ls repeat_masked_VP/*/N.d*/$Strain/*_contigs_softmasked.fa); do
        OutDir=$(dirname $File)
        TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
        OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
        echo "$OutFile"
        bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
        echo "Number of masked bases:"
        cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
        done
    done
```