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