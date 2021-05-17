# Gene Prediction

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls ../../../projects/neonectria_ditissima/repeat_masked_VP/SPAdes_assembly/N.ditissima/$Strain/*unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO_v5/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
        done
    done

    for Strain in Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C2 SVK2 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C3 SVK1 NMaj; do
        for Assembly in $(ls repeat_masked/N.*/$Strain/*unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch -p long $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
        done
    done

      for Strain in Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C2 SVK2 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C3 SVK1 NMaj; do
        for Assembly in $(ls repeat_masked/N.*/$Strain/*unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch -p long $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
        done
    done

    for Assembly in $(ls repeat_masked/Nz_genomes/*/*_nt.fasta); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
        Organism=N.ditissima
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch -p long $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done
```
This was previously done with Busco v3.0.2 and the sordarymyceta_odb9. Re run to update the phylogeny


## Star

Spliced Transcripts Alignment to a Reference. 


```bash
    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls repeat_masked_VP/SPAdes_assembly/N.ditissima/$Strain/*_contigs_unmasked.fa); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev) 
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        echo $Assembly
        echo "$Organism - $Strain"
            for FileF in $(ls -d qc_rna/RNAseq/N.ditissima/Hg199/mycelium/F/Hg199_3_1_trim.fq.gz) ; do
            FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/1_trim.fq/2_trim.fq/g')
            echo $FileF
            echo $FileR
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
            echo $Sample_Name
            OutDir=alignment_VP/star/$Organism/$Strain/Hg199_reads/$Sample_Name
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
            sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir 11
            done
        done
    done
```

### Alignment outputs were concatenated and Braker1 prediction was run

```bash
    for Strain in 118923 118924 226-31 227-31; do
        mkdir -p alignment_VP/star/N.ditissima/$Strain/Hg199_reads/concatenated
        samtools merge -f alignment_VP/star/N.ditissima/$Strain/Hg199_reads/concatenated/concatenated.bam \
        alignment_VP/star/N.ditissima/$Strain/Hg199_reads/Hg199_1/star_aligmentAligned.sortedByCoord.out.bam \
        alignment_VP/star/N.ditissima/$Strain/Hg199_reads/Hg199_2/star_aligmentAligned.sortedByCoord.out.bam \
        alignment_VP/star/N.ditissima/$Strain/Hg199_reads/Hg199_3/star_aligmentAligned.sortedByCoord.out.bam
    done
```

## Braker prediction

Before braker predictiction was performed, I double checked that I had the genemark key in my user area and copied it over from the genemark install directory:

```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/2018/gm_key_64 ~/.gm_key
```

```bash
  for Assembly in $(ls assembly_vAG/canu_1step/N.ditissima/R0905/polished/repeat_masked/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f6 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred_vAG/braker/Ref_Genomes/$Organism/$Strain
    AcceptedHits=alignment_vAG/star/N.ditissima/Ref_Genomes/R0905/Hg199_reads/mycelium/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    #rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```
