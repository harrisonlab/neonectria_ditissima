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
cp /home/gomeza/prog/genemark/gm_key_Nov2020/gm_key_64 ~/.gm_key_64
```

```bash
# If running on a conda env
perl change_path_in_perl_scripts.pl "/usr/bin/env perl" 
for Strain in 118923 118924 226-31 227-31; do
for Assembly in $(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=alignment/star/N.ditissima/$Strain/Hg199_reads/concatenated/concatenated.bam 
GeneModelName="$Organism"_"$Strain"_braker 
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
sbatch -p himem $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
done
#rm -r /home/gomeza/miniconda3/envs/gene_pred/config/species/N.ditissima_*
# Braker needs to be fixed in out HPC
```

Data was copied into the Crop diversity HPC. Braker was installed there in the gene_pred env

```bash
# If running on a conda env
perl change_path_in_perl_scripts.pl "/usr/bin/env perl" 
# Scripts and files in CropDiversity
    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls neonectria_ditissima_WD/repeat_masked/SPAdes_assembly/N.ditissima/$Strain/"$Strain"_contigs_softmasked_repeatmasker_TPSI_appended.fa ); do
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        OutDir=neonectria_ditissima_WD/gene_pred/braker/$Organism/$Strain
        AcceptedHits=neonectria_ditissima_WD/alignment/star/N.ditissima/$Strain/Hg199_reads/concatenated/concatenated.bam
        GeneModelName="$Organism"_"$Strain"_braker 
        ProgDir=apps/scripts/Gene_pred
        sbatch -p long $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
        done
    done
```

## Cufflinks RNA-seq alignments assembler

```bash
    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_contigs_unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
        mkdir -p $OutDir
        AcceptedHits=alignment/star/N.ditissima/$Strain/Hg199_reads/concatenated/concatenated.bam
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        sbatch -p himem $ProgDir/cufflinks.sh $AcceptedHits $OutDir
        done
    done
```

## Codingquarry 

```bash
conda activate antismash_py27

    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_contigs_unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        OutDir=gene_pred/codingquarry/$Organism/$Strain
        mkdir -p $OutDir
        GTF=gene_pred/cufflinks/N.ditissima/$Strain/concatenated_prelim/transcripts.gtf
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
    done
```


### Add additional transcripts to Braker gene models.


Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
for Strain in 118923 118924 226-31 227-31; do
    for BrakerGff in $(ls gene_pred/braker/N.ditissima/$Strain/augustus.hints.gff3); do
        Strain=$(echo $BrakerGff| rev | cut -d '/' -f2 | rev)
        Organism=$(echo $BrakerGff | rev | cut -d '/' -f3 | rev)
        echo "$Organism - $Strain"
        Assembly=$(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
        CodingQuarryGff=gene_pred/codingquarry/N.ditissima/$Strain/out/PredictedPass.gff3
        PGNGff=gene_pred/codingquarry/N.ditissima/$Strain/out/PGN_predictedPass.gff3
        AddDir=gene_pred/codingquary/$Organism/$Strain/additional # Additional transcripts directory
        FinalDir=gene_pred/codingquarry/$Organism/$Strain/final # Final directory
        AddGenesList=$AddDir/additional_genes.txt
        AddGenesGff=$AddDir/additional_genes.gff
        FinalGff=$AddDir/combined_genes.gff
        mkdir -p $AddDir
        mkdir -p $FinalDir
        echo "$CodingQuarryGff - $BrakerGff"
        # Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
        bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
        bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList

        # Creat Gff file with the additional transcripts
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
        $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff

        # Create a final Gff file with gene features
        $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuarry.gff3

        # Create fasta files from each gene feature in the CodingQuarry gff3
        $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuarry.gff3 $FinalDir/final_genes_CodingQuarry

        # Create fasta files from each gene feature in the Braker gff3
        cp $BrakerGff $FinalDir/final_genes_Braker.gff3
        $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

        # Combine both fasta files
        cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuarry.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
        cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuarry.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
        cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuarry.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
        cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuarry.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

        # Combine both gff3 files
        GffBraker=$FinalDir/final_genes_CodingQuarry.gff3
        GffQuarry=$FinalDir/final_genes_Braker.gff3
        GffAppended=$FinalDir/final_genes_appended.gff3
        cat $GffBraker $GffQuarry > $GffAppended
    done
done

  # Check the final number of genes
for Strain in 118923 118924 226-31 227-31; do
for DirPath in $(ls -d gene_pred/codingquarry/N.ditissima/$Strain/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuarry.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
  ```

  ### Remove duplicate and rename genes.

  ```bash
    for Strain in 118923 118924 226-31 227-31; do
    for GffAppended in $(ls gene_pred/codingquarry/N.ditissima/$Strain/final/final_genes_appended.gff3);
    do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/codingquarry/$Organism/$Strain/final
    # Remove duplicated genes
    GffFiltered=gene_pred/codingquarry/N.ditissima/$Strain/final/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    # Rename genes
    GffRenamed=gene_pred/codingquarry/N.ditissima/$Strain/final/final_genes_appended_renamed.gff3
    LogFile=gene_pred/codingquarry/N.ditissima/$Strain/final/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    # Create renamed fasta files from each gene feature   
    Assembly=$(ls repeat_masked/SPAdes_assembly/N.ditissima/$Strain/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
    done 
    done
```


### SignalP and tmhmm 

This section needs to run on gruff
```bash
conda activate annotation
```
```bash
    for Strain in 118923 118924 226-31 227-31; do
    ProgDir=/home/agomez/scratch/apps/scripts/Feature_annotation
    CurPath=$PWD
        for Proteome in $(ls gene_pred/$Strain/final/final_genes_appended_renamed.pep.fasta); do
        Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
        Organism=N.ditissima
        SplitDir=gene_pred/final_genes_split/$Organism/$Strain
        mkdir -p $SplitDir
        BaseName="$Organism""_$Strain"_final_preds
        $ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName 
        for File in $(ls $SplitDir/*_final_preds_*); do
        sbatch $ProgDir/pred_signalP.sh $File signalp-4.1 
        done
        done
    done
```


The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

 ```bash
for Strain in 118923 118924 226-31 227-31; do
  for SplitDir in $(ls -d gene_pred/final_genes_split/*/$Strain); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=final_genes_signalp-4.1
    for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
      InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
      InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
      InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
  done
done
```

Proteins containing a transmembrane domain were identified:

```bash
for Strain in 118923 118924 226-31 227-31; do 
for Proteome in $(ls gene_pred/$Strain/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/agomez/scratch/apps/scripts/Feature_annotation
sbatch $ProgDir/TMHMM.sh $Proteome
done
done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
for File in $(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_neg.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
  cat $File | cut -f1 > $TmHeaders
  SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
  OutDir=$(dirname $SigP)
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
  cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
```



echo "Running using SignalP version: $SigP_Version"
  $SigP_Version -t euk -f summary -c 70 "proteins.fa" > "$OutFile"_sp.txt
  echo '----------------------------------------------------------------------' >> "$OutFile"_sp.txt
  PathToAnnotateSigP=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $PathToAnnotateSigP/sigP_4.1_parser.py --inp_sigP "$OutFile"_sp.txt --out_tab "$OutFile"_sp.tab --out_fasta "$OutFile"_sp.aa --out_neg "$OutFile"_sp_neg.aa --inp_fasta "proteins.fa"
  OutDir=$CurPath/gene_pred/"${Source}_${SigP_Version}"/$Organism/$Strain/split

perl -d:Trace /scratch/software/signalp-4.1/signalp -t euk -f summary -c 70 N.ditissima_118923_final_preds_10000.fa sp.txt >trace 2>&1
  perl -d:Trace /data/scratch/gomeza/prog/signalp/signalp-4.1/signalp-4.1 -t euk -f summary -c 70  sp.txt >trace 2>&1
  perl -d:Trace /scratch/software/signalp-4.1/signalp -t euk -f summary -c 70  sp.txt >trace 2>&1
perl -d:Trace /data/scratch/gomeza/prog/signalp/signalp-4.1/signalp-4.1 -t euk -f summary -c 70 N.ditissima_118923_final_preds_10000.fa sp.txt >trace 2>&1

  scp -r test_sigp/ agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch/neonectria_ditissima_WD
scp -r /data/scratch/gomeza/prog/signalp/signalp-4.1 agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch/apps/prog
```
perl -d:Trace signalp-4.1 -t euk -f summary -c 70  sp.txt >trace 2>&1
 sed -i "s&SIGNALP=.*&SIGNALP=$PWD/signalp-4.1&g" signalp-4.1
  sed -i 's&AWK=nawk&AWK=/usr/bin/awk&g' signalp-4.1
  sed -i 's&AWK=/usr/bin/gawk&AWK=/usr/bin/awk&g' signalp-4.1

srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash


/home/gomeza/.local/bin/signalp6


perl -d:Trace signalp-4.1 -t euk -f summary -c 70  sp.txt >trace 2>&1


mamba install pytorch torchvision cudatoolkit=10.2 -c pytorch
pip install signalp-6-package
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/


signalp6 --fastafile N.ditissima_118923_final_preds_10000.fa --organism eukarya --output_dir ./ --format txt --mode fast


../../apps/prog/signalp-4.1/signalp-4.1 -t euk -f summary -c 70 N.ditissima_118923_final_preds_10000.fa sp.txt

../../apps/scripts/Feature_annotation/sigP_4.1_parser.py --inp_sigP sp.txt --out_tab sp.tab --out_fasta sp.aa --out_neg sp_neg.aa --inp_fasta N.ditissima_118923_final_preds_10000.fa