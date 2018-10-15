Initial analysis of RNA-Seq data was performed in RNA-Seq_analysis.md
this performs co-expression analyses and searches for potential promotor motifs.

# Four key stages of analysis to be performed

A) Building a coexpression network and export it
B) Counting kmers in the 3kb upstream region of all genes in the coexpressed set
C) Test these kmers for enrichment

## A) Building a coexpression network

The threshold for samples depends on the dataset,
make sure to sanity check with the output graph
cluster to keep may also need to be modified

All columns MUST have headers

```bash
OutDir=analysis/coexpression
mkdir -p $OutDir
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
gene_table=gene_pred/annotation/Pathogen/N.ditissima/Hg199/Hg199_gene_table_incl_exp_COMPLETE.tsv
FPKM_start=13
FPKM_end=14
qsub $ProgDir/sub_WGCNA_input_clean.sh $gene_table $OutDir $FPKM_start $FPKM_end
```

```
864 genes were removed due to too many missing samples or zero variance.
Listed in analysis/coexpression/removed_genes.txt
```

Test various values of soft thresholding for building the network

```bash
OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
max_SFT=60
qsub $ProgDir/sub_choose_softthreshold.sh $OutDir 60
```

```
scale-free topology fit index values seem to peak and then decline.
WGCNA tutorial recommends getting as close to 0.9 as possible,
but 48 has a value of 0.8040 and is the highest reached before declining.
```

Build the coexpression network

Merging threshold value may need tweaking


```bash
OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
qsub $ProgDir/sub_create_network.sh $OutDir 48 60 0.1
```

```

With a minimum size of 60:
41 modules were created after merging with a threshold of 0.05

With a minimum size of 30:
56 modules were created after merging with a threshold of 0.05

With a minimum size of 60:
74 modules were created after merging with a threshold of 0.1

```

Export gene lists

Y or N denotes whether unmerged modules will be exported too

```bash
OutDir=analysis/coexpression
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
qsub $ProgDir/sub_Export_Gene_Lists.sh $OutDir Y
```

```
Decided on a minimum module size of 60, with a merging threshold of 0.1
```
I made a list of highly expressed effectors and Cazymes of my dataset. I have checked how they are grouped into the modules.

Secreted effectors
```bash
grep -rnw './' -e 'g4785.t1' -e 'g7947.t1' -e 'g8212.t1' -e 'g12278.t1' -e 'g5856.t1' -e 'g3244.t1' -e 'g8388.t1' -e 'g11580.t1' -e 'g8135.t1' -e 'g5859.t1'
```
./Genes_in_blue.txt:369:"g8135.t1"
./Genes_in_darkred.txt:329:"g3244.t1"
./Genes_in_darkred.txt:510:"g5856.t1"
./Genes_in_green.txt:327:"g8135.t1"
./Genes_in_bisque4.txt:38:"g4785.t1"
./Genes_in_bisque4.txt:66:"g8388.t1"
./Genes_in_darkseagreen2.txt:19:"g5859.t1"
./Genes_in_darkolivegreen1.txt:25:"g8212.t1"
./Genes_in_darkgrey.txt:189:"g8212.t1"
./Genes_in_grey60.txt:681:"g3244.t1"
./Genes_in_grey60.txt:1184:"g5856.t1"
./Genes_in_grey60.txt:1603:"g7947.t1"
./Genes_in_darkslateblue.txt:80:"g4785.t1"
./Genes_in_darkslateblue.txt:99:"g5859.t1"
./Genes_in_darkslateblue.txt:131:"g8388.t1"

CAZY_secreted
```bash
grep -rnw './' -e 'g3929.t1' -e 'g12751.t1' -e 'g921.t1' -e 'g8217.t1' -e 'g12079.t1' -e 'g7755.t1' -e 'g11119.t1' -e 'g2773.t1' -e 'g12485.t1' -e 'g4785.t1'
```
./Genes_in_greenyellow.txt:871:"g11119.t1"
./Genes_in_darkturquoise.txt:449:"g12079.t1"
./Genes_in_darkturquoise.txt:481:"g12751.t1"
./Genes_in_bisque4.txt:38:"g4785.t1"
./Genes_in_blue2.txt:185:"g3929.t1"
./Genes_in_pink.txt:165:"g7755.t1"
./Genes_in_darkseagreen2.txt:2:"g921.t1"
./Genes_in_honeydew1.txt:120:"g8217.t1"
./Genes_in_darkviolet.txt:738:"g12751.t1"
./Genes_in_sienna3.txt:1743:"g11119.t1"
./Genes_in_skyblue3.txt:332:"g8217.t1"
./Genes_in_skyblue3.txt:481:"g12079.t1"
./Genes_in_grey60.txt:587:"g2773.t1"
./Genes_in_darkslateblue.txt:15:"g921.t1"
./Genes_in_darkslateblue.txt:80:"g4785.t1"
./Genes_in_darkseagreen3.txt:213:"g3929.t1"
./Genes_in_black.txt:212:"g7755.t1"
./Genes_in_honeydew.txt:169:"g12485.t1"


## Investigate enrichment of effector class genes in co-expression modules

### Create lists of RxLRs, CRNs, ApoplastP hits, a combined list of all effectors and a list of secreted proteins

Remove statement from beginning of removed_genes file from WGCNA
Create files for Genomic sets of Genes and features

```bash
cat analysis/coexpression2/removed_genes.txt | sed 's/Removing genes: //g' > analysis/coexpression/removed_genes.txt
Removed_Genes=analysis/coexpression/removed_genes.txt

# All genes

Gene_Set=analysis/coexpression/Genes_for_WGCNA.txt
Genome_Fasta=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta
Genome_Headers=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined_Headers.txt
cat $Genome_Fasta | grep '>' | sort | uniq > $Genome_Headers
cat $Genome_Headers | grep -wvf $Removed_Genes | sort | uniq > $Gene_Set

sed -i 's/>//g' Genes_for_WGCNA.txt

# CAZY

CAZY_Set=analysis/coexpression/CAZY_analysed.txt
CAZY_Headers=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
cat $CAZY_Headers | grep -wf $Gene_Set | sort | uniq > $CAZY_Set

# Effector

Effector_Set=analysis/coexpression/Effector_analysed.txt
Effector_Headers=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
cat $Effector_Headers | grep -wf $Gene_Set | sort | uniq > $Effector_Set

# Effectors_Combined

Total_Effector=analysis/coexpression/Total_Effector_Set.txt
cat $CAZY_Set $Effector_Set | sort | uniq > $Total_Effector

# Secreted proteins

Secreted_Set=analysis/coexpression/Secreted_analysed.txt
Augustus_Secreted=gene_pred/final_genes_signalp-4.1/N.ditissima/Hg199_minion/Hg199_minion_final_sp_no_trans_mem_headers.txt
cat $Augustus_Secreted | grep -wf $Gene_Set | sort | uniq > $Secreted_Set
```

Create files for Module sets of Genes and features

```bash
for File in $(ls analysis/coexpression/merged_modules/Genes_in_*.txt)
do
    Filename=$(echo $File | rev | cut -f1 -d "/" | rev)
    Module_ID=$(echo $Filename | cut -f3 -d "_" | cut -f1 -d ".")
    echo "Processing $Module_ID"
    Module_Dir=analysis/coexpression/enrichment/$Module_ID
    mkdir -p $Module_Dir
    Gene_Set_unsorted=$Module_Dir/Gene_Set_unsorted.txt
    Gene_Set=$Module_Dir/Gene_Set.txt
    # Copy entire gene set of module
    cp $File $Gene_Set_unsorted
    cat $Gene_Set_unsorted | sed 's/"//g' | sort | uniq > $Gene_Set
    rm $Gene_Set_unsorted
    echo "Gene set file created for $Module_ID"
    # Create CAZY list
    CAZY_Headers=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
    CAZymes_In_Module=$Module_Dir/CAZY_IDs.txt
    cat $CAZY_Headers | grep -wf $Gene_Set | sort | uniq > $CAZymes_In_Module
    echo "CAZymes file created for $Module_ID"
    # Create list of Effector hits
    Effector_Headers=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
    Effector_In_Module=$Module_Dir/Effector_IDs.txt
    cat $Effector_Headers | grep -wf $Gene_Set | sort | uniq > $Effector_In_Module
    echo "Effector hits file created for $Module_ID"
    # Create combined effector list
    Effectors_Combined_In_Module=$Module_Dir/Effectors_Combined_IDs.txt
    cat $Effector_In_Module $Cazymes_In_Module | sort | uniq > $Effectors_Combined_In_Module
    echo "Effector combined list created for $Module_ID"
    # Create combined secreted proteins list
    Secreted_In_Module=$Module_Dir/Secreted_IDs.txt
    Augustus_Secreted=gene_pred/final_genes_signalp-4.1/N.ditissima/Hg199_minion/Hg199_minion_final_sp_no_trans_mem_headers.txt
    cat $Augustus_Secreted | grep -wf $Gene_Set | sort | uniq > $Secreted_In_Module
    echo "Secreted list created for $Module_ID"
done
```

### Create contigency tables for fishers exact test

```bash
Effector_Headers=analysis/effectorP/N.ditissima/Hg199_minion/N.ditissima_Hg199_minion_EffectorP_headers.txt
Genome_Effector=$(cat $Effector_Headers | sort | uniq | wc -l)
CAZY_Headers=gene_pred/CAZY/N.ditissima/Hg199_minion/Hg199_minion_CAZY_headers.txt
Genome_CAZY=$(cat $CAZY_Headers | sort | uniq | wc -l)
Genome_Effectors_Combined=$(cat $CAZY_Headers $Effector_Headers | sort | uniq | wc -l)
Augustus_Secreted=gene_pred/final_genes_signalp-4.1/N.ditissima/Hg199_minion/Hg199_minion_final_sp_no_trans_mem_headers.txt
Genome_Secreted=$(cat $Augustus_Secreted | sort | uniq | wc -l)
Genes=gene_pred/codingquary/N.ditissima/Hg199_minion/final/final_genes_combined.cdna.fasta
Genome_Genes=$(cat $Genes | sort | uniq | wc -l)

for Module_Dir in $(ls -d analysis/coexpression/enrichment/*)
do
    Module_Name=$(echo $Module_Dir | rev | cut -f1 -d "/" | rev)
    Module_Effector=$Module_Dir/Effector_IDs.txt
    Module_CAZY=$Module_Dir/CAZY_IDs.txt
    Module_Effectors_Combined=$Module_Dir/Effectors_Combined_IDs.txt
    Module_Secreted=$Module_Dir/Secreted_IDs.txt
    Module_Genes=$Module_Dir/Gene_Set.txt
    ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
    echo "Processing $Module_Name"
    python $ProgDir/Fisher_table_prep.py --Module_Effector $Module_Effector --Module_CAZY $Module_CAZY  --Module_Effectors_Combined $Module_Effectors_Combined --Module_Secreted $Module_Secreted --Module_Genes $Module_Genes --Module_Name $Module_Name --Genome_Effector $Genome_Effector --Genome_CAZY $Genome_CAZY --Genome_Effectors_Combined $Genome_Effectors_Combined --Genome_Secreted $Genome_Secreted --Genome_Genes $Genome_Genes --OutDir $Module_Dir
    echo "$Module_Name done"
done
```

### Run Fishers Exact test to test for enrichment

Does not run on head node, throws an error of incorrect version of GLOB
I recommend running in screen due to large numbers of items in loop

```bash
qlogin
cd /data/scratch/gomeza

for Table in $(ls analysis/coexpression/enrichment/*/*Fishertable.txt)
do
    Module_ID=$(echo $Table | rev | cut -f2 -d "/" | rev)
    Filename=$(echo $Table | rev | cut -f1 -d "/" | rev)
    Gene_Type=$(echo $Filename | cut -f2 -d "_")
    OutDir=analysis/coexpression/enrichment/$Module_ID/Fisher_Results
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
    echo "Running $Module_ID" "$Gene_Type"
    /home/adamst/prog/R/R-3.2.5/bin/Rscript --vanilla $ProgDir/fisherstest.R --Input_Table $Table --Output_Directory $OutDir --Module_ID $Module_ID --Gene_Type $Gene_Type
    echo "$Module_ID" "$Gene_Type done"
done
```

Parse fisher results files to fewer files by type and hypothesis being tested

```bash
for Type in up down equal
do
    Files=$(ls analysis/coexpression/enrichment/*/Fisher_Results/enriched_"$Type".txt)
    OutDir=analysis/coexpression/enrichment/Parsed_Fisher_Results/$Type
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/RNAseq_analysis/Co-Expression_Analysis
    python $ProgDir/parse_fisher_results.py --inputs $Files --outdir $OutDir --FDR 0.05 --Types Effector CAZY Effectors_Combined Secreted --Threshold 0.05
done
```

None are listed as equal, so this execution of the script fails
Summary of results

```
RxLRs:
Significantly enriched: 32
Non-Significantly enriched: 61
Significantly lacking: 0
Non-Significantly lacking: 96

CRNs:
Significantly enriched: 5
Non-Significantly enriched: 6
Significantly lacking: 0
Non-Significantly lacking: 158

ApoP:
Significantly enriched: 31
Non-Significantly enriched: 98
Significantly lacking: 0
Non-Significantly lacking: 60

Effectors:
Significantly enriched: 52
Non-Significantly enriched: 90
Significantly lacking: 0
Non-Significantly lacking: 47

Secreted
Significantly enriched: 159
Non-Significantly enriched: 28
Significantly lacking: 0
Non-Significantly lacking: 2
```

Repeat this analysis for a minimum module size of 30 to see if this is more appropriate

Genes from the steelblue module were visually inspected for promotor hunting
This gave 12 high confidence genes and 48 lower confidence genes
Also run the entirety of steelblue, 212 genes so might work better

Extract the 3000 bases upstream of these genes

```bash
for Headers in $(ls promotor_id/large_modules/*genes.txt)
do
    Sequences=gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.upstream3000.fasta
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    OutDir=promotor_id/large_modules
    File_ID=$(echo $Headers | cut -f3 -d '/' | cut -f1 -d '_')
    $ProgDir/extract_from_fasta.py --fasta $Sequences --headers $Headers > $OutDir/"$File_ID"_upstream3000.fasta
done
```

Convert all bases to uppercase

```bash
for input in $(ls promotor_id/large_modules/*.fasta)
do
    filename=$(echo $input | cut -f1 -d '.')
    input_modified="$filename"_modified.fasta
    cat $input | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' \
    > $input_modified
    rm $input
    mv $input_modified $input
done
```

Remove duplicate genes

```bash
for input in $(ls promotor_id/large_modules/*.fasta)
do
    filename=$(echo $input | cut -f1 -d '.')
    input_modified="$filename"_modified.fasta
    awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' $input | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > $input_modified
    rm $input
    mv $input_modified $input
done
```

## B) Identifying novel motifs using DREME

### Set up files

Organise directories and create files for a comparison set to sample from

```bash
for Set in all highconfidence totalmodule
do
    mkdir -p promotor_id/large_modules/$Set
    mv promotor_id/large_modules/"$Set"_genes.txt promotor_id/large_modules/$Set/.
    mv promotor_id/large_modules/"$Set"_upstream3000.fasta promotor_id/large_modules/$Set/.
    cat promotor_id/large_modules/Total_Gene_Set.txt | grep -vf promotor_id/large_modules/$Set/"$Set"_genes.txt > promotor_id/large_modules/$Set/"$Set"_nontarget.txt
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    Fasta=gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.upstream3000.fasta
    Headers=promotor_id/large_modules/$Set/"$Set"_nontarget.txt
    Output=promotor_id/large_modules/$Set/"$Set"_nontarget_upstream3000.fasta
    echo "Extracting fasta for $Set"
    $ProgDir/extract_from_fasta.py --fasta $Fasta --headers $Headers > $Output
    echo "$Set done"
done
```

Split fasta file into 100bp sequences

Before running, install pyfasta

```bash
pip install https://pypi.python.org/packages/be/3f/794fbcdaaa2113f0a1d16a962463896c1a6bdab77bd63f33a8f16aae6cdc/pyfasta-0.5.2.tar.gz --user
```

To use this module, add the following line to your profile

```bash
export PYTHONPATH="$PYTHONPATH:/home/adamst/.local/lib/python2.7/site-packages"
export PATH=${PATH}:/home/adamst/.local/bin
```

Now split target files into 100bp sequences

```bash
for Set in all highconfidence totalmodule
do
    WorkDir=promotor_id/large_modules/$Set
    Fasta=$WorkDir/"$Set"_upstream3000.fasta
    pyfasta split -n 1 -k 100 $Fasta
done
```

Convert all bases to uppercase in non-target files

```bash
for Set in all highconfidence totalmodule
do
    for input in $(ls promotor_id/large_modules/$Set/"$Set"_nontarget_upstream3000.fasta)
    do
        filename=$(echo $input | cut -f1 -d '.')
        input_modified="$filename"_modified.fasta
        cat $input | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' \
        > $input_modified
        rm $input
        mv $input_modified $input
    done
done
```

Remove duplicate genes in non-target files

```bash
for Set in all highconfidence totalmodule
do
    for input in $(ls promotor_id/large_modules/$Set/"$Set"_nontarget_upstream3000.fasta)
    do
        filename=$(echo $input | cut -f1 -d '.')
        input_modified="$filename"_modified.fasta
        awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' $input | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > $input_modified
        rm $input
        mv $input_modified $input
    done
done
```

Now split non-target files into 100bp sequences

```bash
for Set in all highconfidence totalmodule
do
    WorkDir=promotor_id/large_modules/$Set
    Fasta=$WorkDir/"$Set"_nontarget_upstream3000.fasta
    pyfasta split -n 1 -k 100 $Fasta
done
```

### Create negative samples for each set

```bash
for Set in all highconfidence totalmodule
do
    for Rep in {1..100}
    do
        WorkDir=promotor_id/large_modules/$Set
        Pos_Fasta=$WorkDir/"$Set"_upstream3000.split.100mer.fasta
        Neg_Fasta=$WorkDir/"$Set"_nontarget_upstream3000.split.100mer.fasta
        Num_of_Seqs=$(cat $Pos_Fasta | grep '>' | wc -l)
        ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
        OutDir=/data/scratch/adamst/large_modules/$Set
        mkdir -p $OutDir
        Jobs=$(qstat | grep 'sub_fasta' | wc -l)
        while [ $Jobs -gt 20 ]
        do
            sleep 1
            printf "."
            Jobs=$(qstat | grep 'sub_fasta' | wc -l)
        done
        qsub $ProgDir/sub_fasta_subsample.sh $Neg_Fasta $Num_of_Seqs $Rep $OutDir
    done
done
```

### DREME motif analysis

fifth command line argument is e-value, default and recommended is 0.05
Increasing will provide more motifs, but they may not be significantly enriched

```bash
for Rep in {1..100}
do
    for Set in all highconfidence totalmodule
    do
        WorkDir=promotor_id/large_modules/$Set
        NegDir=/data/scratch/adamst/large_modules/$Set
        Positive=$WorkDir/"$Set"_upstream3000.split.100mer.fasta
        Negative=$NegDir/"$Set"_nontarget_upstream3000.split.100mer_random_*_"$Rep".fasta
        ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
        Jobs=$(qstat | grep 'sub_dreme' | wc -l)
        while [ $Jobs -gt 20 ]
        do
            sleep 1
            printf "."
            Jobs=$(qstat | grep 'sub_dreme' | wc -l)
        done
        qsub $ProgDir/sub_dreme.sh $Positive $Negative $NegDir $Rep 0.25
    done
done
```

Combine results from multiple repeats

```bash
for Set in all highconfidence totalmodule
do
    OutDir=promotor_id/large_modules/$Set
    DremeDir=/data/scratch/adamst/large_modules/$Set
    DremeRes=$DremeDir/*dreme*/dreme.txt
    Percentage=90
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
    python $ProgDir/Combine_DREME_Repeats.py --inputs $DremeRes \
    --percentage $Percentage --outdir $OutDir
done
```

Nothing can be reproducibly pulled out, possibly due to the small input sample size
Or due to the lack of a true negative control set

Using larger module sizes

## Investigate enrichment of effector class genes in co-expression modules

### Create lists of RxLRs, CRNs, ApoplastP hits, a combined list of all effectors and a list of secreted proteins

Remove statement from beginning of removed_genes file from WGCNA done previously
Files for Genomic sets of Genes and features done previously

Create files for Module sets of Genes and features

```bash
for File in $(ls analysis/coexpression/high/merged_modules/Genes_in_*.txt)
do
    Filename=$(echo $File | rev | cut -f1 -d "/" | rev)
    Module_ID=$(echo $Filename | cut -f3 -d "_" | cut -f1 -d ".")
    echo "Processing $Module_ID"
    Module_Dir=analysis/coexpression/enrichment_large_modules/$Module_ID
    mkdir -p $Module_Dir
    Gene_Set_unsorted=$Module_Dir/Gene_Set_unsorted.txt
    Gene_Set=$Module_Dir/Gene_Set.txt
    # Copy entire gene set of module
    cp $File $Gene_Set_unsorted
    cat $Gene_Set_unsorted | sed 's/"//g' | sort | uniq > $Gene_Set
    rm $Gene_Set_unsorted
    echo "Gene set file created for $Module_ID"
    # Create RxLR list
    RxLR_Headers=analysis/RxLR_effectors/combined_evidence/P.fragariae/Bc16/Bc16_Total_RxLR_motif_hmm.txt
    RxLRs_In_Module=$Module_Dir/RxLR_IDs.txt
    cat $RxLR_Headers | grep -wf $Gene_Set | sort | uniq > $RxLRs_In_Module
    echo "RxLR file created for $Module_ID"
    # Create CRN list
    CRN_Headers=analysis/CRN_effectors/hmmer_CRN/P.fragariae/Bc16/Bc16_final_CRN.txt
    CRNs_In_Module=$Module_Dir/CRN_IDs.txt
    cat $CRN_Headers | grep -wf $Gene_Set | sort | uniq > $CRNs_In_Module
    echo "CRN file created for $Module_ID"
    # Create list of ApoplastP hits
    ApoP_Headers=analysis/ApoplastP/P.fragariae/Bc16/Bc16_Total_ApoplastP.txt
    ApoP_In_Module=$Module_Dir/ApoP_IDs.txt
    cat $ApoP_Headers | grep -wf $Gene_Set | sort | uniq > $ApoP_In_Module
    echo "ApoplastP hits file created for $Module_ID"
    # Create combined effector list
    Effectors_In_Module=$Module_Dir/Effector_IDs.txt
    cat $RxLRs_In_Module $CRNs_In_Module $ApoP_In_Module | sort | uniq > $Effectors_In_Module
    echo "Effector list created for $Module_ID"
    # Create combined secreted proteins list
    Secreted_In_Module=$Module_Dir/Secreted_IDs.txt
    Augustus_Secreted=gene_pred/combined_sigP_CQ/P.fragariae/Bc16/Bc16_secreted.txt
    ORF_Secreted=gene_pred/combined_sigP_ORF/P.fragariae/Bc16/Bc16_all_secreted_merged.txt
    cat $ORF_Secreted | sed 's/$/.t1/' > tmp_ORF.txt
    cat $Augustus_Secreted tmp_ORF.txt | grep -wf $Gene_Set | sort | uniq > $Secreted_In_Module
    rm tmp_ORF.txt
    echo "Secreted list created for $Module_ID"
done
```

### Create contigency tables for fishers exact test

```bash
RxLR_Headers=analysis/RxLR_effectors/combined_evidence/P.fragariae/Bc16/Bc16_Total_RxLR_motif_hmm.txt
Genome_RxLRs=$(cat $RxLR_Headers | sort | uniq | wc -l)
CRN_Headers=analysis/CRN_effectors/hmmer_CRN/P.fragariae/Bc16/Bc16_final_CRN.txt
Genome_CRNs=$(cat $CRN_Headers | sort | uniq | wc -l)
ApoP_Headers=analysis/ApoplastP/P.fragariae/Bc16/Bc16_Total_ApoplastP.txt
Genome_ApoP=$(cat $ApoP_Headers | sort | uniq | wc -l)
Genome_Effectors=$(cat $RxLR_Headers $CRN_Headers $ApoP_Headers | sort | uniq | wc -l)
Augustus_Secreted=gene_pred/combined_sigP_CQ/P.fragariae/Bc16/Bc16_secreted.txt
ORF_Secreted=gene_pred/combined_sigP_ORF/P.fragariae/Bc16/Bc16_all_secreted_merged.txt
Genome_Secreted=$(cat $Augustus_Secreted $ORF_Secreted | sort | uniq | wc -l)
Genes=gene_pred/annotation/P.fragariae/Bc16/Bc16_genes_incl_ORFeffectors.gene.fasta
Genome_Genes=$(cat $Genes | sort | uniq | wc -l)
for Module_Dir in $(ls -d analysis/coexpression/enrichment_large_modules/*)
do
    Module_Name=$(echo $Module_Dir | rev | cut -f1 -d "/" | rev)
    Module_RxLRs=$Module_Dir/RxLR_IDs.txt
    Module_CRNs=$Module_Dir/CRN_IDs.txt
    Module_ApoP=$Module_Dir/ApoP_IDs.txt
    Module_Effectors=$Module_Dir/Effector_IDs.txt
    Module_Secreted=$Module_Dir/Secreted_IDs.txt
    Module_Genes=$Module_Dir/Gene_Set.txt
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
    echo "Processing $Module_Name"
    python $ProgDir/Fisher_table_prep.py --Module_RxLRs $Module_RxLRs --Module_CRNs $Module_CRNs --Module_ApoP $Module_ApoP --Module_Effectors $Module_Effectors --Module_Secreted $Module_Secreted --Module_Genes $Module_Genes --Module_Name $Module_Name --Genome_RxLRs $Genome_RxLRs --Genome_CRNs $Genome_CRNs --Genome_ApoP $Genome_ApoP --Genome_Effectors $Genome_Effectors --Genome_Secreted $Genome_Secreted --Genome_Genes $Genome_Genes --OutDir $Module_Dir
    echo "$Module_Name done"
done
```

### Run Fishers Exact test to test for enrichment

Does not run on head node, throws an error of incorrect version of GLOB
I recommend running in screen due to large numbers of items in loop

```bash
qlogin
cd /home/groups/harrisonlab/project_files/phytophthora_fragariae

for Table in $(ls analysis/coexpression/enrichment_large_modules/*/*Fishertable.txt)
do
    Module_ID=$(echo $Table | rev | cut -f2 -d "/" | rev)
    Filename=$(echo $Table | rev | cut -f1 -d "/" | rev)
    Gene_Type=$(echo $Filename | cut -f2 -d "_")
    OutDir=analysis/coexpression/enrichment_large_modules/$Module_ID/Fisher_Results
    mkdir -p $OutDir
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
    echo "Running $Module_ID" "$Gene_Type"
    /home/adamst/prog/R/R-3.2.5/bin/Rscript --vanilla $ProgDir/fisherstest.R --Input_Table $Table --Output_Directory $OutDir --Module_ID $Module_ID --Gene_Type $Gene_Type
    echo "$Module_ID" "$Gene_Type done"
done
```

Parse fisher results files to fewer files by type and hypothesis being tested

```bash
for Type in up down equal
do
    Files=$(ls analysis/coexpression/enrichment_large_modules/*/Fisher_Results/enriched_"$Type".txt)
    OutDir=analysis/coexpression/enrichment_large_modules/Parsed_Fisher_Results/$Type
    mkdir -p $OutDir
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/RNA_Seq_scripts
    python $ProgDir/parse_fisher_results.py --inputs $Files --outdir $OutDir --FDR 0.05 --Types RxLR CRN ApoP Effector Secreted --Threshold 0.05
done
```

None are listed as equal, so this execution of the script fails
Summary of results

```
RxLRs:
Significantly enriched: 31
Non-Significantly enriched: 43
Significantly lacking: 0
Non-Significantly lacking: 27

CRNs:
Significantly enriched: 18
Non-Significantly enriched: 9
Significantly lacking: 0
Non-Significantly lacking: 74

ApoP:
Significantly enriched: 31
Non-Significantly enriched: 50
Significantly lacking: 0
Non-Significantly lacking: 20

Effectors:
Significantly enriched: 42
Non-Significantly enriched: 47
Significantly lacking: 0
Non-Significantly lacking: 12

Secreted
Significantly enriched: 98
Non-Significantly enriched: 2
Significantly lacking: 0
Non-Significantly lacking: 1
```
