# Phylogenetic analyses 

```bash
conda activate general_tools
```

### Find single copy busco genes

Create a list of all BUSCO IDs

```bash
    OutDir=analysis_VP/popgen/busco_phylogeny
    mkdir -p $OutDir
    BuscoDb="sordariomycetes_odb10"
    ls -1 /projects/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

```bash
screen -a
srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
# Create a folder for each busco gene
mkdir temp_busco
printf "" > analysis_VP/popgen/busco_phylogeny_2/single_hits.txt
for Busco in $(cat analysis_VP/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
OutDir=analysis_VP/popgen/busco_phylogeny_2/$Busco
mkdir -p $OutDir
# Move all single copy genes to each folder
for Fasta in $(ls ../../data/scratch/gomeza/gene_pred/BUSCO/N.*/*/*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
Strain=$(echo $Fasta | rev | cut -f7 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f8 -d '/' | rev)
FileName=$(basename $Fasta)
contig=$(cat $Fasta | grep '>' | sed 's/ <unknown description>//g' | sed 's/>//g')
echo ">$Busco:$Strain:$contig" > temp_busco/"$Busco"_"$Strain"_new_names.txt
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
python $ProgDir/replace_fasta_records.py -i $Fasta -r temp_busco/"$Busco"_"$Strain"_new_names.txt -o $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
rm temp_busco/"$Busco"_"$Strain"_new_names.txt
#cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
# Create fasta file containing all busco for alignment
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis_VP/popgen/busco_phylogeny_2/single_hits.txt
done

```

```bash
# Check for multiple hits and remove
less analysis_VP/popgen/busco_phylogeny/single_hits.txt | sort -k2 -n
```

```bash
# If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
OutDir=analysis_VP/popgen/busco_phylogeny/alignments
mkdir -p $OutDir
OrganismNum=$(cat analysis_VP/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
for Busco in $(cat analysis_VP/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
HitNum=$(cat analysis_VP/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
if [ $HitNum == $OrganismNum ]; then
cp analysis_VP/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
fi
done
```

### Gene alignments

```bash
# Submit alignment for single copy busco genes with a hit in each organism
AlignDir=analysis_VP/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
sbatch $ProgDir/mafft.sh
cd $CurDir
```


### Trim poor alignments

```bash
# Edit header name keeping BUSCO name and isolate name
cd analysis_VP/popgen/busco_phylogeny/alignments
sed -i 's/:contig.*//g' *_appended_aligned.fasta
sed -i 's/:LD.*//g' *_appended_aligned.fasta
sed -i 's/:NODE.*//g' *_appended_aligned.fasta

# New isolate names contain only numbere. I have renamed them using sed
sed -i 's/118923/US23/g' *_appended_aligned.fasta
sed -i 's/118924/US24/g' *_appended_aligned.fasta
sed -i 's/226-31/DEU26/g' *_appended_aligned.fasta
sed -i 's/227-31/NO27/g' *_appended_aligned.fasta
```

Trimming sequence alignments using Trim-Al. Note - automated1 mode is optimised for ML tree reconstruction

```bash
OutDir=analysis_VP/popgen/busco_phylogeny/trimmed_alignments
mkdir -p $OutDir
for Alignment in $(ls analysis_VP/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
echo $Alignment
trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
done
```

### Randomized Axelerated Maximum Likelihood

```bash
screen -a

    for Alignment in $(ls analysis_VP/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
        sleep 10s
        Prefix=$(basename $Alignment | cut -f1 -d '_')
        OutDir=analysis_VP/popgen/busco_phylogeny/RAxML/$Prefix
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
        sbatch $ProgDir/RAxML.sh $Alignment $Prefix $OutDir
    done
```

### Astral

```bash
OutDir=analysis_VP/popgen/busco_phylogeny/ASTRAL
mkdir -p $OutDir

#Concatenate best trees
#cat analysis_VP/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" > $OutDir/Nd_phylogeny.appended2.tre
cat analysis_VP/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  > $OutDir/Nd_phylogeny.appended3.tre


US24192
# Contract low support brances (below 10% bootsrap support)
nw_ed $OutDir/Nd_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Nd_phylogeny.appended.trimmed.tre

# Calculate combined tree
ProgDir=/scratch/software/ASTRAL/ASTRAL-5.7.1/Astral
java -jar $ProgDir/astral.5.7.1.jar -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.tre 2> $OutDir/Nd_phylogeny.consensus.log
# Score the resulting tree
java -jar $ProgDir/astral.5.7.1.jar -q $OutDir/Nd_phylogeny.consensus.tre -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.scored.tre 2> $OutDir/Nd_phylogeny.consensus.scored.log
```




GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.
* Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

```bash
cat Alt_phylogeny.consensus.scored.geneious.tre | sed 's/:2/:1/g' > Alt_phylogeny.consensus.scored.geneious2.tre
```