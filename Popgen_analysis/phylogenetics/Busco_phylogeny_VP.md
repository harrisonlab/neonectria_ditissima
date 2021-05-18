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