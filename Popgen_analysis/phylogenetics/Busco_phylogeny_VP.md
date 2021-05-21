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
cat analysis_VP/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  > $OutDir/Nd_phylogeny.appended.tre

# Contract low support brances (below 10% bootsrap support)
nw_ed $OutDir/Nd_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Nd_phylogeny.appended.trimmed.tre

# Calculate combined tree
ProgDir=/scratch/software/ASTRAL/ASTRAL-5.7.1/Astral
java -jar $ProgDir/astral.5.7.1.jar -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.tre 2> $OutDir/Nd_phylogeny.consensus.log
# Score the resulting tree
java -jar $ProgDir/astral.5.7.1.jar -q $OutDir/Nd_phylogeny.consensus.tre -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.scored.tre 2> $OutDir/Nd_phylogeny.consensus.scored.log

# Calculate combined tree
ProgDir=/scratch/software/ASTRAL/ASTRAL-5.7.1/Astral
java -jar $ProgDir/astral.5.7.1.jar -i $OutDir/Nd_phylogeny.appended.trimmed.tre -o $OutDir/Nd_phylogeny.trimmed.consensus.tre 2> $OutDir/Nd_phylogeny.trimmed.consensus.log
# Score the resulting tree
java -jar $ProgDir/astral.5.7.1.jar -q $OutDir/Nd_phylogeny.trimmed.consensus.tre -i $OutDir/Nd_phylogeny.appended.trimmed.tre -o $OutDir/Nd_phylogeny.consensus.trimmed.scored.tre 2> $OutDir/Nd_phylogeny.consensus.trimmed.scored.log
```

Manual edition of the final consensus tree is needed 

```
Step 1: Download consensus tree to local machine

Step 2: Import into geneious and export again in newick format to get around polytomy branches having no branch length.

Step 3: Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed
```
```bash
cat Desktop/Astral_tree/Nd_phylogeny.consensus.scored_newick.tre | sed 's/:2/:1/g' > Desktop/Astral_tree/Nd_phylogeny.consensus.scored_newick2.tre
```


## Plot best scored tree

GGtree was used to make a plot. Tutorial tips: https://bioconnector.org/r-ggtree.html

R version > 4.0

```r
setwd("~/Desktop/Astral_tree")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)

packageurl <- "http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.5.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

remotes::install_github("YuLab-SMU/tidytree")


tree <- read.tree("~/Desktop/Astral_tree/Nd_phylogeny.consensus.scored_newick2.tre")
# Basic tree
t<-ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
# Load table with ID, Country, Host and Specie
mydata <- read.csv("~/Desktop/Astral_tree/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$ID
# Sample name should match in tree and csv
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]
# Format nodes by values
nodes <- data.frame(t$data)
# Core tree
t <- ggtree(tree, aes(linetype=nodes$support))
# Add scalebar
t <- t + geom_treescale(offset=-1.0, fontsize = 3) 









# Adjust terminal branch lengths:
# branches <- t$data
# tree$edge.length[branches$isTip] <- 0.1


#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb



#51 is the node of your outgroup?
tree$edge.length[tree$edge.length == 1] <- 0
tree$edge.length[51] <- 0


t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree


# Adjust terminal branch lengths:
branches <- t$data

branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
# tree$edge.length[tree$edge.length == 1] <- 0
# t <- ggtree(tree, aes(linetype=nodes$support))
#Tree <- branches$branch.length

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels


# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, aes(color=Source), size=3, hjust=0, align=T, offset = +0.1) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels
tree_mod <- data.frame(t$data)
tree_mod$label <- tips$pathotype
t <- t + geom_tiplab(data=tree_mod, aes(label=label, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +5.0) +
scale_color_manual(values=c("gray39","black"))

tips$MAT <- factor(tips$MAT)
# t <- t + geom_tippoint(data=tips, aes(shape=MAT), size=2)
t <- t + geom_tiplab(data=tips, aes(label=MAT, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +3.5) +
scale_color_manual(values=c("gray39","black"))



# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=42, label='sect. Alternaria', align=T, colour='black', offset=-1.5)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=51, label='tenuissima clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=45, label='arborescens clade', align=T, colour='black', offset=-4.5)
t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=9.5)
t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='', colour='NA', offset=17.5)

# Save as PDF and force a 'huge' size plot
# t <- ggsave("expanded/Fig3_busco_phylogeny.pdf", width =30, height = 30, units = "cm", limitsize = FALSE)
t <- ggsave("expanded/Fig3_busco_phylogeny.tiff", width =30, height = 30, units = "cm", limitsize = FALSE)

````
