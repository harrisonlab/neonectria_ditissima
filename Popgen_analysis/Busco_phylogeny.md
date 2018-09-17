# 0 Busco Assessment of Neonectria ditissima isolates

New Zealand genomes. Edit headers and convert to fasta.

```bash
cat LDPK01.1.fsa_nt | sed 's/ Neo.*//g' > LDPK01.1.fsa_nt.fasta
cat LDPL01.1.fsa_nt | sed 's/ Neo.*//g' > LDPL01.1.fsa_nt.fasta
```
```bash
  for Assembly in $(ls repeat_masked/Nz_genomes/*/*_nt.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

# 1 Find single copy busco genes in N.ditisima assemblies

```bash
  for Strain in Ag02 Ag04 Ag05 Ag06 Ag08 Ag09_A Ag11_A Ag11_B Ag11_C BGV344 Hg199 ND8 ND9 P112 OPC304 R0905_canu_2017_v2 R37-15 R39-15 R41-15 R42-15 R45-15 R68-17 R6-17-2 R6-17-3; do
    for Assembly in $(ls repeat_masked/N.ditissima/$Strain/*unmasked.fa); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    #ProgDir=/home/gomeza/git_repos/emr_repos/tools/gene_prediction/busco
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco_armita/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
    done
  done
```

For Hg199 and R0905, I used the genome assemblies with the best busco prediction, regardless of contig number. These were the spades assembly of Hg199 and the canu pilon5 of R0905.

Create a list of all BUSCO IDs

```bash
OutDir=analysis/popgen/busco_phylogeny2
mkdir -p $OutDir
BuscoDb="sordariomyceta_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for alignment later.

```bash
printf "" > analysis/popgen/busco_phylogeny2/single_hits.txt
  for Busco in $(cat analysis/popgen/busco_phylogeny2/all_buscos_*.txt); do
  echo $Busco
  OutDir=analysis/popgen/busco_phylogeny2/$Busco
  mkdir -p $OutDir
  for Fasta in $(ls gene_pred/busco/N.*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
  Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
  Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
  FileName=$(basename $Fasta)
  cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
  done
  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny2/single_hits.txt
done
```

Check for multiple hits

```bash
less single_hits.txt | sort -k2 -n
```
Three busco genes gave me unexpected multiple hits. Therefore I removed them for my analysis

```bash
#These three gave multiple hits in the preliminary busco analysis.
#rm gene_pred/busco/N.ditissima/*/assembly/*/single_copy_busco_sequences/EOG093314TJ*
#rm gene_pred/busco/N.ditissima/*/assembly/*/single_copy_busco_sequences/EOG09330AIP*
#rm gene_pred/busco/N.ditissima/*/assembly/*/single_copy_busco_sequences/EOG09330VPU*

rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG093318S0*
rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG093305B4*
rm gene_pred/busco/N.ditissima/*/*/*/single_copy_busco_sequences/EOG0933010Y*
#These are different to the previous ones.
```

If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis/popgen/busco_phylogeny/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
  if [ $HitNum == $OrganismNum ]; then
    cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
  fi
  done
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
AlignDir=analysis/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir
ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
qsub $ProgDir/sub_mafft_alignment.sh
cd $CurDir
```

Trimming sequence alignments using Trim-Al
* Note - automated1 mode is optimised for ML tree reconstruction

```bash
  OutDir=analysis/popgen/busco_phylogeny/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    echo $Alignment
    trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
  done
```
Edit header name keeping BUSCO name and isolate name
```bash
cd analysis/popgen/busco_phylogeny/trimmed_alignments
sed -i 's/_contigs_.*//g' *_appended_aligned_trimmed.fasta
sed -i 's/:LD.*//g' *_appended_aligned_trimmed.fasta
sed -i 's/N.ditissima_//g' *_appended_aligned_trimmed.fasta
```

```bash
for Alignment in $(ls analysis/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
# printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis/popgen/busco_phylogeny/RAxML/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```


Run Astral to build a consensus phylogeny from a collective set of
"best phylogenies" from each BUSCO locus

* Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."

Tutorial tips:
https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-with-unresolved-gene-trees


```bash
OutDir=analysis/popgen/busco_phylogeny/ASTRAL
mkdir -p $OutDir
# cat analysis/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.* > $OutDir/Pcac_phylogeny.appended.tre
# Taxa names were noted to be incorrect at this point and were corrected
cat analysis/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" > $OutDir/Nd_phylogeny.appended.tre
# InTree=$(ls /home/armita/prog/Astral/Astral/test_data/song_primates.424.gene.tre)
# -
# Trimm back branches that have less than 10% bootstrap support for each tree
# in the given file
# -
/home/armita/prog/newick_utilities/newick_utils/src/nw_ed $OutDir/Nd_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Nd_phylogeny.appended.trimmed.tre
# -
# Calculate combined tree
# -
ProgDir=/home/gomeza/prog/Astral/Astral
# java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Pcac_phylogeny.appended.trimmed.tre -o $OutDir/Pcac_phylogeny.consensus.tre | tee 2> $OutDir/Pcac_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.tre | tee 2> $OutDir/Nd_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -q $OutDir/Nd_phylogeny.consensus.tre -i $OutDir/Nd_phylogeny.appended.tre -o $OutDir/Nd_phylogeny.consensus.scored.tre 2> $OutDir/Nd_phylogeny.consensus.scored.log
```


GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.

```r
setwd("/Users/gomeza/Downloads/Ndit/ASTRAL")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("/Users/armita/Downloads/Aalt/ASTRAL/Alt_phylogeny.consensus.scored.geneious2.tre")


mydata <- read.csv("/Users/armita/Downloads/Aalt/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$Isolate
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]



t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree
# Adjust terminal branch lengths:
branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
#Tree <- branches$branch.length
#rescale_tree(t, branches$branch.length)



t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels



# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Format nodes by values
nodes <- data.frame(t$data)
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

# Add in a further set of labels
# tree_mod <- tree
# tree_mod$tip.label <- mydata$Source
# t <- t + geom_tiplab(data=tree_mod, aes(label=label), size=2, offset = +1)

# Annotate a clade with a bar line
t <- t + geom_cladelabel(node=45, label='sect. Alternaria', align=T, colour='black', offset=-1.5)
t <- t + geom_cladelabel(node=68, label='gaisen clade', align=T, colour='black', offset=-4.5)
t <- t + geom_cladelabel(node=53, label='tenuissima clade', align=T, colour='black', offset=-4.5)
t <- t + geom_cladelabel(node=47, label='arborescens clade', align=T, colour='black', offset=-4.5)

# Save as PDF and force a 'huge' size plot
ggsave("tree4.pdf", width =30, height = 30, units = "cm", limitsize = FALSE)

````


<!--

```bash
# For closely related organisms (same species etc.): identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites
# (avoid alignments with low homology and lots of phylogenetically uninformative singletons).
# For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity
# (e.g. 0.1<Pi<0.4), looking for genes with the lowest number of segregating sites.
AlignDir=analysis/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir

# pip install dendropy --user
for Alignment in $(ls *aligned.fasta); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
python $ProgDir/calculate_nucleotide_diversity.py $Alignment
Busco=$(echo $Alignment | cut -f1 -d '_')
mv sequence_stats.txt "$Busco"_seqeunce_stats.txt
mv excel_stats.txt "$Busco"_excel_stats.txt
mkdir -p ../phylogeny
## Copy FASTA files of the aligments into a new directory
cp $Alignment ../phylogeny/.
done

cd $CurDir
```

Visually inspect the alignments of selected genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed
```


##PartitionFinder (nucleotide sequence evolution model)

```bash
cd analysis/popgen/busco_phylogeny/phylogeny

config_template=/home/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
ct=$(basename "$config_template")

mkdir NEXUS

# prepare directory for PartitionFinder run:
for f in $(ls *fasta); do
sed -i 's/:/_/g' $f
c="$(cat $f | awk 'NR%2==0' | awk '{print length($1)}' | head -1)"
p="${f%.fasta}.phy"
n="${f%.fasta}.NEXUS"
dir="${f%.fasta}"

mkdir $dir
cp $config_template $dir/.

# Substitute the name of the alignment file and the sequence length in the config file to become correct for the current run.
sed -i 's,^\(alignment = \).*,\1'"$p;"',' $dir/$ct
sed -i 's,^\(Gene1_pos1 = \).*,\1'"1-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos2 = \).*,\1'"2-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos3 = \).*,\1'"3-$c\\\3;"',' $dir/$ct

# Convert FASTA to phylip for the Partition Finder run
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
$ProgDir/fasta2phylip.pl $f>$p
mv $p $dir

# Convert FASTA to NEXUS for the BEAST run
$ProgDir/Fasta2Nexus.pl $f>$dir/$n

#Problems running PartitionFinder on the cluster. May have to be run locally on your Mac or Windows machine.
# qsub $ProgDir/sub_partition_finder.sh $dir
done
```

Partition finder wasnt run on the cluster. As such fasta alignment files were
downloaded to the local machine where partitionfinder was run
patritionfinder2 was downloaded from:
http://www.robertlanfear.com/partitionfinder/

and the anaconda libraries to support it were downloaded from:
https://www.continuum.io/downloads#macos


copy the fasta files and the partitionfinder config files to
your local computer

```bash
cd Users/armita/Downloads
scp -r cluster:/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny .
```

Alignments were loaded into Geneious where they were visualised and manually sorted into
three categories:
* Good - All sequences present no trimming needed
* Trim - All sequences present short regions may need trimming from the beginning / end of the alignment before use in phylogenetics
* Bad - a region of one or more sequences is missing or the sequences / alignment is not appropriate for phylogenetics

These alignments were then exported from Geneious into the following folders:

```bash
cd Users/armita/Downloads/phylogeny
mkdir good_alignments
mkdir trim_alignments
mkdir bad_alignments
```

Alignments within the "good alignments" directory were taken forward for further
analysis

```bash
  for Dir in $(ls -d *_alignments); do
    for Alignment in $(ls $Dir/*_appended_aligned.phy); do
      Prefix=$(echo $Alignment | cut -f2 -d '/' | sed 's/.phy//g')
      echo $Prefix
      cp $Prefix/$Prefix.NEXUS $Dir/$Prefix/.
      cp -r $Prefix $Dir/.
      /Users/armita/anaconda2/bin/python ../partitionfinder-2.1.1/PartitionFinder.py $Dir/$Prefix --no-ml-tree --force-restart
    done
  done > log.txt
```


Upload partition models back to the cluster:

```bash
ClusterDir=/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny
scp -r bad_alignments cluster:$ClusterDir/.
```


## Preparing to run BEAST


Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
create an XML input file using BEAUTi, with StarBeast template.

Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.

Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub -
as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

StarBeast settings used here:
* Substitution rate: default HKY
* Strict clock
* Species Tree Population Size: Linear with constant root
* Yule prior on species tree
* Chain length: 300 million (this may vary, change run convergence with Tracer during the run to establish the number of iterations required
* Tracer: /home/sobczm/bin/beast/Tracer_v1.6/bin/tracer
some runs may never converge)
* Store every: 10000

```bash

cd /home/groups/harrisonlab/project_files/idris


for File in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/analysis/best_scheme.txt); do
Busco=$(echo $File | cut -f6 -d '/' | cut -f1 -d '_')
Model=$(cat $File | grep -A1 'Best Model' | tail -n1 | cut -f2 -d '|')
printf "$Busco\t$Model\n"
done

# Edit NEXUS files:
for Nexus in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*_appended_aligned.NEXUS); do
  sed -i -r "s/^.*_P\./P./g" $Nexus
  sed -i -r "s/_contig.*\t/\t/g" $Nexus
  sed -i -r "s/_NODE.*\t/\t/g" $Nexus
done

# OUtputs of partitionfinder were used to set models
# of DNA evolution in Beauti, as described on:
# http://www.robertlanfear.com/partitionfinder/faq/#toc-beast
# CHain length was modified from 10000000 to 500000000 as determined
# by a first run of beast where tracer reported the estimated sasmple size to be below 100 (3) - increase by 50 fold.

# Run Beauti
NexusFiles=$(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*.NEXUS | sed -e 's/^/ -nex /g' | tr -d '\n')
OutFile=$(echo $Nexus | sed 's/.NEXUS/.xml/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beauti -template StarBeast.xml $NexusFiles




qlogin -pe smp 8
InXML=analysis/popgen/busco_phylogeny/phylogeny/Pcac_beauti_starBEAST2.xml
OutDir=$(dirname $InXML)"/BEAST4"
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beast -threads 8 -prefix $OutDir $InXML > $OutDir/log.txt
# java -Djava.library.path="C:\Program Files (x86)\Common Files\libhmsbeagle-1.0" -jar "/BEAST175/lib/beast.jar"

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator
for Tree in $(ls $OutDir/*.trees); do
BurnIn=10 # percentage of states to be considered as burnin
SumTree=$(echo $Tree | sed 's/.trees/_summary.tree/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/treeannotator -heights median -burnin $BurnIn $Tree $SumTree
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
FigTree=/home/sobczm/bin/FigTree_v1.4.2/bin/figtree
$FigTree

``` -->






```bash
# For closely related organisms (same species etc.): identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites
# (avoid alignments with low homology and lots of phylogenetically uninformative singletons).
# For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity
# (e.g. 0.1<Pi<0.4), looking for genes with the lowest number of segregating sites.
AlignDir=analysis/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir

# pip install dendropy --user
for Alignment in $(ls *aligned.fasta); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
python $ProgDir/calculate_nucleotide_diversity.py $Alignment
Busco=$(echo $Alignment | cut -f1 -d '_')
mv sequence_stats.txt "$Busco"_sequence_stats.txt
mv excel_stats.txt "$Busco"_excel_stats.txt
mkdir -p ../phylogeny
## Copy FASTA files of the aligments into a new directory
cp $Alignment ../phylogeny/.
done

cd $CurDir
```

Visually inspect the alignments of select genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed
```


##PartitionFinder (nucleotide sequence evolution model)

```bash
cd analysis/popgen/busco_phylogeny/phylogeny

config_template=/home/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
ct=$(basename "$config_template")

mkdir NEXUS

# prepare directory for PartitionFinder run:
for f in $(ls *fasta); do
sed -i 's/:/_/g' $f
c="$(cat $f | awk 'NR%2==0' | awk '{print length($1)}' | head -1)"
p="${f%.fasta}.phy"
n="${f%.fasta}.NEXUS"
dir="${f%.fasta}"

mkdir $dir
cp $config_template $dir/.

# Substitute the name of the alignment file and the sequence length in the config file to become correct for the current run.
sed -i 's,^\(alignment = \).*,\1'"$p;"',' $dir/$ct
sed -i 's,^\(Gene1_pos1 = \).*,\1'"1-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos2 = \).*,\1'"2-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos3 = \).*,\1'"3-$c\\\3;"',' $dir/$ct

# Convert FASTA to phylip for the Partition Finder run
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
$ProgDir/fasta2phylip.pl $f>$p
mv $p $dir

# Convert FASTA to NEXUS for the BEAST run
$ProgDir/Fasta2Nexus.pl $f>$dir/$n

#Problems running PartitionFinder on the cluster. May have to be run locally on your Mac or Windows machine.
# qsub $ProgDir/sub_partition_finder.sh $dir
done
```

Partition finder wasnt run on the cluster. As such fasta alignment files were
downloaded to the local machine where partitionfinder was run
patritionfinder2 was downloaded from:
http://www.robertlanfear.com/partitionfinder/

and the anaconda libraries to support it were downloaded from:
https://www.continuum.io/downloads#macos


copy the fasta files and the partitionfinder config files to
your local computer

```bash
cd Users/armita/Downloads
scp -r cluster:/data/scratch/gomeza/analysis/popgen/busco_phylogeny/phylogeny .
```

Alignments were loaded into Geneious where they were visualised and manually sorted into
three categories:
* Good - All sequences present no trimming needed
* Trim - All sequences present short regions may need trimming from the beginning / end of the alignment before use in phylogenetics
* Bad - a region of one or more sequences is missing or the sequences / alignment is not appropriate for phylogenetics

These alignments were then exported from Geneious into the following folders:

```bash
cd Users/armita/Downloads/phylogeny
mkdir good_alignments
mkdir trim_alignments
mkdir bad_alignments
```

Alignments within the "good alignments" directory were taken forward for further
analysis

```bash
  for Dir in $(ls -d *_alignments); do
    for Alignment in $(ls $Dir/*_appended_aligned.phy); do
      Prefix=$(echo $Alignment | cut -f2 -d '/' | sed 's/.phy//g')
      echo $Prefix
      cp $Prefix/$Prefix.NEXUS $Dir/$Prefix/.
      cp -r $Prefix $Dir/.
      /Users/armita/anaconda2/bin/python ../partitionfinder-2.1.1/PartitionFinder.py $Dir/$Prefix --no-ml-tree --force-restart
    done
  done > log.txt
```


Upload partition models back to the cluster:

```bash
ClusterDir=/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny
scp -r bad_alignments cluster:$ClusterDir/.
```


## Preparing to run BEAST


Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
create an XML input file using BEAUTi, with StarBeast template.

Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.

Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub -
as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

StarBeast settings used here:
* Substitution rate: default HKY
* Strict clock
* Species Tree Population Size: Linear with constant root
* Yule prior on species tree
* Chain length: 300 million (this may vary, change run convergence with Tracer during the run to establish the number of iterations required
* Tracer: /home/sobczm/bin/beast/Tracer_v1.6/bin/tracer
some runs may never converge)
* Store every: 10000

```bash

cd /home/groups/harrisonlab/project_files/idris


for File in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/analysis/best_scheme.txt); do
Busco=$(echo $File | cut -f6 -d '/' | cut -f1 -d '_')
Model=$(cat $File | grep -A1 'Best Model' | tail -n1 | cut -f2 -d '|')
printf "$Busco\t$Model\n"
done

# Edit NEXUS files:
for Nexus in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*_appended_aligned.NEXUS); do
  sed -i -r "s/^.*_P\./P./g" $Nexus
  sed -i -r "s/_contig.*\t/\t/g" $Nexus
  sed -i -r "s/_NODE.*\t/\t/g" $Nexus
done

# OUtputs of partitionfinder were used to set models
# of DNA evolution in Beauti, as described on:
# http://www.robertlanfear.com/partitionfinder/faq/#toc-beast
# CHain length was modified from 10000000 to 500000000 as determined
# by a first run of beast where tracer reported the estimated sasmple size to be below 100 (3) - increase by 50 fold.

# Run Beauti
NexusFiles=$(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*.NEXUS | sed -e 's/^/ -nex /g' | tr -d '\n')
OutFile=$(echo $Nexus | sed 's/.NEXUS/.xml/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beauti -template StarBeast.xml $NexusFiles




qlogin -pe smp 8
InXML=analysis/popgen/busco_phylogeny/phylogeny/Pcac_beauti_starBEAST2.xml
OutDir=$(dirname $InXML)"/BEAST4"
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beast -threads 8 -prefix $OutDir $InXML > $OutDir/log.txt
# java -Djava.library.path="C:\Program Files (x86)\Common Files\libhmsbeagle-1.0" -jar "/BEAST175/lib/beast.jar"

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator
for Tree in $(ls $OutDir/*.trees); do
BurnIn=10 # percentage of states to be considered as burnin
SumTree=$(echo $Tree | sed 's/.trees/_summary.tree/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/treeannotator -heights median -burnin $BurnIn $Tree $SumTree
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
FigTree=/home/sobczm/bin/FigTree_v1.4.2/bin/figtree
$FigTree

```
