PIPELINE for Mauve whole genome alignments and subsequent CLONALFRAME analysis of regions of recombination

# Use move_contigs to order genomes based on reference(Use ssh -Y cluster to ensure trusted X11 connection)

for GENOME in analysis/genome_alignment/mauve/N.ditissima/RS305p/LDPK01.1.fsa_nt.fasta  ; do
GENOME_SHORT=$(basename $GENOME | sed s/".fasta"//g )
echo $GENOME_SHORT
mkdir "$GENOME_SHORT"_dir
REFERENCE=analysis/genome_alignment/mauve/N.ditissima/RS305p/N.ditissima_contigs_unmasked.fa
echo $REFERENCE

java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME
done

java -Xmx500m -cp /home/hulinm/local/src/mauve_2.3.1/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME


##copy all fasta to align directory



# Generate alignment of genome sequences using progressive mauve - remember to use linux-64
GENOMES=$(ls /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/pss/align/removerogue/*.fasta)
/home/hulinm/local/src/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve --output=/home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/pss/align/removerogue/pss.xmfa $GENOMES

# Strip alignment of non-core parts to create core genome alignment. Minimum length of LCB is user defined as 1000
/home/hulinm/local/src/mauve_2.3.1/stripSubsetLCBs /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1.xmfa /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1.xmfa.bbcols /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1_core.xmfa 1000

# Replace file path with nothing so identifier for each strain is short version


sed -e 's/\/home\/hulinm\/pseudomonas_data\/pseudomonas\/analysis\/mauve\/race1\/new\/align\///g' /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1_core.xmfa  > /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1_core_align.xmfa
sed s/".fasta"//g /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1_core_align.xmfa   > /home/hulinm/pseudomonas_data/pseudomonas/analysis/mauve/race1/new/align/race1_core.xmfa



# Convert xmfa to concatenated fasta (all LCBs concatenated)
perl /home/hulinm/git_repos/tools/analysis/python_effector_scripts/alignment_convert.pl -i race1_core.xmfa -o race1.fasta -f fasta -c

# Use Gblocks to strip alignment of ambiguous seq

# Convert fasta to phylip format
perl /home/hulinm/git_repos/tools/analysis/python_effector_scripts/alignment_convert.pl -i race1.fasta -o race1.phy  -f phylip -g fasta

qsub /home/hulinm/git_repos/pseudomonas/orthomcl/sub_raxml.sh race1.phy




### CLONALFRAME

#Remove additional information in headers to make compatible with ClonalFrameML. Then remove blank lines
python /home/hulinm/git_repos/tools/analysis/python_effector_scripts/edit_xmfa.py race1_core.xmfa > race1_core
sed -i '/^$/d' race1_core


#### make sure names are OK - same in xmfa and tree ####

 sed 's/Ps_myricae/Ps_myricaeAZ8448/g' RAxML_bestTree.out | sed 's/9657\/1-486/9657/g' | sed 's/9629\/1-486/9629/g' | sed 's/5300\/1-486/5300/g' | sed 's/Ps_rhaphio/Ps_rhaphiolepidis/g'  | sed 's/9646\/1-486/9646/g' | sed 's/5269\/1-486/5269/g' | sed 's/5244\/1-486/5244/g' | sed 's/2341\/1-486/2341/g' > tree




# ClonalFrameML. Importation status gives you likely regions of recombinant origin


qsub /home/hulinm/git_repos/pseudomonas/orthomcl/sub_clonalframeML.sh tree2 race1_core R1



# Create graph
Rscript /home/hulinm/local/src/ClonalFrameML/src/cfml_results.R R1
