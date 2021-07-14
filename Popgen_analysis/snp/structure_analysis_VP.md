# Structure analysis

```bash
# Set paths
input=Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/NoNMaj_NoAg11_B
scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/Population_genomics
pgdspid=/scratch/software/PGDSpider_2.1.1.5

input_file=R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.vcf 
#Downsample SNPs for Structure analysis as enough information in 10% of the loci (and more not informative because of linkage).
#In certain cases, when small number of markers detected, this step is unnecessary and all can be retained.
vcfrandomsample --rate 0.1 $input_file > R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.vcf
less R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.vcf  | grep -v '^#' | wc -l
#44740
```

```bash
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
##Part of the path is missing here: should it be
input_file=$input/R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.vcf
#Prepare population definition file. Each individual = new population
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $scripts/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$dir/$list_file"',' vcf_to_structure_haploid.spid
#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names


filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
# PGDSpider needs to be installed using conda
PGDSpider2-cli -inputfile $input_file -inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid.spid
#N.ditissima_contigs_unmasked_filtered.recode_annotated.struc in project file so had to copy to SNP_calling folder
cp *.struc Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/NoNMaj_NoAg11_B

#Note - Iteractions can be changed in the mainparams file.
#Important note - I used BURNIN 1000 NUMREPS 10000 for the first time. This gives large variance in lnPD, inconclusive run
# Minimum Burnin reps must be 100000. Number of reps must be between 10000 and 1000000. This will require days.
#Run replicate STRUCTURE runs, with K from 1 to 10
# BURNIN 100000 NUMREPS 100000

#OutDir=Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/NoNMaj_NoAg11_B/structure
cd $input
mkdir -p structure

# inputs are struc file, ploidy, K and number of reps
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 1 10 structure
sbatch -p long $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 2 10 structure
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 3 10 structure
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 4 10 structure
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 5 10 structure
sbatch -p long $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 6 10 structure
sbatch -p long $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 7 10 structure
sbatch -p long $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 8 10 structure
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 9 10 structure
sbatch -p himem $scripts/execute_structure.sh R0905_good_contigs_unmasked_FINAL_filtered.recode_annotated.subsampled.struc 1 10 10 structure
```

## Analyze STRUCTURE output

Generate a folder containing all STRUCTURE output files for all K analyzed

```bash
mkdir structureHarvester
for d in gomeza_7*/*
do
cp $d/*_f $PWD/structureHarvester
done
```

#Tidy working directory
mv structure_* analysis/popgen/SNP_calling/structure

# structureHarvester - summarise the results

```bash
screen -a
srun --partition long --mem-per-cpu 20G --cpus-per-task 10 --pty bash

harvester=/home/gomeza/prog/structureHarvester/structureHarvester.py
$harvester --dir=structureHarvester --out=structureHarvester --evanno --clumpp
```

# CLUMPP - permute the results

```bash
# This will take a few days. It is recommended to run ind and pop in separated runs

cd structureHarvester
clumpp=/home/gomeza/prog/CLUMPP_Linux64.1.1.2
cp $clumpp/paramfile_ind ./
mv paramfile_ind paramfile

#Options fed to CLUMPP
#-i: indfile from StructureHarvester output
#-p: popfile from StructureHarvester output
#-o: output Q matrix for distruct input
#-k: K value (number of clusters tested)

###!!! Options to be changed in each analysis manually
#c: number of individuals (change according to STRUCTURE mainparam file)
#r: number of replicate runs
#s: minimum number of population clusters (K) tested
#f: maximum number of population clusters (K) tested
c=30
r=10
s=1
f=10
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done

# Going to run this in a separate screen!!
screen -a
srun --partition himem --mem-per-cpu 12G --cpus-per-task 24 --pty bash

mkdir CLUMPP_pop
cp *indfile CLUMPP_pop
cp *popfile CLUMPP_pop/
cd CLUMPP_pop/

cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile

c=11
r=10
s=1
f=10
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c $c -r $r
done

```









#Key options in the paramfile
# DISTRUCT to visualise the results
###!!!! Options to be changed in each analysis manually
#-M number of populations assigned in the Structure input file
#-N number of individuals
m=11
n=30
#-K K value
#-p input file (population q's)
#-i input file (individual q's)
#-a input file (labels atop figure)
#-b input file (labels below figure)
#-o output file
distruct=/home/gomeza/prog/distruct1.1
cp $distruct/drawparams ./
for i in $(seq $s $f) #input range of K values tested
do
$distruct/distructLinux1.1 -i K$i.indivq -p K$i.popq -a $names -o K$i.ps -k $i -M $m -N $n -K $i
done

#Output is a number of PostScript files showing the average proportion of each
#individual's genome belonging to a given cluster and allowing some ancestry inference based on
#the most likely true number of population clusters as summarised by StructureHarvester.

-





cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile
c=8
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c $c -r $r
done
#Key options in the paramfile
# DISTRUCT to visualise the results
###!!!! Options to be changed in each analysis manually
#-M number of populations assigned in the Structure input file
#-N number of individuals
m=8
n=26
#-K K value
#-p input file (population q's)
#-i input file (individual q's)
#-a input file (labels atop figure)
#-b input file (labels below figure)
#-o output file
distruct=/home/sobczm/bin/distruct1.1
cp $distruct/drawparams ./
for i in $(seq $s $f) #input range of K values tested
do
$distruct/distructLinux1.1 -i K$i.indivq -p K$i.popq -a $names -o K$i.ps -k $i -M $m -N $n -K $i
done

#Output is a number of PostScript files showing the average proportion of each
#individual's genome belonging to a given cluster and allowing some ancestry inference based on
#the most likely true number of population clusters as summarised by StructureHarvester.
