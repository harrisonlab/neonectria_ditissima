#!/bin/bash
input=/data/scratch/gomeza

input=/data/scratch/gomeza/last_structure
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
#Downsample SNPs for Structure analysis as enough information in 10% of the loci
#(and more not informative because of linkage)
#Already thinned in PopGen.md
#/home/sobczm/bin/vcflib/bin/vcfrandomsample \
#--rate 0.1 Fus2_canu_contigs_unmasked_filtered.vcf > Fus2_canu_contigs_unmasked_filtered_subsampled.vcf
#Run STRUCTURE analysis to test for the presence of
#K (population clusters) 1 to 11, with 5 replicates for each K run consecutively.

#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
##Part of the path is missing here: should it be
input_file=$input/N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
input_file=$input/N.ditissima_contigs_unmasked_filtered.vcf
#Prepare population definition file. Each individual = new population
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$dir/$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
dos2unix $outfile

#N.ditissima_contigs_unmasked_filtered.recode_annotated.struc in project file so had to copy to SNP_calling folder
cp N.ditissima_contigs_unmasked_filtered.recode_annotated.struc analysis/popgen/SNP_calling

#Note - Iteractions can be changed in the mainparams file.
#Important note - I used BURNIN 1000 NUMREPS 10000 for the first time. This gives large variance in lnPD, inconclusive run
# Minimum Burnin reps must be 100000. Number of reps must be between 10000 and 1000000. This will require days.
#Run replicate STRUCTURE runs, with K from 1 to 10
# BURNIN 100000 NUMREPS 1000
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 1 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 2 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 3 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 4 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 5 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 6 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 7 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 8 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 9 5
qsub $scripts/execute_structure.sh $input/N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 10 5

# BURNIN 100000 NUMREPS 10000, with K from 1 to 5
qsub $scripts/execute_structure2.sh N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 1 5
qsub $scripts/execute_structure2.sh N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 2 5
qsub $scripts/execute_structure2.sh N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 3 5
qsub $scripts/execute_structure2.sh N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 4 5
qsub $scripts/execute_structure2.sh N.ditissima_contigs_unmasked_filtered.recode_annotated.struc 1 5 5

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
mkdir structureHarvester
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester
done

#Tidy working directory
mv structure_* analysis/popgen/SNP_calling/structure

# structureHarvester - summarise the results
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=$input/analysis/popgen/SNP_calling/structure/structureHarvester --out=$input/analysis/popgen/SNP_calling/structure/structureHarvester --evanno --clumpp

$harvester --dir=$input/structureHarvester --out=$input/structureHarvester --evanno --clumpp

# CLUMPP - permute the results
cd analysis/popgen/SNP_calling/structure/structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
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
c=24
r=5
s=1
f=5
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done
cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile


for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c $c -r $r
done
#Key options in the paramfile
# DISTRUCT to visualise the results
###!!!! Options to be changed in each analysis manually
#-M number of populations assigned in the Structure input file
#-N number of individuals
m=24
n=24
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

-----------------------------------

Last structure analysis. mainparams and extraparams were edited adding origin information

#!/bin/bash
input=structure_analysis
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
#Move the file to the right directory
cp analysis/popgen/SNP_calling3/Hg199_contigs_unmasked_filtered.recode_annotated.vcf structure_analysis/
#input_file=$input/Hg199_contigs_unmasked_filtered.recode_annotated.vcf
#Remove outgroup. Create a cut-down VCF and filter it
vcflib=/home/sobczm/bin/vcflib/bin
cd $input
$vcflib/vcfremovesamples Hg199_contigs_unmasked_filtered.recode_annotated.vcf NMaj > Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf  --max-missing 0.95 --recode --out Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj

#Downsample SNPs for Structure analysis as enough information in 10% of the loci
#(and more not informative because of linkage). In certain cases, when small number of markers detected,
#this step is unnecessary and all can be retained.
/home/sobczm/bin/vcflib/bin/vcfrandomsample \
--rate 0.1 Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.recode.vcf > Hg199_contigs_unmasked_filtered_noNMaj_subsampled.vcf
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
##Part of the path is missing here: should it be
cd ../
input_file=$input/Hg199_contigs_unmasked_filtered_noNMaj_subsampled.vcf
#Prepare population definition file. Each individual = new population
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
#dos2unix $outfile
#Number of loci in data file. Needed to edit the mainparams file
less Hg199_contigs_unmasked_filtered_noNMaj_subsampled.vcf | grep  -c -v "^#"
#Move file to the right folder
mv Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc structure_analysis
cd structure_analysis
#Note - Iteractions can be changed in the mainparams file.
#Important note - I used BURNIN 1000 NUMREPS 10000 for the first time. This gives large variance in lnPD, inconclusive run
#Minimum Burnin reps must be 100000. Number of reps must be between 10000 and 1000000. This will require days.
#mainparams file configuration: Label 1 (required), Popdata 1 (required but uninformative), Markernames 1 (required), Onerowperind 1 (possibly useless).
#The top row of the data file contains a list of L names corresponding to the markers used.
#Marker row and ind rows are separated in several row. This is possibly due to the large number of SNPs but I think this is consireded as a one single column for structure.
# BURNIN 100000 NUMREPS 10000, with K from 1 to 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 1 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 2 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 3 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 4 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 5 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 6 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 7 5
qsub -R y $scripts/execute_structure2.sh Hg199_contigs_unmasked_filtered_noNMaj_subsampled.struc 1 8 5

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
mkdir structureHarvester2
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester2
done

# structureHarvester - summarise the results
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=structureHarvester2 --out=structureHarvester2 --evanno --clumpp

# CLUMPP - permute the results
cd structure_analysis/structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
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
c=27
r=5
s=1
f=8
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done

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
n=27
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

-----------------------------------

The R0905 genome was used as reference for SNP calling since it is the most contiguous and complete genome.

#!/bin/bash
input=analysis/structure/R0905
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
#Move the file to the right directory
cp analysis/popgen/SNP_calling_R0905/R0905_good_contigs_unmasked_filtered.recode_annotated.vcf analysis/structure/R0905
#Remove outgroup and low-coverage isolates. Create a cut-down VCF and filter it
vcflib=/home/sobczm/bin/vcflib/bin
cd $input
$vcflib/vcfremovesamples R0905_good_contigs_unmasked_filtered.recode_annotated.vcf NMaj Ag11_B > R0905_isolates_filtered.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf R0905_isolates_filtered.vcf  --max-missing 0.95 --recode --out R0905_isolates_filtered

#Downsample SNPs for Structure analysis as enough information in 10% of the loci
#(and more not informative because of linkage). In certain cases, when small number of markers detected,
#this step is unnecessary and all can be retained.
/home/sobczm/bin/vcflib/bin/vcfrandomsample \
--rate 0.1 R0905_isolates_filtered.recode.vcf > R0905_isolates_filtered.recode_subsampled.vcf
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
##Part of the path is missing here: should it be
cd ../../../
input_file=$input/R0905_isolates_filtered.recode_subsampled.vcf
#Prepare population definition file. Each individual = new population (Integer value per sample defining populations).
#Popdata column can be edited at this point or later editing the struc file (probably easier).
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
#dos2unix $outfile
#Move file to the right folder
mv R0905_isolates_filtered.recode_subsampled.struc $input
cd $input
#struc file can be edited adding PopData using nano. Give a integer value on the second column depending of the isolates origin, p. e. UK will be 1.
nano R0905_isolates_filtered.recode_subsampled.struc
#File with edited popdata column save as R0905_isolates_filtered.recode_subsampled_ori.struc
#Number of loci in data file. Needed to edit the mainparams file
less R0905_isolates_filtered.recode_subsampled.vcf | grep  -c -v "^#"
#Note - Iteractions can be changed in the mainparams file.
#Important note - I used BURNIN 1000 NUMREPS 10000 for the first time. This gives large variance in lnPD, inconclusive run
#Minimum Burnin reps must be 100000. Number of reps must be between 10000 and 1000000. This will require days.
#mainparams file configuration: Label 1 (required), Popdata 1 (required but uninformative), Markernames 1 (required), Onerowperind 1 (possibly useless).
#The top row of the data file contains a list of L names corresponding to the markers used.
#Marker row and ind rows are separated in several row. This is possibly due to the large number of SNPs but I think this is consireded as a one single column for structure.
# BURNIN 100000 NUMREPS 10000, with K from 1 to 8
#If PopData is used, manual suggest to use the number of sampling locations + 3 for the K, but 8 should be more than enough since some sampling locations are very close.
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 1 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 2 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 3 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 4 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 5 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 6 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 7 5
qsub -R y $scripts/execute_structure2.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 8 5

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
mkdir structureHarvester
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester
done

# structureHarvester - summarise the results
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=structureHarvester --out=structureHarvester --evanno --clumpp

# CLUMPP - permute the results
cd structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
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
c=26
r=5
s=1
f=8
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done

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


-----------------------------------

The R0905 genome was used as reference for SNP calling since it is the most contiguous and complete genome.
No popdata used

#!/bin/bash
input=analysis_vAG/structure/R0905_NoPop
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
pgdspid=/home/sobczm/bin/PGDSpider_2.1.0.3
#Move the file to the right directory
cp analysis_vAG/SNPs/SNP_calling_R0905/R0905_good_contigs_unmasked_filtered.recode_annotated.vcf $input
#Remove outgroup and low-coverage isolates. Create a cut-down VCF and filter it
vcflib=/home/sobczm/bin/vcflib/bin
cd $input
$vcflib/vcfremovesamples R0905_good_contigs_unmasked_filtered.recode_annotated.vcf NMaj Ag11_B > R0905_isolates_filtered.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf R0905_isolates_filtered.vcf  --max-missing 0.95 --recode --out R0905_isolates_filtered

After filtering, kept 626425 out of a possible 626425 Sites

#Downsample SNPs for Structure analysis as enough information in 10% of the loci
#(and more not informative because of linkage). In certain cases, when small number of markers detected,
#this step is unnecessary and all can be retained.
/home/sobczm/bin/vcflib/bin/vcfrandomsample \
--rate 0.1 R0905_isolates_filtered.recode.vcf > R0905_isolates_filtered.recode_subsampled.vcf
#Prepare STRUCTURE input (PGDSpider does not work when wrapped up in a bash script, grrr)
#haploid (for diploid change the conversion script to vcf_to_structure_diploid.spid)
#!!!! Need to change the path to file with population definitions !!!
##Part of the path is missing here: should it be
cd ../../../
input_file=analysis_vAG/structure/R0905_NoPop/R0905_isolates_filtered.recode_subsampled.vcf
#Prepare population definition file. Each individual = new population (Integer value per sample defining populations).
#Popdata column can be edited at this point or later editing the struc file (probably easier).
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i " " $i "_pop"}' >"${input_file%.vcf}.lst"
#Copy the configuration file and change the path to the population definition file.
#For haploid organisms:
config=vcf_to_structure_haploid_pop.spid
#For diploid organisms:
#config=vcf_to_structure_diploid_pop.spid
cp $pgdspid/$config ./
dir=$PWD
list_file=$(echo "${input_file%.vcf}.lst")
sed -i 's,^\(VCF_PARSER_POP_FILE_QUESTION=\).*,\1'"$list_file"',' vcf_to_structure_haploid_pop.spid

#Also, create a label reference file to be used in the final step by distruct to label indidviduals in the output
names="${input_file%.vcf}.label"
grep "#CHROM" $input_file | head -1 | awk '{for(i=10;i<=NF;++i)print $i }' >temp
nl temp | sed 's/^ *//' | sed 's/\t/ /g' >$names
filename=$(basename "$input_file")
outfile="${filename%.*}.struc"
#Execute VCF to .struc (input format for STRUCTURE) conversion
java -jar $pgdspid/PGDSpider2-cli.jar -inputfile $input_file \
-inputformat VCF -outputfile $outfile -outputformat STRUCTURE -spid vcf_to_structure_haploid_pop.spid
#dos2unix $outfile
#Move file to the right folder
mv R0905_isolates_filtered.recode_subsampled.struc $input
cd $input
#struc file can be edited adding PopData using nano. Give a integer value on the second column depending of the isolates origin, p. e. UK will be 1.
#nano R0905_isolates_filtered.recode_subsampled.struc
#File with edited popdata column save as R0905_isolates_filtered.recode_subsampled_ori.struc
#Number of loci in data file. Needed to edit the mainparams file
less R0905_isolates_filtered.recode_subsampled.vcf | grep  -c -v "^#"
#Note - Iteractions can be changed in the mainparams file.
#Important note - I used BURNIN 1000 NUMREPS 10000 for the first time. This gives large variance in lnPD, inconclusive run
#Minimum Burnin reps must be 100000. Number of reps must be between 10000 and 1000000. This will require days.
#mainparams file configuration: Label 1 (required), Popdata 1 (required but uninformative), Markernames 1 (required), Onerowperind 1 (possibly useless).
#The top row of the data file contains a list of L names corresponding to the markers used.
#Marker row and ind rows are separated in several row. This is possibly due to the large number of SNPs but I think this is consireded as a one single column for structure.
# BURNIN 100000 NUMREPS 10000, with K from 1 to 8
#If PopData is used, manual suggest to use the number of sampling locations + 3 for the K, but 8 should be more than enough since some sampling locations are very close.
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 1 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 2 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 3 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 4 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 5 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 6 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 7 5
qsub -R y $scripts/execute_structure.sh R0905_isolates_filtered.recode_subsampled_ori.struc 1 8 5

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
mkdir structureHarvester
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester
done

# structureHarvester - summarise the results
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=structureHarvester --out=structureHarvester --evanno --clumpp

# CLUMPP - permute the results
cd structureHarvester
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
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
c=26
r=5
s=1
f=8
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done

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
