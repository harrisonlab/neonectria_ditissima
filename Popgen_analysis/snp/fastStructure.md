#Evaluation of population structure using fastStructure

This uses only biallelic SNP sites

###Set up initial variables

```bash
input=/data/scratch/gomeza/analysis/fastStructure2
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```
##### Remove outgroup. Create a cut-down VCF and filter it

```bash
cd $input
mkdir NoNMaj_logistic/
cd NoNMaj_logistic/
cp ../../popgen/SNP_calling3/Hg199_contigs_unmasked_filtered.recode_annotated.vcf ./
vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples Hg199_contigs_unmasked_filtered.recode_annotated.vcf NMaj > Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf  --max-missing 0.95 --recode --out Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj
```

###Convert from VCF to Plink's PED format

```bash
input_file=Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.recode.vcf
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```

####Test various K values (K represents model complexity)

This will run with prior logistic

```bash
# Set minimum number of considered clusters
s=1
# Set maximum number of considered clusters
f=10
for i in $(seq $s $f)
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i
done
```

####Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
input_vcf_file=${input_file%.vcf}

python $structure/chooseK.py --input $input_vcf_file > ${input_file%.vcf}_K_choice
```

####Visualise expected admixture proportions with Distruct plots

This uses the mean of variational posterior distribution over admixture proportions

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# X11 forwarding is required, set up on an OSX local machine running OSX v10.13.4 using:
# https://stackoverflow.com/questions/39622173/cant-run-ssh-x-on-macos-sierra
# -X option needed for ssh to allow X11 forwarding. -Y option work for me.

input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
input_vcf_file=${input_file%.vcf}
Popfile=${input_file%.vcf}.lab
s=1
f=10
structure=/home/sobczm/bin/fastStructure
for i in $(seq $s $f)
do
    Output=${input_file%.vcf}_${i}.svg
    python $structure/distruct_mod.py -K $i --input $input_vcf_file --output $Output --title K$i --popfile $Popfile
done
```
Works, Throws error for being unable to load a Fontconfig file (non-critical)



### No prior option used

##### Set up initial variables

```bash
input=/data/scratch/gomeza/analysis/fastStructure2
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```

##### Remove outgroup. Create a cut-down VCF and filter it

```bash
cd $input
mkdir NoNMaj/
cd NoNMaj/
vcflib=/home/sobczm/bin/vcflib/bin
$vcflib/vcfremovesamples ../Hg199_contigs_unmasked_filtered.recode_annotated.vcf NMaj > Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf  --max-missing 0.95 --recode --out Hg199_filtered_recode_noNMaj
```

##### Convert from VCF to Plink's PED format

```bash
input_file=Hg199_filtered_recode_noNMaj.recode.vcf
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```

####Test various K values (K represents model complexity)

```bash
# Set minimum number of considered clusters
s=1
# Set maximum number of considered clusters
f=10
for i in $(seq $s $f)
do
qsub $scripts/sub_fast_structure_noprior.sh ${input_file%.vcf} $i
done
```

####Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
input_file=Hg199_filtered_recode_noNMaj.recode.vcf
input_vcf_file=${input_file%.vcf}

python $structure/chooseK.py --input $input_vcf_file > ${input_file%.vcf}_K_choice
```

####Visualise expected admixture proportions with Distruct plots

This uses the mean of variational posterior distribution over admixture proportions

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# X11 forwarding is required, set up on an OSX local machine running OSX v10.13.4 using:
# https://stackoverflow.com/questions/39622173/cant-run-ssh-x-on-macos-sierra
# -X option needed for ssh to allow X11 forwarding. -Y option works for me.

input_file=Hg199_filtered_recode_noNMaj.recode.vcf
input_vcf_file=${input_file%.vcf}
Popfile=${input_file%.vcf}.lab
s=1
f=10
structure=/home/sobczm/bin/fastStructure
for i in $(seq $s $f)
do
    Output=${input_file%.vcf}_${i}.svg
    python $structure/distruct_mod.py -K $i --input $input_vcf_file --output $Output --title K$i --popfile $Popfile
done
```
Works, Throws error for being unable to load a Fontconfig file (non-critical)


## Logistic prior option

###Convert from VCF to Plink's PED format

```bash
cd $input
mkdir -p prior_logistic
cp analysis/popgen/SNP_calling/N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf prior_logistic/
cd prior_logistic
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
#plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```
####Test various K values (K represents model complexity)

By default, fastStructure runs using F-prior. Alternatively, logistic prior can be used when structure is difficult to resolve.

```bash
# Set minimum number of considered clusters
s=1
# Set maximum number of considered clusters
f=10
for i in $(seq $s $f)
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i --prior logistic
done
```

For some reason K3 is taking too long to run, so I run the following commands in blacklace11.

input=${input_file%.vcf}
structure=/home/sobczm/bin/fastStructure
python $structure/structure.py -K 3 --input $input --output $input --prior logistic

This is running on screen in blacklace11. the next is running from test folder.


```bash
mkdir prior_logistic_all
cp N.ditissima_contigs_unmasked_filtered.recode_annotated.* prior_logistic_all/
#I run K3 in different nodes at the same time since the first run was taking too long. The final results are copied together into the prior_logistic_all folder.
mv N.ditissima_contigs_unmasked_filtered.recode_annotated.3* prior_logistic_all/
```

####Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
cd prior_logistic_all/
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf  
input_vcf_file=${input_file%.vcf}

python $structure/chooseK.py --input $input_vcf_file > ${input_file%.vcf}_K_choice
```

####Visualise expected admixture proportions with Distruct plots

This uses the mean of variational posterior distribution over admixture proportions

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab
#Draw output
# X11 forwarding is required, set up on an OSX local machine running OSX v10.13.4 using:
# https://stackoverflow.com/questions/39622173/cant-run-ssh-x-on-macos-sierra
# -X option needed for ssh to allow X11 forwarding. -Y option work for me.
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
input_vcf_file=${input_file%.vcf}
Popfile=${input_file%.vcf}.lab
s=1
f=10
structure=/home/sobczm/bin/fastStructure
for i in $(seq $s $f)
do
    Output=${input_file%.vcf}_${i}.svg
    python $structure/distruct_mod.py -K $i --input $input_vcf_file --output $Output --title K$i --popfile $Popfile
done
```

##Create a cut-down vcf that only include strains of interest.

This commands can be used to remove samples from our VCF files. Since I had clearly a different population corresponding to the isolate R68-17, I remove this one and then I repeated the fast_structure analysis with the logistic prior.

```bash
scripts=/home/sobczm/bin/popgen/snp
vcflib=/home/sobczm/bin/vcflib/bin

$vcflib/vcfremovesamples N.ditissima_contigs_unmasked_filtered.vcf  R68-17 > N.ditissima_contigs_unmasked_filtered_noR68.vcf

input_file=N.ditissima_contigs_unmasked_filtered_noR68.vcf
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log

#I remove a few samples from different origins for comparison
$vcflib/vcfremovesamples N.ditissima_contigs_unmasked_filtered.vcf  R68-17 ND8 ND9 P112 BGV304 OPC304 > N.ditissima_contigs_unmasked_filtered_noEsBrIt.vcf
```

###Convert from VCF to Plink's PED format

```bash
cd $input/prior_logistic
mkdir -p NoR68
cp N.ditissima_contigs_unmasked_filtered_noR68.vcf NoR68/
cd NoR68
input_file=N.ditissima_contigs_unmasked_filtered_noR68.vcf
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```

####Test various K values (K represents model complexity)

```bash
# Set minimum number of considered clusters
s=1
# Set maximum number of considered clusters
f=10
for i in $(seq $s $f)
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i
done
```

####Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
input_file=N.ditissima_contigs_unmasked_filtered_noR68.vcf
input_vcf_file=${input_file%.vcf}

python $structure/chooseK.py --input $input_vcf_file > ${input_file%.vcf}_K_choice
```

####Visualise expected admixture proportions with Distruct plots

This uses the mean of variational posterior distribution over admixture proportions

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# X11 forwarding is required, set up on an OSX local machine running OSX v10.13.4 using:
# https://stackoverflow.com/questions/39622173/cant-run-ssh-x-on-macos-sierra
# -X option needed for ssh to allow X11 forwarding. -Y option work for me.

input_file=N.ditissima_contigs_unmasked_filtered_noR68.vcf
input_vcf_file=${input_file%.vcf}
Popfile=${input_file%.vcf}.lab
s=1
f=10
structure=/home/sobczm/bin/fastStructure
for i in $(seq $s $f)
do
    Output=${input_file%.vcf}_${i}.svg
    python $structure/distruct_mod.py -K $i --input $input_vcf_file --output $Output --title K$i --popfile $Popfile
done
```


###Convert from VCF to Plink's PED format

```bash
mkdir -p $input
cd $input
cp ../popgen/SNP_calling3/Hg199_contigs_unmasked_filtered.recode_annotated.vcf  .
input_file=Hg199_contigs_unmasked_filtered.recode_annotated.vcf
#plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```


#Evaluation of population structure using fastStructure

R0905 genome used as reference. This uses only biallelic SNP sites

###Set up initial variables

```bash
input=/data/scratch/gomeza/analysis/fastStructure/R0905
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
```
##### Remove outgroup. Create a cut-down VCF and filter it

```bash
cd $input
mkdir logistic_prior/
mkdir default
cd default/
#Outgroup and low coverage sequenced isolates were removed and non informative SNPs were filtered. This was done in the structure script.
cp ../../../structure/R0905/R0905_isolates_filtered.recode_subsampled.vcf ./

#vcflib=/home/sobczm/bin/vcflib/bin
#$vcflib/vcfremovesamples Hg199_contigs_unmasked_filtered.recode_annotated.vcf NMaj > Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf
#vcftools=/home/sobczm/bin/vcftools/bin
#$vcftools/vcftools --vcf Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj.vcf  --max-missing 0.95 --recode --out Hg199_contigs_unmasked_filtered.recode_annotated_noNMaj
```

###Convert from VCF to Plink's PED format

```bash
input_file=R0905_isolates_filtered.recode_subsampled.vcf
plink --indep-pairwise 100000 1 0.5 --allow-extra-chr --const-fid 0 --vcf $input_file --make-bed --recode --out ${input_file%.vcf} > ${input_file%.vcf}.log
```

####Test various K values (K represents model complexity)

This will run with prior logistic

```bash
# Set minimum number of considered clusters
s=1
# Set maximum number of considered clusters
f=10
for i in $(seq $s $f)
do
qsub $scripts/sub_fast_structure.sh ${input_file%.vcf} $i
done
```

####Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
input_vcf_file=${input_file%.vcf}

python $structure/chooseK.py --input $input_vcf_file > ${input_file%.vcf}_K_choice
```
