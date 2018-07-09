#Evaluation of population structure using fastStructure

This uses only biallelic SNP sites

###Set up initial variables

```bash
input=/data/scratch/gomeza/fastStructure
scripts=/home/sobczm/bin/popgen/snp
```

###Convert from VCF to Plink's PED format

```bash
mkdir -p $input
cd $input
cp analysis/popgen/SNP_calling/N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf .
input_file=N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
#plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log
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
Works, Throws error for being unable to load a Fontconfig file (non-critical)
```
