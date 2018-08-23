# Further analysis of linkage disequilibrium using the LDhat suite

## Count number of contigs to be used for LDhat

```bash
for num in 1
do
  input=analysis/LDhat/Hg199
  mkdir -p $input
  cd $input
  #This needs to be done with the newest genome. But I used this one for the snp calling.
  Input_file=../../../Hg199_genome/repeat_masked/N.ditissima/Hg199_minion/N.ditissima_contigs_unmasked.fa
  Output_file=cut_down_assembly.fa
  ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools
  Minimum_contig_size=500000
  python $ProgDir/Size_selection_assembly.py --input $Input_file --output $Output_file --min_size $Minimum_contig_size
  echo "The number of contigs larger than 500Kb is:"
  cat $Output_file | grep '>' | wc -l
done
```

```
The number of contigs larger than 500kb is:
30
```

## Convert vcf file to correct format for LDhat

Treat each contig as a separate "chromosome"

```bash
vcftools=/home/sobczm/bin/vcftools/bin
input_vcf=../../popgen/SNP_calling/N.ditissima_contigs_unmasked_filtered.recode_annotated.vcf
Out_prefix=merged_contigs_unmasked_Hg199
cut_down_assembly=cut_down_assembly.fa

contigs=$(cat $cut_down_assembly | grep '>' | tr -d '>' | cut -f2 -d "_")
for num in $(echo $contigs)
do
  contig_name=$(echo contig_"$num")
  Out_prefix=$(echo ldhat_"$contig_name")
  mkdir -p $contig_name
  $vcftools/vcftools --vcf $input_vcf --out $Out_prefix --chr $contig_name --ldhat
  mv ldhat_"$contig_name".log $contig_name/.
  mv ldhat_"$contig_name".ldhat.locs $contig_name/.
  mv ldhat_"$contig_name".ldhat.sites $contig_name/.
done
```

## Use pairwise to build a lookup table and calculate some statistics

The following commands must be run in a screen session running a qlogin job
Requires some user input on the command line

My recommendations:

```
Unless another source of evidence available, use the suggested value of theta
Minor changes in theta do not seem to have an effect
4Ner should range between 20 - 100
No. of points on grid should range between 21 - 201
Larger values of 4Ner & No. of points take longer, but are more accurate
I recommend number of points to be twice the max 4Ner + 1
(eg. 4Ner max = 200, number of points = 401)

After generation of the table, use default grid value for recombination rate
If estimates at extreme of grid, change the defaults - can exceed above limits:
'Number of classes to estimate for' = No. of points

Ignore sliding window analyses - interval is better

Use option 2 for rmin to write output to a file not just screen

Allow estimation of 4Ner by moment method

Test for recombination
```

The main thing needed from pairwise is the lookup table
The rest are mostly improved upon by interval & LDhot

```bash
for input_dir in $(ls -d contig_*)
do
    cd $input_dir
    seq_file=ldhat_"$input_dir".ldhat.sites
    loc_file=ldhat_"$input_dir".ldhat.locs
    echo "$input_dir"
    pairwise -seq $seq_file -loc $loc_file
    cd ../
done
```

Iterative process, follow the ls -d with grep -v -w 'contig_id' to narrow selection
Results can vary due to the monte carlo simulation being used

## Use complete to refine the lookup table

Use the same 4Ner_max, no of points and theta as in pairwise, stored in log file
Done separately for each contig

```bash
for input_dir in $(ls -d contig_*)
do
    n=$(cat $input_dir/ldhat_contig_*.ldhat.sites | grep '>' | wc -l)
    four_Ner_max=$(cat $input_dir/outfile.txt | grep -e 'Maximum at 4Ner(region)' | cut -f1 -d ':' | cut -f2 -d '=')
    no_points=201
    theta=$(cat $input_dir/outfile.txt | grep -e 'Theta' | cut -f2 -d '=')
    ProgDir=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis
    Jobs=$(qstat | grep 'sub_comple' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_comple' | grep 'qw' | wc -l)
    done
    printf "\n"
    qsub $ProgDir/sub_complete.sh $input_dir $n $four_Ner_max $no_points $theta
done
```

## Run interval to estimate the recombination rate

This only uses the crossing over model with a bayesian rjMCMC approach

```bash
for input_dir in $(ls -d contig_*)
do
    mv $input_dir/new_lk.txt $input_dir/exhaustive_lk.txt
    sequence_data=ldhat_"$input_dir".ldhat.sites
    location_data=ldhat_"$input_dir".ldhat.locs
    lookup_table=exhaustive_lk.txt
    #10,000,000 is the recommended number of iterations in the manual
    iterations=10000000
    # Sample rate keeps every xth result of the markov chain
    sample_rate=2000
    # Values from 0 - 50 are recommended, 5 is okay for humans
    # If you have data on expected recombination rate
    # run simulations with different penalities to test which is least 'noisy'
    block_penalty=5
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    Jobs=$(qstat | grep 'sub_interv' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_interv' | grep 'qw' | wc -l)
    done
    printf "\n"
    qsub $ProgDir/sub_interval.sh $input_dir $sequence_data $location_data $lookup_table $iterations $block_penalty $sample_rate
done
```

### Visualise output of interval using stat

Average values are stored in the output log files from sge
Ensure all the jobs in one loop are finished before starting the other
They both produce intermediate files of the same name initially due to source code

```bash
for input_dir in $(ls -d contig_*)
do
    cd $input_dir
    rates_file=rates.txt
    # Specify number of iterations to discard
    # 100,000 iterations is recommended min
    # So divide 100,000 by sampling rate for burn-in value
    burn_in=50
    location_file=ldhat_"$input_dir".ldhat.locs
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    Jobs=$(qstat | grep 'sub_stat' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_stat' | grep 'qw' | wc -l)
    done
    qsub $ProgDir/sub_stat_rates.sh $rates_file $burn_in $location_file
    cd ../
done

for input_dir in $(ls -d contig_*)
do
    cd $input_dir
    bounds_file=bounds.txt
    # Specify number of iterations to discard
    # 100,000 iterations is recommended min
    # So divide 100,000 by sampling rate for burn-in value
    burn_in=50
    location_file=ldhat_"$input_dir".ldhat.locs
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    Jobs=$(qstat | grep 'sub_stat' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_stat' | grep 'qw' | wc -l)
    done
    qsub $ProgDir/sub_stat_bounds.sh $bounds_file $burn_in $location_file
    cd ../
done
```

Whilst rhomap is included in LDhat for identifying hotspots, this uses a model
based on human recombination. LDhat authors recommend using LDhot.

## Investigate recombination hotspots using LDhot

### Installation and profile modifications

Not required for other users

```bash
cd  /home/adamst/prog
git clone https://github.com/auton1/LDhot.git

cd LDhot
make
```

### To use LDhot, add the following line to your profile

```bash
PATH=${PATH}:/home/adamst/prog/LDhot
```

This uses the res_rates file from LDhat program stat to identify hotspots
This program requires a minimum of 100 SNPs per contig
A failure will be recorded in the output log file
If running in screen, first run:

```bash
. ~/.profile
```

## Initial identification of recombination hotspots

```bash
for input_dir in $(ls -d contig_*)
do
    cd $input_dir
    sequence_file=ldhat_"$input_dir".ldhat.sites
    location_file=ldhat_"$input_dir".ldhat.locs
    lookup_table=exhaustive_lk.txt
    rates_file=res_rates.txt
    # 1000 is the recommended minimum value for number of simulations
    num_simulations=1000
    Out_prefix=LDhot_results
    Jobs=$(qstat | grep 'sub_ldhot' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_ldhot' | grep 'qw' | wc -l)
    done
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    qsub $ProgDir/sub_ldhot.sh $sequence_file $location_file $lookup_table $rates_file $num_simulations $Out_prefix
    cd ../
done
```

```
LDhot ran on the following contigs:
contig 8
```

### Summarise LDhot results

```bash
for input_dir in $(ls -d contig_* | grep -e 'contig_8')
do
    cd $input_dir
    rates_file=res_rates.txt
    hotspot_file=LDhot_results.hotspots.txt
    Out_prefix=LDhot_summary
    Jobs=$(qstat | grep 'sub_ldhot_' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_ldhot_' | grep 'qw' | wc -l)
    done
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    qsub $ProgDir/sub_ldhot_summary.sh $rates_file $hotspot_file $Out_prefix
    cd ../
done
```

Location results are in kilobases

## Plot results of interval and LDhot as a line graph

Plot recombination rate over the contig
With hotspots highlighted (where identified)

```bash
for input_dir in $(ls -d contig_*)
do
    Rho_file=$input_dir/res_rates.txt
    Out_file=$input_dir/Rho_plot.pdf
    Hotspot_file=$input_dir/LDhot_summary.hot_summary.txt
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    if [ -f $Hotspot_file ]
    then
        Rscript --vanilla $ProgDir/plot_LD.R --out_file $Out_file --res_in $Rho_file --hotspot_in $Hotspot_file
    else
        Rscript --vanilla $ProgDir/plot_LD_no_hotspot.R --out_file $Out_file --res_in $Rho_file
    fi
done
```
