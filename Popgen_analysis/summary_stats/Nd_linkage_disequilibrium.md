# Commands to run analysis of linkage disequilibrium

```bash
input=/data/scratch/gomeza/analysis/summary_stats
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/summary_stats
vcftools=/home/sobczm/bin/vcftools/bin
```

Calculate D, D' and r^2 for SNPs separated by between 1 and 100 kbp
in non-pathogens (program calculates the stats using only the individuals
listed after "--indv" switch, if that option not present, all the individuals in the file)
In order to calculate all versus all SNP comparison, remove options --ld-window-bp-min
and --ld-window-bp

## Analysis on all N. ditissima isolates less R68-17 based on fastStructure results

```bash
$vcftools/vcftools --vcf N.ditissima_contigs_unmasked_noR68_filtered.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv Ag02 --indv Ag04 --indv Ag05 --indv Ag06 --indv Ag08 --indv Ag09_A --indv Ag11_A --indv Ag11_B --indv Ag11_C --indv BGV344 --indv Hg199 --indv ND8 --indv ND9 --indv P112 --indv OPC304 --indv R0905 --indv R37-15 --indv R39-15 --indv R41-15 --indv R42-15 --indv R45-15 --indv R68-17 --indv R6-17-2 --indv R6-17-3
mv out.hap.ld ld.Nd

qsub $scripts/sub_plot_ld.sh ld.Nd

mkdir -p NdLd
mv ld* NdLd/.
```

### LD plot (heatmap) for r2 values per contig
```bash
qsub $scripts/sub_ld_plot.sh ld.Nd
```

### Calculate LD decay and plot curve

```
Size is:
The number of individuals sampled * number of chromosomes sequenced * ploidy
I have three as Pf n = 10 - 12 (https://doi.org/10.1073/pnas.96.10.5878)
```

```bash
for Size in Nd8 Nd9
do
    Out_File_Fitted=UK123/r^2_decay_"$Size"_fitted.pdf
    Out_File_Unfitted=UK123/r^2_decay_"$Size"_unfitted.pdf
    LD_file=UK123/ld.UK123_no_min
    units=bp
    window_size=100000
    bin_size=1000
    Cstart=0.1
    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
    echo "Sample size of $Size:"
    Rscript --vanilla $ProgDir/plot_LD_decay.R --out_file_fitted $Out_File_Fitted --out_file_unfitted $Out_File_Unfitted --Chromosome_number $Size --LD_statistics $LD_file --units $units --window_size $window_size --bin_size $bin_size --Cstart $Cstart
    printf "\n"
done
```
