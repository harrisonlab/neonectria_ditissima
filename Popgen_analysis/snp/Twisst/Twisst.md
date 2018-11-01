
General tools for genomic analyses. https://github.com/simonhmartin/genomics_general

#Twisst

Topology weightying by iterative sampling of sub-trees.

Literature: Martin and Van Belleghem 2017 and Van Belleghem et al. 2017.

Firstly, I used a tree based on all the SNPs in newick format.

```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t Hg199_contigs_unmasked_filtered_nj.nwk -w output.weights.csv.gz -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 -g NM NMaj --method complete
```

##Pipeline to generate the input trees file

##Parsing VCF files

Transform simple vcf file to a file containing scaffold, position and genotype for each sample.

```bash
python $scripts/parseVCF.py -i Hg199_contigs_unmasked_filtered.vcf --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.gz
```
If you have whole genome sequence data, it is recommended to infer trees for narrow genomic intervals. 50 SNPs proved a useful window size in various simulations.

```bash
python $scripts/phyml_sliding_windows.py -T 10 -g output.geno.gz --prefix output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n
```

##Run Twisst

```bash
screen -a
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t output.phyml_bionj.w50.trees.gz -w output.weights2.csv --outputTopos topologies.trees -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 -g NM NMaj --method complete
```

#Analysis of individuals Contigs

Contig 1

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst

python $scripts/parseVCF.py -i Hg199_contigs_unmasked_filtered.vcf --include contig_1 --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.contig1.gz

python $scripts/phyml_sliding_windows.py -T 10 -g output.geno.contig1.gz --prefix output.phyml_bionj.w50.contig1 -w 50 --windType sites --model GTR --optimise n

screen -a
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t output.phyml_bionj.w50.contig1.trees.gz -w output.weights.contig1.csv --outputTopos topologies.trees.contig1 -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 -g NM NMaj --method complete

Contig 2

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst

python $scripts/parseVCF.py -i Hg199_contigs_unmasked_filtered.vcf --include contig_2 --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > contig_2/output.geno.contig2.gz

python $scripts/phyml_sliding_windows.py -T 10 -g contig_2/output.geno.contig2.gz --prefix contig_2/output.phyml_bionj.w50.contig2 -w 50 --windType sites --model GTR --optimise n

screen -a
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t contig_2/output.phyml_bionj.w50.contig2.trees.gz -w contig_2/output.weights.contig2.csv --outputTopos contig_2/topologies.trees -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 -g NM NMaj --method complete

Contig 3

scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst

python $scripts/parseVCF.py -i Hg199_contigs_unmasked_filtered.vcf --include contig_3 --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > contig_3/output.geno.gz

python $scripts/phyml_sliding_windows.py -T 10 -g contig_3/output.geno.gz --prefix contig_3/output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n

screen -a
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t contig_3/output.phyml_bionj.w50.trees.gz -w contig_3/output.weights.csv --outputTopos contig_3/topologies.trees -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 -g NM NMaj --method complete





4 taxons.

```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t Hg199_contigs_unmasked_filtered_nj.nwk -w output.weights.csv.gz -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 --method complete
```


```bash
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
OutDir=contig_13

python $scripts/parseVCF.py -i Hg199_contigs_unmasked_filtered.vcf --include contig_13,contig_14,contig_15,contig_16,contig_17 --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > $OutDir/output.geno.gz

python $scripts/phyml_sliding_windows.py -T 10 -g $OutDir/output.geno.gz --prefix $OutDir/output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n

screen -a

OutDir=contig_13
scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/Twisst
python $scripts/twisst.py -t $OutDir/output.phyml_bionj.w50.trees.gz -w $OutDir/output.weights.csv --outputTopos $OutDir/topologies.trees -g N Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3,Ag08,Ag09_A,Ag11_A,Ag11_B,Ag11_C,R37-15,R39-15,R41-15,R42-15,R45-15 -g W ND8,ND9 -g S P112,OPC304,BGV344 -g E R68-17-C2,R68-17-C3,SVK1,SVK2 --method complete
```
