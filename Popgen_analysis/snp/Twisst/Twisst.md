
General tools for genomic analyses. https://github.com/simonhmartin/genomics_general

#Twisst

Topology weightying by iterative sampling of sub-trees.

Literature: Martin and Van Belleghem 2017 and Van Belleghem et al. 2017.

Firstly, I used a tree based on all the SNPs in newick format.

```bash
python twisst.py -t N.ditissima_contigs_unmasked_filtered_nj.nwk -w output.weights.csv.gz -g ENG Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3 -g NI Ag11_A,Ag11_B,Ag11_C -g IRL Ag08,Ag09_A -g BR ND8,ND9 -g BE R37-15,R39-15 -g NL R41-15,R42-15,R45-15 -g ES P112,OPC304,BGV344 -g IT R68-17 --method complete
```

##Pipeline to generate the input trees file

##Parsing VCF files

Transform simple vcf file to a file containing scaffold, position and genotype for each sample.

```bash
python /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/parseVCF.py -i N.ditissima_contigs_unmasked_filtered.vcf --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.gz
```
If you have whole genome sequence data, it is recommended to infer trees for narrow genomic intervals. 50 SNPs proved a useful window size in various simulations.

```bash
python /home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp/phyml_sliding_windows.py -T 10 -g output.geno --prefix output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n
```

Run Twisst

```bash
screen -a
python twisst.py -t output.phyml_bionj.w50.trees.gz -w output.weights2.csv.gz -o topologies.trees -g ENG Ag02,Ag04,Ag05,Ag06,Hg199,R0905,R6-17-2,R6-17-3 -g NI Ag11_A,Ag11_B,Ag11_C -g IRL Ag08,Ag09_A -g BR ND8,ND9 -g BE R37-15,R39-15 -g NL R41-15,R42-15,R45-15 -g ES P112,OPC304,BGV344 -g IT R68-17 --method complete
```
