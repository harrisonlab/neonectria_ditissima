Identification of potential contaminants in the raw reads

```bash
for Strain in Ag11_B R68-17; do
Fusarium=$(ls ../fusarium/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa)
Stenotrophomonas=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/stenotrophomonas/all_complete_GbStenotrophomonas.fasta)
for Reference in $(ls $Fusarium $Stenotrophomonas); do
for StrainPath in $(ls -d qc_dna/paired/*/$Strain); do
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix=$(basename $Reference | sed 's/.fasta//g' | sed 's/.fa//g' | sed 's/.fna//g')
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_${Prefix}
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
done
done
```
