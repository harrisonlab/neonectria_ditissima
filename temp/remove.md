```bash
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.csv
# printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly" > tmp.csv
printf \
"Trim:
Sequence name,	length,	span(s),	apparent source
contig_2	3685016	2630642..2631051
contig_7	2122932	14290..15776
contig_15	1318801	1303809..1306436
contig_17	1243348	1098076..1099539
contig_21	602725	120117..121484,162018..164688,289092..289649
" \t
> tmp.csv
for Assembly in $(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/*_contigs_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
mkdir -p $OutDir
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
done
rm tmp.csv
```

printf \
"Trim:
Sequence name,\tlength,\tspan(s),\tapparent source
contig_2\t3685016\t2630642..2631051\tlow coverage
contig_7\t2122932\t14290..15776\tlow coverage
contig_15\t1318801\t1303809..1306436\tlow coverage
contig_17\t1243348\t1098076..1099539\tlow coverage
contig_21\t602725\t120117..121484,162018..164688,289092..289649\tlow coverage
"
> tmp.csv

```bash
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.csv
# printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly" > tmp.csv
printf \
"Trim:
Sequence name,\tlength,\tspan(s),\tapparent source
contig_2\t3685016\t2630642..2631051\tlow coverage
contig_7\t2122932\t14290..15776\tlow coverage
contig_15\t1318801\t1303809..1306436\tlow coverage
contig_17\t1243348\t1098076..1099539\tlow coverage
contig_21\t602725\t120117..121484,162018..164688,289092..289649\tlow coverage
"
> tmp.csv
for Assembly in $(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/*_contigs_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
mkdir -p $OutDir
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
done
rm tmp.csv
```


```bash
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
touch tmp.csv
# printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly" > tmp.csv
printf \
"Trim:
Sequence name,\tlength,\tspan(s),\tapparent source
contig_2\t3685016\t2630642..2631051\tlow coverage
contig_7\t2122932\t14290..15776\tlow coverage
contig_15\t1318801\t1303809..1306436\tlow coverage
contig_17\t1243348\t1098076..1099539\tlow coverage
contig_21\t602725\t120117..121484,162018..164688,289092..289649\tlow coverage
" \
> tmp.csv
for Assembly in $(ls assembly/merged_canu_spades/*/R0905/filtered_contigs/*_contigs_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
mkdir -p $OutDir
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
done
rm tmp.csv
```
