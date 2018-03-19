PIPELINE for Mauve whole genome alignments and subsequent CLONALFRAME analysis of regions of recombination

# Use move_contigs to order genomes based on reference(Use ssh -Y cluster to ensure trusted X11 connection)

```bash
for GENOME in analysis/genome_alignment/mauve/N.ditissima/RS305p/LDPK01.1.fsa_nt.fasta  ; do
GENOME_SHORT=$(basename $GENOME | sed s/".fasta"//g )
echo $GENOME_SHORT
mkdir "$GENOME_SHORT"_dir
REFERENCE=analysis/genome_alignment/mauve/N.ditissima/RS305p/N.ditissima_contigs_unmasked.fa
echo $REFERENCE

java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME
done

java -Xmx500m -cp /home/hulinm/local/src/mauve_2.3.1/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME

for GENOME in analysis/genome_alignment/mauve/N.ditissima/RS324p/LDPL01.1.fsa_nt.fasta  ; do
GENOME_SHORT=$(basename $GENOME | sed s/".fasta"//g )
echo $GENOME_SHORT
mkdir "$GENOME_SHORT"_dir
REFERENCE=analysis/genome_alignment/mauve/N.ditissima/RS324p/N.ditissima_contigs_unmasked.fa
echo $REFERENCE

java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME
done

java -Xmx500m -cp /home/hulinm/local/src/mauve_2.3.1/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS324p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME
```

```bash
progressiveMauve --output=RS305p.xmfa N.ditissima_contigs_unmasked.fa LDPK01.1.fsa_nt.fasta
java -cp Mauve.jar org.gel.mauve.analysis.SnpExporter -f RS305p.xmfa -o RS305p.snps

java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output analysis/genome_alignment/mauve/N.ditissima/RS305p/"$GENOME_SHORT"_dir -ref $REFERENCE -draft $GENOME
```
