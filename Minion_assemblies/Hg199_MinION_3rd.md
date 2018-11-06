# Hg199 minion assembly using miniasm

Thisis fast approach for assembling and correcting PacBio and MinION data using miniasm and racon.

## Assembly

### Assembly will with full trimming of reads:

The next 2 steps were done in the first assembly of the minion data, but repeated here.

#### Splitting reads and trimming adapters using porechop

```bash
    RawReads=raw_dna/minion/N.ditissima/Hg199/*allfiles.fastq.gz
    OutDir=qc_dna/minion/N.ditissima/Hg199_2ndround
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
```

#### Read correction using Canu

```bash
  for TrimReads in $(ls qc_dna/minion/N.ditissima/Hg199_2ndround/*allfiles_trim.fastq.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu_minion2/N.ditissima/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub -R y $ProgDir/sub_canu_correction.sh $TrimReads 46m $Strain $OutDir
  done
```

#### Assembly

```bash
# If reads have same names or same part splited by space, fix them using rename.sh from bbtools.
# This can only be done in blacklace11, since it needs a specific library for blasr.
ssh blacklace11.blacklace
/home/gomeza/prog/bbtools/bbmap/rename.sh in=Hg199.correctedReads.fasta.gz out=CorrecteOut.fasta prefix=Hg199

# Fast all-against-all overlap of raw reads
# Overlap for MinION reads (or use "-x ava-pb" for Pacbio read overlapping)
/home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 CorrecteOut.fasta CorrecteOut.fasta | gzip -1 > Hg199_fastq_allfiles.paf.gz

# Concatenate pieces of read sequences to generate the final sequences.
# Thus the per-base error rate is similar to the raw input reads. Make sure you correct your reads.
# Layout
/home/gomeza/prog/miniasm/miniasm -f CorrecteOut.fasta Hg199_fastq_allfiles.paf.gz > reads.gfa

# Convert gfa file to fasta file.
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > Hg199_miniasm.fa
```

## Error correction using racon:

```bash
  for Assembly in $(ls Hg199_miniasm2/Hg199_miniasm.fa); do
    Strain=Hg199
    Organism=N.ditissima
    echo "$Organism - $Strain"
    ReadsFq=$(ls qc_dna/minion/N.ditissima/Hg199/*allfiles_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/racon
    qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
```
