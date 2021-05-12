# Quickmerge

## Merging long and short read assemblies

```bash
# Best Hg199 assembly
    for MinIONAssembly in $(ls assembly_VP/flye/N.ditissima/Hg199/racon_10/medaka/medaka/pilon/pilon10_renamed.fasta); do
        Organism=$(echo $MinIONAssembly | rev | cut -f7 -d '/' | rev)
        Strain=$(echo $MinIONAssembly | rev | cut -f6 -d '/' | rev)
        Assembler=$(echo $MinIONAssembly | rev | cut -f8 -d '/' | rev)
        HybridAssembly=$(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/Hg199_hybrid_renamed.fasta)
        AnchorLength=384699
        OutDir=assembly_VP/merged_assemblies/"$Assembler"_spades/$Organism/"$Strain"_minion_380k
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    done
```

```bash
# Best R0905 assembly
    for MinIONAssembly in $(ls assembly_VP/canu/N.ditissima/R0905/medaka/medaka/pilon/pilon10_renamed.fasta); do
        Organism=$(echo $MinIONAssembly | rev | cut -f6 -d '/' | rev)
        Strain=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
        Assembler=$(echo $MinIONAssembly | rev | cut -f7 -d '/' | rev)
        HybridAssembly=$(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/R0905_hybrid_renamed.fasta)
        AnchorLength=524431
        OutDir=assembly_VP/merged_assemblies/"$Assembler"_spades/$Organism/"$Strain"_minion_520k
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/quickmerge.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    done
```

## Polish merged genomes

```bash
for Assembly in $(ls assembly_VP/merged_assemblies/canu_spades/N.ditissima/R0905_minion_520k/merged_pacbio_hybrid_merged.fasta); do
    Strain=R0905
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done

for Assembly in $(ls assembly_VP/merged_assemblies/flye_spades/N.ditissima/Hg199_minion_380k/merged_pacbio_hybrid_merged.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

## Quast and BUSCO

```bash
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
  for Assembly in $(ls assembly_VP/merged_assemblies/*/N.ditissima/*/pilon/pilon_10.fasta); do
      OutDir=$(dirname $Assembly)
      sbatch $ProgDir/quast.sh $Assembly $OutDir
  done
```

```bash
    for Assembly in $(ls assembly_VP/merged_assemblies/*/N.ditissima/*/pilon/pilon_10.fasta); do
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
    sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done
```

## Hybridgenomes polished

```bash
# These genomes were polished by mistake. They might be used for quickmerge.

for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/R0905_hybrid_renamed.fasta); do
    Strain=R0905
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/R0905_all)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done

for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/Hg199_hybrid_renamed.fasta); do
    Strain=Hg199
    Organism=N.ditissima
    IlluminaDir=$(ls -d qc_dna/paired/N.ditissima/Hg199)
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
    OutDir=$(dirname $Assembly)/pilon
    Iterations=10
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
    sbatch -p himem $ProgDir/pilon_1_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```
```bash
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    # If split or remove contigs is needed, provide FCSreport file by NCBI.
    touch tmp.txt
    for Assembly in $(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/pilon/pilon_10.fasta ); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/R0905_hybrid_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
```



## Merging long and short read assemblies. Test script


```bash
# Best Hg199 assembly
    for MinIONAssembly in $(ls assembly_VP/flye/N.ditissima/Hg199/racon_10/medaka/medaka/pilon/pilon10_renamed.fasta); do
        Organism=$(echo $MinIONAssembly | rev | cut -f7 -d '/' | rev)
        Strain=$(echo $MinIONAssembly | rev | cut -f6 -d '/' | rev)
        Assembler=$(echo $MinIONAssembly | rev | cut -f8 -d '/' | rev)
        HybridAssembly=$(ls assembly_VP/hybridSPAdes/N.ditissima/Hg199/Hg199_hybrid_renamed.fasta)
        AnchorLength=384699
        OutDir=assembly_VP/merged_assemblies_vtest/"$Assembler"_spades/$Organism/"$Strain"_minion_380k
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/quickmerge_test_script.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    done
```

```bash
# Best R0905 assembly
    for MinIONAssembly in $(ls assembly_VP/canu/N.ditissima/R0905/medaka/medaka/pilon/pilon10_renamed.fasta); do
        Organism=$(echo $MinIONAssembly | rev | cut -f6 -d '/' | rev)
        Strain=$(echo $MinIONAssembly | rev | cut -f5 -d '/' | rev)
        Assembler=$(echo $MinIONAssembly | rev | cut -f7 -d '/' | rev)
        HybridAssembly=$(ls assembly_VP/hybridSPAdes/N.ditissima/R0905/R0905_hybrid_renamed.fasta)
        AnchorLength=524431
        OutDir=assembly_VP/merged_assemblies_vtest/"$Assembler"_spades/$Organism/"$Strain"_minion_520k
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
        sbatch $ProgDir/quickmerge_test_script.sh $MinIONAssembly $HybridAssembly $OutDir $AnchorLength
    done
```

```bash
    for Assembly in $(ls assembly_VP/merged_assemblies_vtest/*/N.ditissima/*/merged_pacbio_hybrid_merged.fasta); do
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
    OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
    sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done
```
