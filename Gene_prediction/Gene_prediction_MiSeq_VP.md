# Gene Prediction

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
    for Strain in 118923 118924 226-31 227-31; do
        for Assembly in $(ls repeat_masked_VP/SPAdes_assembly/N.ditissima/$Strain/*unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
        done
    done

    for Strain in Ag02 Ag05 Ag08 Ag11_A Ag11_C Hg199 ND9 P112 R37-15 R41-15 R45-15 R6-17-3 R68-17-C2 SVK2 Ag04 Ag06 Ag09_A Ag11_B BGV344 ND8 OPC304 R0905 R39-15 R42-15 R6-17-2 R68-17-C3 SVK1 NMaj; do
        for Assembly in $(ls repeat_masked/N.*/$Strain/*unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
        Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev) 
        echo "$Organism - $Strain"
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
        BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
        OutDir=gene_pred/BUSCO/$Organism/$Strain/busco_sordariomycetes_obd10
        sbatch -p long $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
        done
    done
```
This was previously done with Busco v3.0.2 and the sordarymyceta_odb9. Re run to update the phylogeny

