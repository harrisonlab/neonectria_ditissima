# Submission Commands

### Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .

Data is stored on NIAB HPC

```bash
# Metadata 
/projects/neonectria_ditissima/NCBI/SRA_submission
```

```bash
# Submission commands
mkdir -p /projects/neonectria_ditissima/NCBI/SRA_reads/Ndit_PRJNA369735

for isolate in 118923 226-31 Ag02 Ag06 Ag09_A Ag11-B AgN04 ND9 OPC304 R37-15 R41-15 R45-15 R6-17-3 R68-17-C3 SVK2 118924 227-31 Ag05 Ag08 Ag11_A Ag11_C BGV344 ND8 P112 R39-15 R42-15 R6-17-2 R68-17-C2 SVK1; do
for Read in $(ls raw_dna/paired/*/$isolate/*/*.fastq.gz); do
echo $Read;
cp $Read /projects/neonectria_ditissima/NCBI/SRA_reads/Ndit_PRJNA369735
done
done

cd $SubFolder

# Enter user name and password. Copy and paste allowed. It has to be quick.
ncftp -u subftp -p w4pYB9VQ ftp-private.ncbi.nlm.nih.gov
cd uploads/antonio.gomez_emr.ac.uk_IhsSqEFy
mput *
```

### Submission genomes

