#README for popgen analyses


Tools for population genetics analyses. The order in which main bash scripts in the folder are executed:

1. pre_SNP_calling_cleanup.md - Inital set up commands before running SNP calling
2. SNP_calling_multithreaded.md - Runs the core SNP calling pipeline
3. sub_SNP_calling_multithreaded.md - Script that is qsubbed for core SNP calling
4. SNP_analysis.md - Inital sets of analysis on resulting .vcf file
