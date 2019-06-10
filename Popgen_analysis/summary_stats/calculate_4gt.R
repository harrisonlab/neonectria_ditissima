library("PopGenome")
library(ggplot2)

Nd <- c("Ag02", "Ag04", "Ag05", "Ag06", "Ag08", "Ag09_A", "Ag11_A", "Ag11_B", "Ag11_C", "BGV344", "ND8", "ND9", "Hg199", "OPC304", "P112", "R0905", "R37-15", "R39-15", "R41-15", "R42-15", "R45-15", "R6-17-2", "R6-17-3","R68-17-C2", "R68-17-C3", "SVK1", "SVK2")

populations <- list(Nd)
population_names <- c(Nd)
population_no <- length(populations)
pairs <- choose(population_no, 2)

interval <-  10000
jump_size <-  interval / 10

gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)

#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig-containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""])
{
  contig_folder <- paste("contigs/", dir, sep="")
  GENOME.class <- readData(contig_folder, gffpath=gff, include.unknown = TRUE)
  GENOME.class <- set.populations(GENOME.class, populations)

  GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
  GENOME.class.slide <- sliding.window.transform(GENOME.class,width=interval,jump=jump_size,type=2, whole.data=TRUE)
  #per gene
  GENOME.class.split <- recomb.stats(GENOME.class.split)
  fourgamete_split <- get.recomb(GENOME.class.split)
  #per interval
  GENOME.class.slide <- recomb.stats(GENOME.class.slide)
  fourgamete_slide <- get.recomb(GENOME.class.slide)
  ids <- length(GENOME.class.slide@region.names)
  xaxis <- seq(from = 1, to = ids, by = 1)

  #Loop over each population: print figure and table with raw data to file
  for (i in seq_along(population_names))
  {
    fgt <- unlist(fourgamete_split[i])
    file_table = paste(dir, "_", population_names[i], "_4GT_per_gene.txt", sep="")
    file_table2 = paste("genome_", population_names[i], "_4GT_per_gene.txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    fgt_table <- as.data.frame(fourgamete_split[i])
    write.table(fgt_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(fgt_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }

  for (i in seq_along(population_names))
  {
    fgt <- unlist(fourgamete_slide[i])
    #write table with raw data
    slide_table <- paste(dir, "_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
    slide_table2 <- paste("genome_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
    write.table(as.data.frame(fourgamete_slide[i]), file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(as.data.frame(fourgamete_slide[i]), file=slide_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }

}

##Print files with combined results across the entire genome
for (i in seq_along(population_names))
{
  #Four gamete test
  file_table2 = paste("genome_", population_names[i], "_4GT_per_gene.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_4GT_per_gene.pdf", sep="")
  fgt_plot <- ggplot(x, aes(x=x[,2])) + geom_histogram(colour="black", fill="cornsilk") + xlab("Four gamete test") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,2], n = 10))
  ggsave(file_hist, fgt_plot)

  file_table2 = paste("genome_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_4GT_per_sliding_window.pdf", sep="")
  fgt_plot <- ggplot(x, aes(x=x[,2])) + geom_histogram(colour="black", fill="cornsilk") + xlab("Four gamete test") + ylab("Number of intervals") + scale_x_continuous(breaks = pretty(x[,2], n = 10))
  ggsave(file_hist, fgt_plot)
}
