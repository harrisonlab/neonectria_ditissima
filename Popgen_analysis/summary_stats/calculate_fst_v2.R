library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pop1 <- c("Ag02", "Ag04", "Ag05", "Ag06", "Ag08", "Ag09_A", "Ag11_A", "Ag11_B", "Ag11_C", "BGV344", "ND8", "ND9", "Hg199", "OPC304", "P112", "R0905", "R37-15", "R39-15", "R41-15", "R42-15", "R45-15", "R6-17-2", "R6-17-3")
Pop2 <- c("R68-17-C2", "R68-17-C3")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pop1,Pop2)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pop1","Pop2")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
interval <-  1000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""]){
  contig_folder <- paste("contigs/", dir, sep = "")
  GENOME.class <- readData(contig_folder, gffpath = gff, include.unknown = TRUE)
  GENOME.class <- set.populations(GENOME.class, populations)

###############################################################################

#### Gene based analysis
GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
GENOME.class.split <- F_ST.stats(GENOME.class.split)
get.F_ST(GENOME.class.split)

FST_all <- GENOME.class.split@nuc.F_ST.vs.all
FST_all_d <- as.data.frame(FST_all)
FST_pairwise <- GENOME.class.split@nuc.F_ST.pairwise
Hudson_KST <- GENOME.class.split@Hudson.K_ST

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_total_FST_per_gene.pdf",
  sep = "")
  fst_plot <- ggplot(FST_all_d, aes(x = FST_all_d[, i])) +
  geom_histogram(colour = "black", fill = "darkseagreen") +
  xlab(expression(paste("Total F"[ST], " per gene"))) +
  ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FST_all_d[, i], n = 10)) +
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  axis.text = element_text(size = 14), axis.title = element_text(size = 18))
  ggsave(file_hist, fst_plot)
  file_table <- paste(dir, "_", population_names[i], "_total_FST_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_total_FST_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  fst_table <- cbind(gene_ids, FST_all[, i])
  write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

for (i in seq(pairs)){
  FST_pairwise_d <- as.data.frame(as.vector(FST_pairwise[i, ]))
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
  file_hist <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".pdf",
  sep = "")
  fst_plot <- ggplot(FST_pairwise_d, aes(x = FST_pairwise_d[, 1])) +
  geom_histogram(colour = "black", fill = "cadetblue") +
  xlab(expression(paste("Pairwise F"[ST], " per gene"))) +
  ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FST_pairwise_d[, 1], n = 10)) +
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  axis.text = element_text(size = 14), axis.title = element_text(size = 18))
  ggsave(file_hist, fst_plot)
  file_table <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  fst_table <- cbind(gene_ids, FST_pairwise_d[, 1])
  write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE)
  write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE, append = TRUE)
}

### Plot Hudson KST
Hudson_KST_d <- as.data.frame(as.vector(Hudson_KST[1, ]))
file_hist <- paste(dir, "_Hudson_KST_per_gene", ".pdf", sep = "")
fst_plot <- ggplot(Hudson_KST_d, aes(x = Hudson_KST_d[, 1])) +
geom_histogram(colour = "black", fill = "springgreen") +
xlab(expression(paste("Hudson K"[ST], " per gene"))) +
ylab("Number of genes") +
scale_x_continuous(breaks = pretty(Hudson_KST_d[, 1], n = 10)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
ggsave(file_hist, fst_plot)
file_table <- paste(dir, "_Hudson_KST_per_gene", ".txt", sep = "")
file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
current_gff <- paste(gff, "/", dir, ".gff", sep = "")
gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
    feature = FALSE, extract.gene.names = TRUE)
fst_table <- cbind(gene_ids, Hudson_KST_d[, 1])
write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
col.names = FALSE, row.names = FALSE)
write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
col.names = FALSE, row.names = FALSE, append = TRUE)

#### Sliding window analysis (interval)
GENOME.class.slide <- sliding.window.transform(GENOME.class, width = interval,
    jump = jump_size, type = 2, whole.data = TRUE)
GENOME.class.slide <- F_ST.stats(GENOME.class.slide, mode = "nucleotide")
FST_all_slide <- GENOME.class.slide@nuc.F_ST.vs.all
FST_pairwise_slide <- GENOME.class.slide@nuc.F_ST.pairwise
FST_all_slide_d <- as.data.frame(FST_all_slide)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names)){
  file_slide <- paste(dir, "_", population_names[i],
  "_total_FST_per_sliding_window.pdf", sep = "")
  slide_plot <- ggplot(FST_all_slide_d, aes(x = xaxis,
      y = FST_all_slide_d[, i])) +
      geom_smooth(colour = "black", fill = "plum") +
      xlab("Contig coordinate (kbp)") +
      ylab(expression(paste("Total F"[ST], " per interval"))) +
      scale_x_continuous(breaks = pretty(xaxis, n = 10)) +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 18))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i],
  "_total_FST_per_sliding_window.txt", sep = "")
  write.table(FST_all_slide[, i], file = slide_table, sep = "\t",
  quote = FALSE, col.names = FALSE)
}

#Plot pairwise FST
for (i in seq(pairs)){
  FST_pairwise_slide_d <- as.data.frame(as.vector(FST_pairwise_slide[i, ]))
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise_slide)[i])
  file_hist <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
  ".pdf", sep = "")
  slide_plot <- ggplot(FST_pairwise_slide_d, aes(x = xaxis,
      y = FST_pairwise_slide_d[, 1])) +
      geom_smooth(colour = "black", fill = "slateblue") +
      xlab("Contig coordinate (kbp)") +
      ylab(expression(paste("Pairwise F"[ST], " per interval"))) +
      scale_x_continuous(breaks = pretty(xaxis, n = 10)) +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 18))
  ggsave(file_hist, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
  ".txt", sep = "")
  fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[, 1])
  write.table(fst_table, file = slide_table, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE)
}

}

##Genomewide results
for (i in seq_along(population_names)){
#Total FST
file_table2 <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 3])) +
geom_histogram(colour = "black", fill = "darkseagreen") +
xlab(expression(paste("Total F"[ST], " per gene"))) +
ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 3], n = 10)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
ggsave(file_hist, fst_plot)
}

#Hudson KST
file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_Hudson_KST_per_gene_all", ".pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 2])) +
geom_histogram(colour = "black", fill = "springgreen") +
xlab(expression(paste("Hudson K"[ST], " per gene"))) +
ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 2], n = 10)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
ggsave(file_hist, fst_plot)


for (i in seq(pairs)){
  #Pairwise FST
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
  file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_pairwise_FST_per_gene_", labelling, ".pdf",
  sep = "")
  fst_plot <- ggplot(x, aes(x = x[, 2])) +
  geom_histogram(colour = "black", fill = "cadetblue") +
  xlab(expression(paste("Pairwise F"[ST], " per gene"))) +
  ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(x[, 2], n = 10)) +
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  axis.text = element_text(size = 14), axis.title = element_text(size = 18))
  ggsave(file_hist, fst_plot)
}
