```R
#Plot coverage of each isolate against the R0905 reference genome.
library(readr)
setwd("~/Desktop/Grouped")
df_R0905 <- read_delim("vs_R0905_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("Ag02", "Ag04", "Ag05", "Ag06", "Ag08", "Ag09_A", "Ag11_A", "Ag11_B", "Ag11_C", "BGV344", "Hg199", "ND8", "ND9", "OPC304", "P112", "R0905", "R37-15", "R39-15", "R41-15", "R42-15", "R45-15", "R6-17-2", "R6-17-3", "R68-17-C2", "R68-17-C3", "SVK1", "SVK2", "NMaj"))), trim_ws = TRUE)
myFun <- function(x) {
c(min = min(x), max = max(x),
mean = mean(x), median = median(x),
std = sd(x))
}
colnames(df_R0905) <- c("contig","position", "depth", "strain")
df_R0905$treatment <- paste(df_R0905$strain , df_R0905$contig)
tapply(df_R0905$depth, df_R0905$treatment, myFun)
df2 <- cbind(do.call(rbind, tapply(df_R0905$depth, df_R0905$treatment, myFun)))
write.csv(df2, 'R0905_contig_coverage.csv')
df_R0905$depth <- ifelse(df_R0905$depth > 100, 100, df_R0905$depth)
# install.packages("ggplot2")
library(ggplot2)
require(scales)
for (i in 1:50){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=df_R0905[df_R0905$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 28, ncol = 1, strip.position = "left")
outfile = paste("R0905_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}

#Plot coverage of each isolate against the Hg199 CSAR genome.
library(readr)
setwd("~/Desktop/Grouped_Hg199")
df_Hg199 <- read_delim("vs_Hg199_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("Ag02", "Ag04", "Ag05", "Ag06", "Ag08", "Ag09_A", "Ag11_A", "Ag11_B", "Ag11_C", "BGV344", "Hg199", "ND8", "ND9", "OPC304", "P112", "R0905", "R37-15", "R39-15", "R41-15", "R42-15", "R45-15", "R6-17-2", "R6-17-3", "R68-17-C2", "R68-17-C3", "SVK1", "SVK2", "NMaj"))), trim_ws = TRUE)
myFun <- function(x) {
c(min = min(x), max = max(x),
mean = mean(x), median = median(x),
std = sd(x))
}
colnames(df_Hg199) <- c("contig","position", "depth", "strain")
df_Hg199$treatment <- paste(df_Hg199$strain , df_Hg199$contig)
tapply(df_Hg199$depth, df_Hg199$treatment, myFun)
df2 <- cbind(do.call(rbind, tapply(df_Hg199$depth, df_Hg199$treatment, myFun)))
write.csv(df2, 'Hg199_contig_coverage.csv')
df_Hg199$depth <- ifelse(df_Hg199$depth > 100, 100, df_Hg199$depth)
# install.packages("ggplot2")
library(ggplot2)
require(scales)
for (i in 1:50){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=df_Hg199[df_Hg199$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 28, ncol = 1, strip.position = "left")
outfile = paste("Hg199_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}
```
