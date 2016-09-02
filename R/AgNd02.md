seedling=read.csv("seedling4.csv")
seedling[,11]=rep(1:3,154)
colnames(seedling)[11]=c("pseudo")
library (reshape2)


melt_seedling=melt(seedling, id.vars=c("Plant.ID","Cross.No.","pseudo"), measure.vars=c("X1_11","X2_14","X3_17","X4_21","X5_25","X6_28","X7_31"))
melt_seedling[,1]=as.factor(melt_seedling[,1])
melt_seedling[,2]=as.factor(melt_seedling[,2])
melt_seedling[,3]=as.factor(melt_seedling[,3])
melt_seedling[,4]=as.factor(melt_seedling[,4])
write.csv(melt_seedling,"melted_seedling2.csv")


library(agricolae)
days<-c(11,14,17,21,25,28,31)


seedling[,1]=as.factor(seedling[,1])
seedling[,2]=as.factor(seedling[,2])
seedling[,11]=as.factor(seedling[,11])
seedling[,12]=audpc(seedling[,3:9],days)
seedling[,13]=audpc(seedling[,3:9],days,"relative")


colnames(seedling)[12]=c("audpc")
colnames(seedling)[13]=c("audpc_r")

library(RColorBrewer)
myColours <- brewer.pal(6,"Blues")

my.settings <- list(
  superpose.polygon=list(col=myColours[2:5], border="transparent"),
  strip.background=list(col=myColours[6]),
  strip.border=list(col="black")
)

library(lattice)
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|Cross.No.,data=seedling)


pdf("seedling.pdf")
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|Cross.No.,data=seedling,  xlab=("AUDPC") , ylab=("% Disease") )
dev.off()



tapply(seedling$audpc, seedling$Cross.No.,median)
tapply(seedling$audpc, seedling$Cross.No.,IQR)


#MEDIAN
MDX051   MDX052   MDX053   MDX054   MDX057   MDX060   MDX061   MDX063   MDX064   MDX065   MDX068
28.5825  60.4125  25.9200 181.4150 266.6050 234.0600 180.5250 244.5875  85.8000 105.0525 253.3900
#IQR
MDX051   MDX052   MDX053   MDX054   MDX057   MDX060   MDX061   MDX063   MDX064   MDX065   MDX068
101.4975 108.8000 148.4050 248.9513 227.0325 185.8450 251.2275 292.9613 242.6900 158.5325 251.9625
>




aov_fit=aov(seedling$audpc,seedling$Cross.No)
TukeyHSD(aov_fit)

library(pgirmess)
kruskalmc(seedling$audpc~seedling$Cross.No)
