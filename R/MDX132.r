install.packages("agricolae")
library(agricolae)

# V1 foler. All timepoints.

AG01=read.csv("C20.csv")
days<-c(12,19,26,33,40,54,76)
AG01[,12]=audpc(AG01[,5:11],days)
AG01[,13]=audpc(AG01[,5:11],days,"relative")
colnames(AG01)[12]=c('audpc')
colnames(AG01)[13]=c('audpc_r')
write.csv(AG01,'C20.csv')

AG02=read.csv("C18.csv")
days<-c(12,20,34,48,55,66)

AG02[,1]=as.factor(AG02[,1])
AG02[,2]=as.factor(AG02[,2])
AG02[,3]=as.factor(AG02[,3])
AG02[,4]=as.factor(AG02[,4])
AG02[,11]=audpc(AG02[,5:10],days)
AG02[,12]=audpc(AG02[,5:10],days,"relative")
colnames(AG02)[11]=c('audpc')
colnames(AG02)[12]=c('audpc_r')
write.csv(AG02,'C18.csv')

AG03=read.csv("Cut1.csv")
days<-c(15,22,30,37,44,52,59)
AG03[,17]=audpc(AG03[,10:16],days)
AG03[,18]=audpc(AG03[,10:16],days,"relative")
colnames(AG03)[17]=c('audpc')
colnames(AG03)[18]=c('audpc_r')
write.csv(AG03,'Cut1.csv')

AG04=read.csv("Cut2.csv")
days<-c(15,22,30,37,44,52,59)
AG04[,17]=audpc(AG04[,10:16],days)
AG04[,18]=audpc(AG04[,10:16],days,"relative")
colnames(AG04)[17]=c('audpc')
colnames(AG04)[18]=c('audpc_r')
write.csv(AG04,'Cut2.csv')

library(lattice)
histogram(~audpc | 'Plant_ID', data=AG01)
histogram(~audpc | 'Plant_ID', data=AG02)
histogram(~audpc | 'Plant_ID', data=AG03)
histogram(~audpc | 'Plant_ID', data=AG04)

trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG01, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG02, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG03, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG04, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))


aggregate(audpc~Plant_ID,AG01,mean)

aggregate(audpc~Plant_ID,AG01,FUN=(function(x){ifelse(sum(x==0)>0 & sum(x !=0) >0, mean(x[x>0]), mean(x))}))

AG01[,14] <- rowMeans(AG01[c('audpc')], na.rm=TRUE)






library(PerformanceAnalytics)
datos<-read.csv('FinalData4.csv')
my_datar <- datos[, c(2:5)] # select just pheno collumn
chart.Correlation(my_data, histogram=TRUE, col=“darkgrey”)

library(agricolae)







 





AGNd09=read.csv('analysis Exp4.csv')
  #AGNd09[,13]=rep(1:2,176)
  #colnames(AGNd08)[15]=c('pseudo')
library(agricolae)

AGNd08[,1]=as.factor(AGNd08[,1])
AGNd08[,2]=as.factor(AGNd08[,2])
AGNd08[,3]=as.factor(AGNd08[,3])
AGNd08[,4]=as.factor(AGNd08[,4])
AGNd08[,5]=as.factor(AGNd08[,5])

days<-c(15,23,30,37,44,51,58)

AGNd08[,13]=audpc(AGNd08[,6:12],days)
AGNd08[,14]=audpc(AGNd08[,6:12],days,"relative")
colnames(AGNd08)[13]=c('audpc')
colnames(AGNd08)[14]=c('audpc_r')

AGNd08[,15]=mean(audpc,Plant_ID)

summary(AGNd08)
write.csv(AGNd08,'analysis Exp3.csv')

library(lattice)
histogram(~audpc | 'Plant ID', data=AGNd08)

library(RColorBrewer)
myColours <- brewer.pal(6,"Blues")

my.settings <- list(
  superpose.polygon=list(col=myColours[2:5], border="transparent"),
  strip.background=list(col=myColours[6]),
  strip.border=list(col="black")
)

library(lattice)
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant ID',data=AGNd08, xlab=("AUDPC") , ylab=("% Disease"))


library (lme4)
scion <- lm(audpc ~ Replicate+Order.in.tray.*Genotype.no.*pseudo + (1| Order.in.tray/Genotype.no.),data=AGNd08)
summary(scion)

library(car)
Anova(scion)

library(effects)
ef <- effect("Genotype.no.", pseudo)
sc <- as.data.frame(ef)
sc <-sc[1:169,]
library(ggplot2)
ggplot(sc, aes(reorder(Genotype.no.,fit),fit)) + geom_point() + coord_fixed(ratio=0.001) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Genotype.no.")

pdf("scion.pdf")

rep1_means=aggregate(AGNd08$audpc~AGNd08$Plant.ID, FUN=mean)

rep1=read.csv('C20audcp.csv')
rep2=read.csv('C18audcp.csv')
rep1_means=aggregate(rep1$audpc~rep1$Plant.ID, FUN=mean)
rep2_means=aggregate(rep2$audpc~rep2$Plant.ID, FUN=mean)
colnames(rep1_means)=c("Plant.ID","audpc")
colnames(rep2_means)=c("Plant.ID","audpc")
merged_data=merge(rep1_means,rep2_means,by="Plant.ID")


plot(merged_data[,2],merged_data[,3])
cor.test(merged_data[,2],merged_data[,3])

Pearson's product-moment correlation

data:  merged_data[, 2] and merged_data[, 3]
t = 3.3063, df = 148, p-value = 0.001186
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.1064804 0.4054849
sample estimates:
      cor
0.2622664
