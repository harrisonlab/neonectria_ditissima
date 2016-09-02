library(agricolae)

leafscar=read.csv("AGND03 Analysis data2.csv")
leafscar[,13]=rep(1:5,40)
days<-c(70,77,90,104,113,122,133,139,153)
colnames(leafscar)[13]=c('pseudo')


leafscar[,1]=as.factor(leafscar[,1])
leafscar[,2]=as.factor(leafscar[,2])
leafscar[,3]=as.factor(leafscar[,3])
leafscar[,13]=as.factor(leafscar[,13])
leafscar[,14]=audpc(leafscar[,4:12],days)
leafscar[,15]=audpc(leafscar[,4:12],days,"relative")
colnames(leafscar)[14]=c('audpc')
colnames(leafscar)[15]=c('audpc_r')
library (lattice)

histogram(~audpc | Cultivar, data=leafscar)


scion <- lm(audpc ~ Block+Cultivar  ,data=leafscar)
summary(scion)
library(car)
Anova(scion)

library(effects)
ef <- effect("Cultivar", scion)
sc <- as.data.frame(ef)
sc <-sc[1:10,]
library(ggplot2)
ggplot(sc, aes(reorder(Cultivar,fit),fit)) + geom_point() + coord_fixed(ratio=0.01) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")




pdf("scion.pdf")
ggplot(sc, aes(reorder(Cultivar,fit),fit)) + geom_point() + coord_fixed(ratio=0.002) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()


#PLOTS
png("scion.png",height=400,width=600)
ggplot(sc, aes(reorder(cultivar,fit),fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()

pdf("scion.pdf")
ggplot(sc, aes(reorder(cultivar,fit),fit)) + geom_point() + coord_fixed(ratio=0.01) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()
