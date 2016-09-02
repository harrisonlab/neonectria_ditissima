
AGNd04=read.csv('canker_data.csv')
AGNd04[,15]=rep(1:2,132)
colnames(AGNd04)[15]=c('pseudo')
library(agricolae)

AGNd04[,1]=as.factor(AGNd04[,1])
AGNd04[,2]=as.factor(AGNd04[,2])
AGNd04[,3]=as.factor(AGNd04[,3])
AGNd04[,4]=as.factor(AGNd04[,4])
AGNd04[,5]=as.factor(AGNd04[,5])
AGNd04[,15]=as.factor(AGNd04[,15])

#days<-c(14,18,21,27,34,39,45,49,54)
days<-c(14,18,21,27,34)
AGNd04[,16]=audpc(AGNd04[,6:10],days)
AGNd04[,17]=audpc(AGNd04[,6:10],days,"relative")
colnames(AGNd04)[16]=c('audpc')
colnames(AGNd04)[17]=c('audpc_r')
library(lattice)

histogram(~audpc | Cultivar*Treatment, data=AGNd04)

#colnames(AGNd04)
# [1] "Plant.ID"  "Cultivar"  "Tray"      "Cabinet"   "Treatment" "X1"        "X2"        "X3"
# [9] "X4"        "X5"        "X6"        "X7"        "X8"        "X9"        "pseudo"    "audpc"
# [17] "audpc_r"


# Fixed effect model =  cabinet+inoculum*cultivar*pseudorep
# Random effect model = cabinet/tray/cultivar
library (lme4)
scion <- lmer(audpc ~ Cabinet+Treatment*Cultivar*pseudo + (1 |Tray/Cultivar) ,data=AGNd04)
summary(scion)

library(car)
Anova(scion)

library(effects)
ef <- effect("Cultivar", scion)
sc <- as.data.frame(ef)
sc <-sc[c(1,2,3,4,5,6,7,8,10,11),]
library(ggplot2)
ggplot(sc, aes(reorder(Cultivar,fit),fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")

pdf("scion_35.pdf")
ggplot(sc, aes(reorder(Cultivar,fit),fit)) + geom_point()+ coord_fixed(ratio=0.02) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar") + expand_limits( y = 0)
dev.off()


library(effects)
ef <- effect("Treatment:Cultivar", scion)
sc <- as.data.frame(ef)
library(ggplot2)
ggplot(sc, aes(reorder(Treatment:Cultivar,fit),fit,fill=Treatment,colour=Treatment)) + geom_point(stat="identity",position="dodge")+ geom_errorbar(aes(ymin=fit-se, ymax=fit+se, postion="dodge",stat="identity",position="dodge"), width=0.4)  + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")

pdf("scionbyisolate.pdf")
ggplot(sc, aes(reorder(Treatment:Cultivar,fit),fit,fill=Treatment,colour=Treatment)) + geom_point(stat="identity",position="dodge")+ coord_fixed(ratio=0.08)+ geom_errorbar(aes(ymin=fit-se, ymax=fit+se, postion="dodge",stat="identity",position="dodge"), width=0.4)  + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()


write.csv(AGNd04,'AGNd04.csv')
