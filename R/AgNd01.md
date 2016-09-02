library(agricolae)


#ROOTSTOCK ANALYSIS

rootstock_data=read.csv("rootstocks.csv")
rootstock_data[,1]=as.factor(rootstock_data[,1])
rootstock_data[,2]=as.factor(rootstock_data[,2])

days<-c(1,9,24,48,70)
evaluation1<-unname(rootstock_data[1:65,c(1,2,3,4,5,6,9,12,15,18)])
evaluation2<-unname(rootstock_data[1:65,c(1,2,3,4,5,7,10,13,16,19)])
evaluation3<-unname(rootstock_data[1:65,c(1,2,3,4,5,8,11,14,17,20)])
colnames(evaluation1)=c("tray","cabinet","plant_id","cultivar","inoculum","t1","t9","t24","t48","t70")
colnames(evaluation2)=c("tray","cabinet","plant_id","cultivar","inoculum","t1","t9","t24","t48","t70")
colnames(evaluation3)=c("tray","cabinet","plant_id","cultivar","inoculum","t1","t9","t24","t48","t70")
combined_rootstocks=rbind(rbind(evaluation1,evaluation2),evaluation3)
combined_rootstocks[,11]=as.factor(c(rep("A",65),rep("B",65),rep("C",65)))

colnames(combined_rootstocks)=c("tray","cabinet","plant_id","cultivar","inoculum","t1","t9","t24","t48","t70","pseudorep")
combined_rootstocks[,12]=audpc(combined_rootstocks[,6:10],days)
colnames(combined_rootstocks)=c("tray","cabinet","plant_id","cultivar","inoculum","t1","t9","t24","t48","t70","pseudorep","abs")
write.csv(combined_rootstocks,"combined_roots.csv")
par(mar=c(8,4,1,1))
plot(abs~cultivar,data=combined_rootstocks,las=2,xlab=NULL)

library (lme4)
rootstock <- lmer(abs ~ cultivar*pseudorep  + (1 |tray/cultivar) ,data=combined_rootstocks)
summary(rootstock)

library(car)
Anova(rootstock)


library(effects)
ef <- effect("cultivar", rootstock)
rs <- as.data.frame(ef)
library(ggplot2)
ggplot(rs, aes(cultivar,fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ ylab("AUDPC") + xlab("Rootstock")





###SCION###
canker_data=read.csv("canker_assessment2.csv")
canker_data[,1]=as.factor(canker_data[,1])
canker_data[,2]=as.factor(canker_data[,2])

#SEPARATE OUT PSEUDOREPS
days<-c(12,16,22,27,31,35)
evaluation1<-unname(canker_data[1:144,c(1,2,3,4,5,6,9,12,15,18,21)])
evaluation2<-unname(canker_data[1:144,c(1,2,3,4,5,7,10,13,16,19,22)])
evaluation3<-unname(canker_data[1:144,c(1,2,3,4,5,8,11,14,17,20,23)])

#SORT OUT PSEUDOREPS
colnames(evaluation1)=c("tray","cabinet","plant_id","cultivar","inoculum","t12","t16","t22","t27","t31","t35")
colnames(evaluation2)=c("tray","cabinet","plant_id","cultivar","inoculum","t12","t16","t22","t27","t31","t35")
colnames(evaluation3)=c("tray","cabinet","plant_id","cultivar","inoculum","t12","t16","t22","t27","t31","t35")
combined_canker=rbind(rbind(evaluation1,evaluation2),evaluation3)
combined_canker[,12]=as.factor(c(rep("A",144),rep("B",144,),rep("C",144)))
colnames(combined_canker)=c("tray","cabinet","plant_id","cultivar","inoculum","t12","t16","t22","t27","t31","t35","pseudorep")

#AUDPC
combined_canker[,13]=audpc(combined_canker[,6:11],days)
colnames(combined_canker)=c("tray","cabinet","plant_id","cultivar","inoculum","t12","t16","t22","t27","t31","t35","pseudorep","abs")

write.csv(combined_canker,"combined_scion.csv")

#PLOT OF RAW DATA
par(mar=c(8,4,1,1))
plot(abs~cultivar,data=combined_canker,las=2,xlab=NULL)


#MODEL
#REML variance components analysis

# Fixed effect model =  cabinet+inoculum*cultivar*pseudorep
# Random effect model = cabinet/tray/cultivar
library (lme4)
#scion <- lmer(abs ~ tray + cultivar + pseudorep  + tray:pseudorep + cultivar:pseudorep + (1 |tray:cultivar) ,data=combined_canker)
#THIS SHOULD BE THE SAME AS GENSTAT
#scion <- lmer(abs ~ cabinet+inoculum*cultivar*pseudorep + (1 |cabinet/tray/cultivar/pseudorep) ,data=combined_canker)
#BUT ONLY THIS ONE RUNS
scion <- lmer(abs ~ cabinet+inoculum*cultivar*pseudorep + (1 |tray/cultivar) ,data=combined_canker)
summary(scion)

library(car)
Anova(scion)


library(effects)
ef <- effect("cultivar", scion)
sc <- as.data.frame(ef)
sc <-sc[2:12,]
library(ggplot2)
ggplot(sc, aes(reorder(cultivar,fit),fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")



#PLOTS
png("scion.png",height=400,width=600)
ggplot(sc, aes(reorder(cultivar,fit),fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()

png("rootstock.png",height=400,width=600)
ggplot(rs, aes(reorder(cultivar,fit),fit)) + geom_point() + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ ylab("AUDPC") + xlab("Rootstock")
dev.off()

pdf("scion.pdf")
ggplot(sc, aes(reorder(cultivar,fit),fit)) + geom_point() + coord_fixed(ratio=0.01) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Cultivar")
dev.off()
pdf("rootstock.pdf")
ggplot(rs, aes(reorder(cultivar,fit),fit)) + geom_point() + coord_fixed(ratio=0.0005) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ ylab("AUDPC") + xlab("Rootstock")
dev.off()
