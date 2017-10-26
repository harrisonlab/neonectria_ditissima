
```R
AGC=read.csv('Data.csv')
library(agricolae)

AGC[,1]=as.factor(AGC[,1])
AGC[,2]=as.factor(AGC[,2])
AGC[,3]=as.factor(AGC[,3])

days<-c(10,16,28,36,51,62,71)

AGC[,11]=audpc(AGC[,4:10],days)
AGC[,12]=audpc(AGC[,4:10],days,"relative")
colnames(AGC)[11]=c('audpc')
colnames(AGC)[12]=c('audpc_r')

summary(AGC)
write.csv(AGC,'Data.csv')
```
```R
sampling<-read.table("Isolates.txt", sep="",header=T,na.string="*")

attach(sampling)

library(ggplot2)
library(grid)
library(sciplot) ## for the se function

stderr <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
lowsd <- function(x){return(mean(x)-stderr(x))}
highsd <- function(x){return(mean(x)+stderr(x))}

p<-ggplot(data=sampling,aes(Treatment,AUDPC))
p<-p+stat_summary(data=sampling, fun.y=mean, geom="point")
p<-p+stat_summary(data=sampling, fun.y=mean, fun.ymin=lowsd, fun.ymax=highsd, geom="errorbar", width=0.5)
p<-p+aes(reorder(Treatment,Audpc),Audpc)
p<-p+theme_bw(base_size=12)
p<-p+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p<-p+ylab("AUDPC")+xlab("Isolate")
p<-p+theme(aspect.ratio = 2/(1+sqrt(5)))
p<-p+scale_y_continuous(breaks=seq(0,1800,200))
p<-p+coord_cartesian(ylim=c(0,1800))
p
```
