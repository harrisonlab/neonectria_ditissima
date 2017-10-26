```R
sampling<-read.table("E:\\Isolates.txt", sep="",header=T,na.string="*")

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
p<-p+scale_y_continuous(breaks=seq(0,450,50))
p<-p+coord_cartesian(ylim=c(0,450))
p
```
