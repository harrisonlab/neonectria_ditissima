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
histogram(~audpc | 'Plant_ID', data=AG03)

trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG01, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG02, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG03, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG04, breaks=100, xlab=("AUDPC") , ylab=("% Disease"))

#Replace 0 values with NA.
AG01$audpc[AG01$audpc==0]<-NA
AG02$audpc[AG02$audpc==0]<-NA
AG03$audpc[AG03$audpc==0]<-NA
AG04$audpc[AG04$audpc==0]<-NA

#Mean of pseudoreps audpc values
AG011=aggregate(audpc~Plant_ID,AG01,FUN=mean, na.rm=TRUE)
AG021=aggregate(audpc~Plant_ID,AG02,FUN=mean, na.rm=TRUE)
AG031=aggregate(audpc~Plant_ID,AG03,FUN=mean, na.rm=TRUE)
AG041=aggregate(audpc~Plant_ID,AG01,FUN=mean, na.rm=TRUE)

histogram(~audpc | 'Plant_ID', data=AG011)
histogram(~audpc | 'Plant_ID', data=AG021)
histogram(~audpc | 'Plant_ID', data=AG031)
histogram(~audpc | 'Plant_ID', data=AG041)

trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG01, breaks=20, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG021, breaks=20, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG031, breaks=20, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG041, breaks=20, xlab=("AUDPC") , ylab=("% Disease"))


# V2 foler. Remove last timepoint. Possible more accurate values.

library(agricolae)

AG01=read.csv("C20.csv")
days<-c(12,19,26,33,40,54)
AG01[,11]=audpc(AG01[,5:10],days)
AG01[,12]=audpc(AG01[,5:10],days,"relative")
colnames(AG01)[11]=c('audpc')
colnames(AG01)[12]=c('audpc_r')
write.csv(AG01,'C20.csv')

AG02=read.csv("C18.csv")
days<-c(12,20,34,48,55)
AG02[,10]=audpc(AG02[,5:9],days)
AG02[,11]=audpc(AG02[,5:9],days,"relative")
colnames(AG02)[10]=c('audpc')
colnames(AG02)[11]=c('audpc_r')
write.csv(AG02,'C18.csv')

AG03=read.csv("Cut1.csv")
days<-c(15,22,30,37,44,52)
AG03[,16]=audpc(AG03[,10:15],days)
AG03[,17]=audpc(AG03[,10:15],days,"relative")
colnames(AG03)[16]=c('audpc')
colnames(AG03)[17]=c('audpc_r')
write.csv(AG03,'Cut1.csv')

AG04=read.csv("Cut2.csv")
days<-c(15,22,30,37,44,52)
AG04[,16]=audpc(AG04[,10:15],days)
AG04[,17]=audpc(AG04[,10:15],days,"relative")
colnames(AG04)[16]=c('audpc')
colnames(AG04)[17]=c('audpc_r')
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

#Replace 0 values with NA.
AG01$audpc[AG01$audpc==0]<-NA
AG02$audpc[AG02$audpc==0]<-NA
AG03$audpc[AG03$audpc==0]<-NA
AG04$audpc[AG04$audpc==0]<-NA

#Mean of pseudoreps audpc values
AG011=aggregate(audpc~Plant_ID,AG01,FUN=mean, na.rm=TRUE)
AG021=aggregate(audpc~Plant_ID,AG02,FUN=mean, na.rm=TRUE)
AG031=aggregate(audpc~Plant_ID,AG03,FUN=mean, na.rm=TRUE)
AG041=aggregate(audpc~Plant_ID,AG04,FUN=mean, na.rm=TRUE)

histogram(~audpc | 'Plant_ID', data=AG011)
histogram(~audpc | 'Plant_ID', data=AG021)
histogram(~audpc | 'Plant_ID', data=AG031)
histogram(~audpc | 'Plant_ID', data=AG041)

trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG011, breaks=30, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG021, breaks=30, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG031, breaks=30, xlab=("AUDPC") , ylab=("% Disease"))
trellis.device(width=7, height=5, new=FALSE, color=FALSE)
histogram(~audpc|'Plant_ID',data=AG041, breaks=30, xlab=("AUDPC") , ylab=("% Disease"))

#Merged data of each experiment
merged_data=merge(AG011,AG021,by="Plant_ID")
merged_cut=merge(AG031,AG041,by="Plant_ID")

plot(merged_data[,2],merged_data[,3])
cor.test(merged_data[,2],merged_data[,3])
plot(merged_cut[,2],merged_cut[,3])
cor.test(merged_cut[,2],merged_cut[,3])

#Merged replicate data in a csv file. Not needed
write.csv(merged_data,'Merged1.csv')
write.csv(merged_cut,'Merged2.csv')

install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)

my_data<-merged_data[,c(2:3)] # select just pheno collumn
chart.Correlation(my_data, histogram=TRUE)
my_data2<-merged_cut[,c(2:3)] # select just pheno collumn
chart.Correlation(my_data2, histogram=TRUE)
