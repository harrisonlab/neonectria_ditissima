# Path test analyses

### Apple cv vs Nd populations

Path tests with 4 apple cultivars and 5 isolates of Nd from different poupulations

```r
# Load libraries
library(data.table)
library(tidyverse) 
library(nlme) 
library(lme4) 
library(emmeans) 
#library(metafuncs)
library(agricolae)

# Load data
dat=read.csv('data.csv')
days<-c(11,25,42,58,66)

# Calculate audpc
dat[,10]=audpc(dat[,5:9],days)
dat[,11]=audpc(dat[,5:9],days,"relative")
colnames(dat)[10]=c('audpc')
colnames(dat)[11]=c('audpc_r')
summary(dat)
write.csv(dt,'dt.csv')

# Get character (or not numeric) columns
cols <- names(dat)[!grepl("numeric",sapply(dat,class))] 
# coerce cols to factors
dat[,(cols):=lapply(.SD,as.factor),.SDcols=cols]

days<-c(11,25,42,58,66)

dat[,10]=audpc(dat[,5:9],days)
dat[,11]=audpc(dat[,5:9],days,"relative")
colnames(dat)[10]=c('audpc')
colnames(dat)[11]=c('audpc_r')


dat[,1]=as.factor(dat[,1])
dat[,2]=as.factor(dat[,2])
dat[,3]=as.factor(dat[,3])
dat[,4]=as.factor(dat[,4])
dat[,10]=as.factor(dat[,10])
dat[,11]=as.factor(dat[,11])

dat2 <- fread("data2.csv")

# get character (or not numeric) columns
cols <- names(dat2)[!grepl("numeric",sapply(dat2,class))] 
# coerce cols to factors
dat2[,(cols):=lapply(.SD,as.factor),.SDcols=cols]


dat2[,group:=as.factor(paste(Plant,Isolate,Block,sep="_"))]



# add experimental unit factor
dat3 <- copy(dat2)

cols <- names(dat3)[grep("numeric",sapply(dat,class))]
dat3 <- dat3[,lapply(.SD,mean),by=c("group","Block","Isolate","Plant"),.SDcols=cols]

ggplot(dat3, aes(x=Plant,y=audpc)) + geom_boxplot() + theme_bw(base_size=12)
ggplot(dat3, aes(x=Isolate,y=audpc)) + geom_boxplot() + theme_bw(base_size=12)



 # using mean of "pseudo"
#m1.1 <- lm(audpc~Block + Cultivar,dat)
m1.1 <- lme(audpc~Block + pseudo + Cultivar + pseudo:Cultivar ,random=~1|group,data=dat2)
anova(m1.1)

m1.1 <- lme (audpc~Block+pseudo+Isolate+Plant + pseudo:Plant + pseudo:Isolate + Isolate:Plant, random=~1|group,data=dat)
m1.2 <- lme (audpc~Block+pseudo+Isolate+Plant + pseudo:Plant random=~1|group,data=dat2)
m1.3 <- lme (audpc~Block+pseudo+Isolate+Plant + pseudo:Plant, random=~1|group,data=dat2)

require(ggplot2)
ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Label))




lmer<_audpc - rootstock + isolate + pseudoreplicate + (1|bloque)


ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Label))


6:02
ranef (model)




library(agricolae)
dt=read.csv('data_c.csv')
days<-c(11,25,42,58,66)

dt[,10]=audpc(dt[,5:9],days)
dt[,11]=audpc(dt[,5:9],days,"relative")
colnames(dt)[10]=c('audpc')
colnames(dt)[11]=c('audpc_r')

summary(dat)
write.csv(dat,'data2.csv')

ggplot(data = dt, aes(x=Isolate, y=audpc)) + geom_boxplot(aes(fill=Plant))  + theme_bw(base_size=12)


# get character (or not numeric) columns
cols <- names(dt)[!grepl("numeric",sapply(dt,class))] 
# coerce cols to factors
dt[,(cols):=lapply(.SD,as.factor),.SDcols=cols]

# get character (or not numeric) columns
cols <- names(dt)[!grepl("numeric",sapply(dt,class))] 
# coerce cols to factors
dt[,(cols):=lapply(.SD,as.factor),.SDcols=cols]


dt[,group:=as.factor(paste(Plant,Isolate,Block,sep="_"))]


m1.1 <- lme (audpc~Block+pseudo+Isolate+Plant + pseudo:Plant + pseudo:Isolate + Isolate:Plant, random=~1|group,data=dt)


ggplot(data = dat, aes(x=Isolate, y=audpc)) + geom_boxplot(aes(fill=Plant))  + theme_bw(base_size=12)




bat=read.csv('data2.csv')
days<-c(11,25,42,58,66)

bat[,10]=audpc(bat[,5:9],days)
bat[,11]=audpc(bat[,5:9],days,"relative")
colnames(bat)[10]=c('audpc')
colnames(bat)[11]=c('audpc_r')

#summary(dat)
write.csv(bat,'bat.csv')

bat <- fread("bat.csv")
cols <- names(bat)[!grepl("numeric",sapply(bat,class))]
bat[,(cols):=lapply(.SD,as.factor),.SDcols=cols]

bat[,group:=as.factor(paste(Plant,Isolate,Block,sep="_"))]
bat2 <- copy(bat)
# get numeric colums
cols <- names(bat)[grep("numeric",sapply(bat,class))]
# get mean of reps for each exp unit
bat <- bat2[,lapply(.SD,mean),by=c("group","Block","Isolate","Plant"),.SDcols=cols]

b<-ggplot(data = bat, aes(x=Isolate, y=audpc)) + geom_boxplot(aes(fill=Plant))  + theme_bw(base_size=12)

b1.1 <- lme (audpc~Block+pseudo+Isolate+Plant + pseudo:Plant + pseudo:Isolate + Isolate:Plant, random=~1|group,data=bat2)

> anova(b1.1)
               numDF denDF  F-value p-value
(Intercept)        1    92 623.4940  <.0001
Block              4    76   0.9792  0.4240
pseudo             1    92  20.7964  <.0001
Isolate            4    76  12.1143  <.0001
Plant              3    76   6.2787  0.0007
pseudo:Plant       3    92   0.3029  0.8232
pseudo:Isolate     4    92   1.4756  0.2160
Isolate:Plant     12    76   0.5343  0.8857


emmeans(b1.1,pairwise~Isolate)

> emmeans(b1.1,pairwise~Isolate)
NOTE: Results may be misleading due to involvement in interactions
$emmeans
 Isolate emmean   SE df lower.CL upper.CL
 Hg199      909 78.2 76      753     1064
 P112      1282 78.2 76     1126     1437
 SVK1       866 78.2 76      711     1022
 U11823     527 78.2 76      371      683
 U11824     782 78.2 76      626      937

Results are averaged over the levels of: Block, pseudo, Plant 
Degrees-of-freedom method: containment 
Confidence level used: 0.95 

$contrasts
 contrast        estimate  SE df t.ratio p.value
 Hg199 - P112      -373.0 111 76 -3.373  0.0100 
 Hg199 - SVK1        42.1 111 76  0.381  0.9954 
 Hg199 - U11823     381.4 111 76  3.450  0.0080 
 Hg199 - U11824     127.0 111 76  1.149  0.7802 
 P112 - SVK1        415.1 111 76  3.755  0.0030 
 P112 - U11823      754.4 111 76  6.823  <.0001 
 P112 - U11824      500.0 111 76  4.522  0.0002 
 SVK1 - U11823      339.3 111 76  3.068  0.0242 
 SVK1 - U11824       84.8 111 76  0.767  0.9392 
 U11823 - U11824   -254.4 111 76 -2.301  0.1559 

Results are averaged over the levels of: Block, pseudo, Plant 
Degrees-of-freedom method: containment 
P value adjustment: tukey method for comparing a family of 5 estimates 

> emmeans(b1.1,pairwise~Plant)
NOTE: Results may be misleading due to involvement in interactions
$emmeans
 Plant emmean   SE df lower.CL upper.CL
 Cox     1008 69.9 76      869     1147
 Gala     819 69.9 76      680      958
 GD       649 69.9 76      509      788
 M9      1017 69.9 76      877     1156

Results are averaged over the levels of: Block, pseudo, Isolate 
Degrees-of-freedom method: containment 
Confidence level used: 0.95 

$contrasts
 contrast   estimate   SE df t.ratio p.value
 Cox - Gala   188.83 98.9 76  1.909  0.2328 
 Cox - GD     359.24 98.9 76  3.633  0.0028 
 Cox - M9      -8.89 98.9 76 -0.090  0.9997 
 Gala - GD    170.41 98.9 76  1.723  0.3190 
 Gala - M9   -197.72 98.9 76 -1.999  0.1973 
 GD - M9     -368.13 98.9 76 -3.722  0.0021 

Results are averaged over the levels of: Block, pseudo, Isolate 
Degrees-of-freedom method: containment 
P value adjustment: tukey method for comparing a family of 4 estimates 