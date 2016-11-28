AGNd05=read.csv('ash_data.csv')

library(agricolae)

AGNd05[,1]=as.factor(AGNd05[,1])
AGNd05[,2]=as.factor(AGNd05[,2])
AGNd05[,3]=as.factor(AGNd05[,3])
AGNd05[,4]=as.factor(AGNd05[,4])

days<-c(22,36,57,70,80,90)

AGNd05[,11]=audpc(AGNd05[,5:10],days)
AGNd05[,12]=audpc(AGNd05[,5:10],days,"relative")
colnames(AGNd05)[11]=c('audpc')
colnames(AGNd05)[12]=c('audpc_r')

```bash
#library(lattice)

#histogram(~audpc | Accession, data=AGNd05)
```

colnames(AGNd05)

scion <- lm(audpc ~ Accession  ,data=AGNd05)
summary(scion)

```bash
Call:
lm(formula = audpc ~ Accession, data = AGNd05)

Residuals:
    Min      1Q  Median      3Q     Max
-3923.0  -782.2     0.0   725.7  5870.0

Coefficients:
                       Estimate Std. Error t value Pr(>|t|)
(Intercept)              1126.3     1026.7   1.097 0.276795
Accession108014          -903.8     1452.0  -0.622 0.535866
Accession109001          2641.0     1452.0   1.819 0.073677 .
Accession203014          5070.3     1452.0   3.492 0.000882 ***
Accession203015          1591.8     1452.0   1.096 0.277106
Accession204003          -901.3     2053.4  -0.439 0.662199
Accession204009          2174.3     1452.0   1.498 0.139253
Accession204010          2843.8     1452.0   1.959 0.054588 .
Accession301006          -335.8     1623.3  -0.207 0.836773
Accession302007          1432.7     1623.3   0.883 0.380839
Accession302008          1824.5     1452.0   1.257 0.213547
Accession302020           465.8     1452.0   0.321 0.749401
Accession302025          1034.0     1452.0   0.712 0.479007
Accession303002          4085.5     1452.0   2.814 0.006524 **
Accession303003         -1126.3     1452.0  -0.776 0.440810
Accession304007          1049.9     1623.3   0.647 0.520135
Accession304008           583.7     1452.0   0.402 0.689056
Accession305002           856.7     1623.3   0.528 0.599549
Accession305008          -176.2     1452.0  -0.121 0.903816
Accession305023          6954.2     2053.4   3.387 0.001223 **
Accession401006          -375.5     1452.0  -0.259 0.796775
Accession402004          3296.2     2053.4   1.605 0.113443
Accession402014          4766.8     1452.0   3.283 0.001678 **
Accession403001          1353.9     1623.3   0.834 0.407415
Accession403049          1887.3     1452.0   1.300 0.198388
Accession404004           646.3     1452.0   0.445 0.657741
Accession404018           460.8     1452.0   0.317 0.752000
Accession404027          -488.6     1623.3  -0.301 0.764426
Accession405001           576.4     1623.3   0.355 0.723716
Accession405017           827.8     1452.0   0.570 0.570605
Accession406001           460.3     1452.0   0.317 0.752260
Accession406004          4563.8     1452.0   3.143 0.002548 **
AccessionF.americana    -1126.3     1452.0  -0.776 0.440810
AccessionF.mandsch.     -1126.3     1452.0  -0.776 0.440810
AccessionF.mariesii      1187.0     1452.0   0.818 0.416716
AccessionF.ornus         -339.8     1452.0  -0.234 0.815704
AccessionGreat g. w.      391.8     1452.0   0.270 0.788146
AccessionPoll Tail Co.   3008.7     1452.0   2.072 0.042351 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1778 on 63 degrees of freedom
Multiple R-squared:  0.6174,	Adjusted R-squared:  0.3927
F-statistic: 2.748 on 37 and 63 DF,  p-value: 0.0002004
```

library(car)
Anova(scion)

```bash
Anova Table (Type II tests)

Response: audpc
             Sum Sq Df F value    Pr(>F)
Accession 321497802 37  2.7477 0.0002004 ***
Residuals 199224309 63
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

library(effects)
ef <- effect("Accession", scion)
sc <- as.data.frame(ef)
sc <-sc[1:37,]
library(ggplot2)
ggplot(sc, aes(reorder(Accession,fit),fit)) + geom_point() + coord_fixed(ratio=0.0001) + geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("AUDPC") + xlab("Accession")

pdf("scion.pdf")
