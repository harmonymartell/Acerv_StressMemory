##USAGE "Thermal Priming Costs Outweigh Benefits in the Staghorn Coral Acropora cervicornis (Lamarck 1816)"
# written by Harmony Martell July 2022

# Stress Memory Experiment script
# loads libraries
# loads in mean Fv/Fm data (n=3 per fragment) from 6 timepoints (days 0,1,2,5,10,11)
# checks assumptions using Shapiro Wilk for normality and Levene's Tests for homogeneity
# performs LME model hypothesis tests to identify significant differences in response vars due to fixed effects (trt within tp) while accounting for random effect of genet
# transforms data when necessary
# performs multiple comparison tests to identify different trts at each timepoint, when appropriate
# creates plots

#### load libraries
library(tibble)
library(ggplot2)
library(stats)
library(psych)
library(lme4)
library(car)
library(gtools)
library(plyr)
library(lattice)
library(DescTools)
library(data.table)
library(dplyr)
library(PMCMRplus)
library(onewaytests)
library(userfriendlyscience)
library(multcomp)
####################################

#### Set the working directory for all analyses & load the chapter 3 data
setwd("~/Documents/Documents - Harmony’s MacBook Pro - 1/Dissertation/Chapter3/analyses/revisionJune2022/")
fvfm=read.csv("StressMemory_allFvFm.csv", header=TRUE, stringsAsFactors=TRUE) # all values and timepoints for fv/fm

# reorder factors to appear as desired
fvfm$genet=factor(fvfm$genet,levels=c("A","B","C","D","E","F","G","H","I","J"))
fvfm$trt=factor(fvfm$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
fvfm$tp=factor(fvfm$tp, levels=c("1","2","3","4","5","6"), labels=c("day 0","day 1","day 2", "day 5", "day 10", "day 11"))

#### Examine the complete dataset
names(fvfm) # what variables are there
str(fvfm) # look at the data structure
headTail(fvfm) # look at the first few rows
summary(fvfm) # look for possible missing values or unbalanced design
colSums(is.na(fvfm))
# no NAs


#### Perform summary statistics for boxplots

fvfm_summary<- ddply(fvfm, c("day","trt"), summarise,
                         n = length(meanfvfm),
                         mean = mean(meanfvfm),
                         median = median(meanfvfm),
                         min = min(meanfvfm),
                         max = max(meanfvfm),
                         iqr = IQR(meanfvfm, na.rm=TRUE),
                         sd = sd(meanfvfm),
                         se = sd/ sqrt(n),
                         ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
fvfm_summary
fwrite(fvfm_summary,"summaryFvFm_x_Trt_Time.csv")

#### Parse data with into timepoints for plotting (or facet by tp)
fvfmtime1<-fvfm[fvfm$tp=="day 0",]; dim(fvfmtime1)
fvfmtime1$trt=factor(fvfmtime1$trt,levels=c("C","LH","HH"), labels=c("C","LH","HH"))
nor.test(meanfvfm~trt, data=fvfmtime1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time1<-homog.test(meanfvfm~trt, data=fvfmtime1, method=c("Levene"))

fvfmtime2<-fvfm[fvfm$tp=="day 1",]; dim(fvfmtime2)
fvfmtime2$trt=factor(fvfmtime2$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
nor.test(meanfvfm~trt, data=fvfmtime2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time2<-homog.test(meanfvfm~trt, data=fvfmtime2, method=c("Levene"))

fvfmtime3<-fvfm[fvfm$tp=="day 2",]; dim(fvfmtime3)
fvfmtime3$trt=factor(fvfmtime3$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
nor.test(meanfvfm~trt, data=fvfmtime3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time3<-homog.test(meanfvfm~trt, data=fvfmtime3, method=c("Levene"))

fvfmtime4<-fvfm[fvfm$tp=="day 5",]; dim(fvfmtime4)
fvfmtime4$trt=factor(fvfmtime4$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
nor.test(meanfvfm~trt, data=fvfmtime4, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time4<-homog.test(meanfvfm~trt, data=fvfmtime4, method=c("Levene"))

fvfmtime5<-fvfm[fvfm$tp=="day 10",]; dim(fvfmtime5)
fvfmtime5$trt=factor(fvfmtime5$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
nor.test(meanfvfm~trt, data=fvfmtime5, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time5<-homog.test(meanfvfm~trt, data=fvfmtime5, method=c("Levene"))

fvfmtime6<-fvfm[fvfm$tp=="day 11",]; dim(fvfmtime6)
fvfmtime6$trt=factor(fvfmtime6$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
nor.test(meanfvfm~trt, data=fvfmtime6, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
levene.time6<-homog.test(meanfvfm~trt, data=fvfmtime6, method=c("Levene"))

Mydotplot <- function(DataSelected){
  
  P <- dotplot(as.matrix(as.matrix(DataSelected)),
               groups=FALSE,
               strip = strip.custom(bg = 'white',
                                    par.strip.text = list(cex = 1.2)),
               scales = list(x = list(relation = "free", draw = TRUE),
                             y = list(relation = "free", draw = FALSE)),
               col=1, cex  = 0.5, pch = 16,
               xlab = list(label = "Value of the variable", cex = 1.5),
               ylab = list(label = "Order of the data in the file", cex = 1.5))
  print(P) 
}

MyVar <- c("meanfvfm")

# TIME 1
Mydotplot(fvfmtime1[,MyVar])
Mydotplot(fvfmtime2[,MyVar])
Mydotplot(fvfmtime3[,MyVar])
Mydotplot(fvfmtime4[,MyVar])
Mydotplot(fvfmtime5[,MyVar])
Mydotplot(fvfmtime6[,MyVar])



par(mfrow=c(1,1))
boxplot(time1$cells ~ time1$trt)
boxplot(time1$chl ~ time1$trt)
boxplot(time1$chlpcell ~ time1$trt)
boxplot(time1$prot ~ time1$trt)
boxplot(time1$protpcell ~ time1$trt)



#### Plot the complete dataset
black.bold.text<-element_text(face="bold",color="black", size=14)
black.italic.text<-element_text(face="bold.italic",color="black", size=18)
italic.text<-element_text(face="italic",color="black", size=14)
plain.text<-element_text(face="plain",color="black", size=12)

# create a mean for fvfm from all three trt samples on day 0
day0fvfm <- ddply(fvfm, "day", summarise,
                  n = length(meanfvfm),
                  mean = mean(meanfvfm),
                  median = median(meanfvfm),
                  min = min(meanfvfm),
                  max = max(meanfvfm),
                  iqr = IQR(meanfvfm, na.rm=TRUE),
                  sd = sd(meanfvfm),
                  se = sd/ sqrt(n),
                  ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
day0fvfm

day0fvfm = day0fvfm[1,] #take only the top row
day0fvfm<-add_column(day0fvfm, trt = "C", .after = "day")
fvfm_tps = fvfm_summary[4:29,]
fvfm_alltp_summary = rbind(day0fvfm, fvfm_tps)
fwrite(fvfm_alltp_summary, "fvfmSummary.csv") 

## Linear mixed-effects models

## Day 0 - Timepoint 1
par(mfrow = c(1,1))
fvfm.tp1 <- lmer(meanfvfm ~ trt + (1|genet), data=fvfmtime1, REML = TRUE)
summary(fvfm.tp1)
boxplot(meanfvfm ~ trt, data=fvfmtime1)
summary(glht(fvfm.tp1, linfct=mcp(trt="Tukey")))	# Use summary just to get differences and CIs
j <-Anova(fvfm.tp1,test.statistic = "F")
AovTbl<-anova(fvfm.tp1) # produces type I sum of squares
adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
AovTbl$Holm<-adj #this is adding adjusted pvalues
AovTbl
summary(fvfm.tp1) # Produces overall p-value, parameter estimates
plot(fvfm.tp1)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp1)) # histogram of residuals which should be normal
plot(fitted(fvfm.tp1), residuals(fvfm.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
## NS - lump together as controls and all starting values for all corals!

VarCorr(fvfm.tp1) # low variance by genet likely because all same zoox

p<-ggplot(fvfmtime1, aes(x=primeTank, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p + labs(x="Tank", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 0")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

## Day 1 - Timepoint 2
par(mfrow = c(1,1))
fvfm.tp2 <- lmer(meanfvfm ~ trt + (1|genet), data=fvfmtime2, REML = TRUE)
summary(fvfm.tp2)
boxplot(meanfvfm ~ trt, data=fvfmtime2)
summary(glht(fvfm.tp2, linfct=mcp(trt="Tukey")))
j <-Anova(fvfm.tp2,test.statistic = "F")
AovTbl<-anova(fvfm.tp2) 
adj<-p.adjust(j[,4],method="holm") 
AovTbl$Holm<-adj
AovTbl
summary(fvfm.tp2)
plot(fvfm.tp2)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp2))
plot(fitted(fvfm.tp2), residuals(fvfm.tp2))
## NS

VarCorr(fvfm.tp2) # more variance here possible attributable to overall thermal performance?

p<-ggplot(fvfmtime2, aes(x=trt, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p + labs(x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 1")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

## Day 2 - Timepoint 3
par(mfrow = c(1,1))
fvfm.tp3 <- lmer(meanfvfm ~ trt + (1|genet), data=fvfmtime3, REML = TRUE)
summary(fvfm.tp3)
boxplot(meanfvfm ~ trt, data=fvfmtime3)
summary(glht(fvfm.tp3, linfct=mcp(trt="Tukey")))
j <-Anova(fvfm.tp3,test.statistic = "F")
AovTbl<-anova(fvfm.tp3) 
adj<-p.adjust(j[,4],method="holm")
AovTbl$Holm<-adj 
AovTbl
summary(fvfm.tp3) 
plot(fvfm.tp3)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp3)) 
plot(fitted(fvfm.tp3), residuals(fvfm.tp3))
## NS 

VarCorr(fvfm.tp3)

p<-ggplot(fvfmtime3, aes(x=trt, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p + labs(x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 2")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

## Day 5 - Timepoint 4
par(mfrow=c(1,2))
hist(fvfmtime4$meanfvfm); qqnorm(fvfmtime4$meanfvfm)
fvfmtime4$cubemeanfvfm <- fvfmtime4$meanfvfm^3 # data have been cube transformed to correct left skewness
par(mfrow=c(1,2))
hist(fvfmtime4$cubemeanfvfm); qqnorm(fvfmtime4$cubemeanfvfm)
par(mfrow = c(1,1))
fvfm.tp4 <- lmer(cubemeanfvfm ~ trt + (1|genet), data=fvfmtime4, REML = TRUE)
summary(fvfm.tp4)
boxplot(meanfvfm ~ trt, data=fvfmtime4)
summary(glht(fvfm.tp4, linfct=mcp(trt="Tukey")))
j <-Anova(fvfm.tp4,test.statistic = "F")
AovTbl<-anova(fvfm.tp4)
adj<-p.adjust(j[,4],method="holm") 
AovTbl$Holm<-adj 
AovTbl
summary(fvfm.tp4)
plot(fvfm.tp4)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp4)) 
plot(fitted(fvfm.tp4), residuals(fvfm.tp4))

# multiple comparisons: trt is greater than the control
dunnett.fvfm.time4<-dunnettTest(meanfvfm ~ trt, fvfmtime4, alternative = "greater")
dt<-DunnettTest(meanfvfm ~ trt, fvfmtime4) # ignore p-values, they are not adjusted!
dt = data.frame(dt$`C`) # convert dt to a dataframe
pval.adj<- p.adjust(dunnett.fvfm.time4$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
dt$pval.adj <- pval.adj # add it to the dt df
dt	

VarCorr(fvfm.tp4)

p<-ggplot(fvfmtime4, aes(x=trt, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p + labs(x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 5")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

## Day 10 - Timepoint 5
library(LambertW) # https://stats.stackexchange.com/questions/85687/how-to-transform-leptokurtic-distribution-to-normality
citation("LambertW")
par(mfrow=c(1,2))
hist(fvfmtime5$meanfvfm); qqnorm(fvfmtime5$meanfvfm)
mod.Lh<- MLE_LambertW(fvfmtime5$meanfvfm, distname = "normal", type = "s") 
summary(mod.Lh)
fvfmtime5$probintmeanfvfm <- get_input(mod.Lh) # data have been probability integral transformed to correct Leptokurtic distribution to meet assumptions
par(mfrow=c(1,2))
hist(fvfmtime5$probintmeanfvfm); qqnorm(fvfmtime5$probintmeanfvfm)

par(mfrow = c(1,1))
fvfm.tp5 <- lmer(probintmeanfvfm ~ trt + (1|genet), data=fvfmtime5, REML = TRUE)
summary(fvfm.tp5)
boxplot(meanfvfm ~ trt, data=fvfmtime5)
summary(glht(fvfm.tp5, linfct=mcp(trt="Tukey")))
j <-Anova(fvfm.tp5,test.statistic = "F")
AovTbl<-anova(fvfm.tp5)
adj<-p.adjust(j[,4],method="holm")
AovTbl$Holm<-adj 
AovTbl
summary(fvfm.tp5)
plot(fvfm.tp5)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp5)) 
plot(fitted(fvfm.tp5), residuals(fvfm.tp5))

# multiple comparisons: trt is less than the control
dunnett.fvfm.time5<-dunnettTest(meanfvfm ~ trt, fvfmtime5, alternative=c("less"))
dt<-DunnettTest(meanfvfm ~ trt, fvfmtime5) 
dt = data.frame(dt$`C`) 
pval.adj<- p.adjust(dunnett.fvfm.time5$p.value,method="holm"); 
dt$pval.adj <- pval.adj 
dt

p<-ggplot(fvfmtime5, aes(x=trt, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  ylim(0,.7) +
  theme(legend.position="none")
p.labs <- p + labs(x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 10")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

## Day 11 - Timepoint 6

par(mfrow=c(1,2))
hist(fvfmtime6$meanfvfm); qqnorm(fvfmtime6$meanfvfm)
mod.Lh<- MLE_LambertW(fvfmtime6$meanfvfm, distname = "normal", type = "s") 
summary(mod.Lh)
fvfmtime6$probintmeanfvfm <- get_input(mod.Lh) # data have been probability integral transformed to correct Leptokurtic distribution to meet assumptions
par(mfrow=c(1,2))
hist(fvfmtime6$probintmeanfvfm); qqnorm(fvfmtime6$probintmeanfvfm)

par(mfrow = c(1,1))
fvfm.tp6 <- lmer(probintmeanfvfm ~ trt + (1|genet), data=fvfmtime6, REML = TRUE)
summary(fvfm.tp6)
boxplot(meanfvfm ~ trt, data=fvfmtime6)
summary(glht(fvfm.tp6, linfct=mcp(trt="Tukey")))
j <-Anova(fvfm.tp6,test.statistic = "F")
AovTbl<-anova(fvfm.tp6) 
adj<-p.adjust(j[,4],method="holm") 
AovTbl$Holm<-adj 
AovTbl
summary(fvfm.tp6)
plot(fvfm.tp6)
# check assumptions of the model
par(mfrow=c(1,2))
hist (residuals(fvfm.tp6)) 
plot(fitted(fvfm.tp6), residuals(fvfm.tp6))

# multiple comparisons: trt is less than the control
dunnett.fvfm.time6<-dunnTest(meanfvfm ~ trt, fvfmtime6)
dt<-DunnettTest(meanfvfm ~ trt, fvfmtime6) 
dt = data.frame(dt$`C`) 
pval.adj<- p.adjust(dunnett.fvfm.time6$p.value,method="holm");
dt$pval.adj <- pval.adj
dt
summary(glht(fvfm.tp6, linfct=mcp(trt="Tukey")), test = adjusted("holm"))


VarCorr(fvfm.tp6)

p<-ggplot(fvfmtime6, aes(x=trt, y=meanfvfm)) + 
  geom_boxplot(width=.5, lwd=1) +
  theme_bw() +
  ylim(0,.7) +
  theme(legend.position="none")
p.labs <- p + labs(x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), subtitle="Day 11")
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

###### PLOTS of FVFM SERIES

fvfmsummod=read.csv("fvfmSummary_modified.csv", header=TRUE, stringsAsFactors=TRUE) # controls on day 0 lumped together
names(fvfmsummod)
head(fvfmsummod)
# reorder factors to appear as desired
fvfmsummod$trt=factor(fvfmsummod$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

fvfmsummod2=read.csv("fvfmSummary_modified2.csv", header=TRUE, stringsAsFactors=TRUE) # controls and naives on all days lumped together for a faceted plot of control and naive on one panel and then other groups faceted?
names(fvfmsummod2)
head(fvfmsummod2)

# reorder factors to appear as desired
#fvfmsummod2$trt=factor(fvfmsummod2$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("Control","Naïve","LL (1d at 28°C)","LH (1d at 30.5°C)","HL (2d at 28°C","HH (2d at 30.5°C)"))
fvfmsummod2$trt=factor(fvfmsummod2$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
controls<- fvfmsummod2[fvfmsummod2$trt == c('C','N'),] 
treatments<- fvfmsummod2[!fvfmsummod2$trt == c('C','HH'),] 
hh<-fvfmsummod2[fvfmsummod2$trt == 'HH',] 
fvfmsummod2noHH<- fvfmsummod2[!fvfmsummod2$trt == 'HH',] 

black.italic.text<- element_text(family="Arial", face="bold.italic", color="black", size=24)
black.bold.text<- element_text(family="Arial", face="bold", color="black", size=18)
italic.text<- element_text(family="Arial", face="italic",color="black", size=18)
plain.text<- element_text(family="Arial", face="plain",color="black", size=18)
cbPalette <- c("#000000", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00") #color blindness friendly palette

### Fv/Fm across entire timeseries
p<-ggplot(fvfmsummod2, aes(y=mean, x=day, color=trt)) + 
  geom_line(lwd=.5,position=position_dodge(width=.7)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(width=.7), width=0)  +
  geom_point(size=3,position=position_dodge(width=.7))  +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) + 
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank())
  p.labs <- p + labs(x="Time (days)", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))  )
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

# specify new dfs for altering plots to reveal controls as 1 grouping on Day 0 or as Naives
fvfm<-read.csv("fvfm_df_modified_for_plotting_v2.csv", header=TRUE, stringsAsFactors = TRUE)
fvfm$trt=factor(fvfm$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
fvfm$new.trt.C=factor(fvfm$new.trt.C,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
fvfm$new.trt.N=factor(fvfm$new.trt.N,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
fvfmsummod$trt=factor(fvfmsummod$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

p<-ggplot(fvfm, aes(y=meanfvfm, x=day, fill = new.trt.C)) + 
  geom_boxplot(aes(group = interaction(day, new.trt.C)), size=.5,position=position_dodge(width=.8))  +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank())
p.labs <- p + labs(x="Time (days)", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))  )
plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

controls<-read.csv("controlNaive_day10-11.csv", header=TRUE, stringsAsFactors = TRUE)
controls$trt=factor(controls$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
controls$new.trt.C=factor(controls$new.trt.C,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
controls$new.trt.N=factor(controls$new.trt.N,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

## plots of controls and naives on days 10-11
p<-ggplot(controls, aes(x=day, y=meanfvfm)) +
  geom_line(data=fvfmsummod2[c(25,26,31,32),], aes(y=median, x=day, color=trt),lwd=1) +
  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) + #,position=position_dodge(width=.8))  +
  scale_x_continuous(breaks=c(9,10,11,12)) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position = "right")
p.labs <- p +  labs(x = "Time (days)", y = expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))) )
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
