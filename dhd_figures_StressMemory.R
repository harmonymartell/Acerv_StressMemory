## USAGE: "Thermal Priming and Hormesis in the Staghorn Coral Acropora cervicornis (Lamarck 1816)
# written by Harmony Martell 11/2021

# loads libraries
# loads in SM data (all timepoints and values)
# plots stuff
# loads hormesis data
# makes figures for priming intensity vs. cells, chl post priming and post bleaching

library(multtest)
library(ggplot2)
library(stats)
library(reshape)
library(psych)
library(easyGgplot2)
library(devtools) 
library(lmmfit)
library(lme4)
library(labdsv)
library(vegan)
library(plotrix)
library(pgirmess)
library(gridExtra)
library(pbkrtest)
library(RVAideMemoire)
library(car)
library(gtools)
library(quantreg)
library(calibrate)
library(MASS)
library(e1071)
library(FSA)
library(desc)
library(plyr)
library(lattice)
library(DescTools)
library(data.table)
library(Hmisc)
library(dplyr)
library(PMCMRplus)
library(AICcmodavg)
library(nlme)
library(exactRankTests)
library(MCMCglmm)
library(onewaytests)
library(userfriendlyscience)
library(multcomp)
library(sandwich)
source("~/Documents/Documents - Harmony’s MacBook Pro - 1/ODU/BARSHIS_LAB_ALL/Chapter1_Recent Thermal History/Data/RsquaredGLMM.R")

####################################

#### Set the working directory for all analyses & load the chapter 3 data
setwd("~/Desktop/PUBS/CoralReefs_StressMemoryPaper/CoralReefsSubmission2/")
data=read.csv("SM2018_allSamples.csv", header=TRUE, stringsAsFactors=TRUE) # all values and timepoints

# reorder factors to appear as desired
data$genet=factor(data$genet,levels=c("A","B","C","D","E","F","G","H","I","J"))
data$trt=factor(data$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
data$tp=factor(data$tp, levels=c("1","2","3"), labels=c("1","2","3"))

#### Examine the complete dataset
names(data) # what variables are there
str(data) # look at the data structure
headTail(data) # look at the first few rows
summary(data) # look for possible missing values or unbalanced design
colSums(is.na(data))
# no NAs

#subset data for plotting
#### Parse data with into timepoints for plotting

time1<-data[data$tp=="1",]; dim(time1)
time2<-data[data$tp=="2",]; dim(time2)
time3<-data[data$tp=="3",]; dim(time3)

head(time1)
time1$trt[time1$trt=="N"]<-"C" # for Time 1, take mean of naive and controls combined to represent all heat stress free corals
summary(time1)

### Time 1 cell stats
time1stats <- ddply(time1, "trt", summarise,
                    n = length(cells),
                    min = min(cells),
                    max = max(cells),
                    mean = mean(cells),
                    median = median(cells),
                    iqr = IQR(cells, na.rm=TRUE),
                    sd = sd(cells),
                    se = sd/ sqrt(n),
                    ci = se * qt(.95/2 + .5, n-1), # Confidence at the 95% interval
)
time1stats
time1stats$cum.dhd = c(0.49, 2.31, 4.44, 5.11, 9.26)
time1stats

write.csv(time1stats, "cells_time1stats.csv")

###################
#LINEAR REGRESSION#
###################
# General equation: y = ax + b
### y is response variable (i.e., chl, zoox, Fv/Fm)
### x is the predictor variable (i.e., dhd)
### a and b are constants, the coefficients
# Specifics
### plot bleaching vs. dhd
# First, set x and y
x <- time1$dhd
y <- time1$cells
# Perform the regression
model <- lm(y ~ x)
# Examine regression summary details
model
print(summary(model))
# Examine the names of the summary output table
names(model)
names(summary(model))
AovTbl<-anova(model)
AovTbl
# Write out the Anova Table to a file
write.csv(AovTbl, file = "cells_v_dhdPriming.csv", na = " ")

library(ggplot2)

black.bold.text<-element_text(face="bold",color="black", size=18)
small.text<-element_text(face="bold",color="black", size=14)
legend.big.text<-element_text(face="italic",color="black", size=16)
legend.text<-element_text(face="bold",color="black", size=14)
plain.text<-element_text(face="plain",color="black", size=14)

p<-ggplot(data=time1stats, aes(x=cum.dhd, y=median, colour=trt, fill=trt)) +
  geom_point(size=5) +
  geom_abline(intercept=c(3605055), slope=(-146946), size=1) +
  geom_errorbar(aes(ymax=median+se, ymin=median-se), width=.5) +
  ylim(1500000,4000000) +
  theme_bw() +
  theme(legend.position="none") #c(.1,.2), legend.title=element_blank(), strip.text.x = element_text(size=14, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) 
p.labs <- p + labs(title=" ", x =expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y=expression(bold(paste('Symbiont Density (cells cm'^-2*") "))))
plot(p.labs + theme(axis.title= black.bold.text, axis.text=small.text, legend.text = legend.text))


### Chlorophyll

### Time 1 chl stats
time1stats <- ddply(time1, "trt", summarise,
                    n = length(chl),
                    min = min(chl),
                    max = max(chl),
                    mean = mean(chl),
                    median = median(chl),
                    iqr = IQR(chl, na.rm=TRUE),
                    sd = sd(chl),
                    se = sd/ sqrt(n),
                    ci = se * qt(.95/2 + .5, n-1), # Confidence at the 95% interval
)
time1stats
time1stats$cum.dhd = c(0.49, 2.31, 4.44, 5.11, 9.26)
time1stats

write.csv(time1stats, "chl_time1stats.csv")

###################
#LINEAR REGRESSION#
###################
# General equation: y = ax + b
### y is response variable (i.e., loss chl, loss zoox, loss Fv/Fm)
### x is the predictor variable (i.e., permutations of temperature at each site for each assay time)
### a and b are constants, the coefficients
# Specifics
### plot bleaching vs. temp permutation of 5 assays x 5 individuals (25 responses)
# First, set x and y
x <- time1$dhd
y <- time1$chl
# Perform the regression
model <- lm(y ~ x)
# Examine regression summary details
model
print(summary(model))
# Examine the names of the summary output table
names(model)
names(summary(model))
AovTbl<-anova(model)
AovTbl
# Write out the Anova Table to a file
write.csv(AovTbl, file = "chl_v_dhdPriming.csv", na = " ")

p<-ggplot(data=time1stats, aes(x=cum.dhd, y=median, colour=trt, fill=trt)) +
  geom_point(size=5) +
  geom_abline(intercept=c(2.08715), slope=(-0.06808), size=1) +
  geom_errorbar(aes(ymax=median+se, ymin=median-se), width=.5) +
  theme_bw() +
  theme(legend.position=c(.1,.2), legend.title=element_blank(), strip.text.x = element_text(size=14, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) 
p.labs <- p + labs(title=" ", x =expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y=expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") "))))
plot(p.labs + theme(axis.title= black.bold.text, axis.text=small.text, legend.text = legend.text))


#####################################################################################
hormesis=read.csv("hormesis.csv", header=TRUE, stringsAsFactors=TRUE) # only the timepoint 3 values
## This dataset has dhd.prime, which is the cumulative heat stress (including the field) up to the end of the priming exposure (Timepoint 1)

# reorder factors to appear as desired
hormesis$genet=factor(hormesis$genet,levels=c("A","B","C","D","E","F","G","H","I","J"))
hormesis$trt=factor(hormesis$trt,levels=c("N","LL","LH","HL","HH"), labels=c("N","LL","LH","HL","HH"))
hormesis$tp=factor(hormesis$tp, levels=c("1","2","3"), labels=c("1","2","3"))

#### Examine the complete dataset
names(hormesis) # what variables are there
str(hormesis) # look at the data structure
headTail(hormesis) # look at the first few rows
summary(hormesis) # look for possible missing values or unbalanced design
colSums(is.na(hormesis))

library(ggplot2)

black.bold.text<-element_text(face="bold",color="black", size=18)
small.text<-element_text(face="bold",color="black", size=14)
legend.big.text<-element_text(face="italic",color="black", size=16)
legend.text<-element_text(face="bold",color="black", size=14)
plain.text<-element_text(face="plain",color="black", size=14)

#Now plot mean with error bars in both x and y directions
### For Cells
hormesisStatsCells <- ddply(hormesis, "trt", summarise,
                            n = length(cells),
                            min = min(cells),
                            max = max(cells),
                            mean = mean(cells),
                            median = median(cells),
                            iqr = IQR(cells, na.rm=TRUE),
                            sd = sd(cells),
                            se = sd/ sqrt(n),
                            ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
hormesisStatsCells
hormesisStatsCells$prime.dhd <- tapply(hormesis$dhd.prime, hormesis$trt, mean)
hormesisStatsCells$cum.dhd <- tapply(hormesis$dhd.cum, hormesis$trt, mean)
hormesisStatsCells$cum.dhd.se <- tapply(hormesis$dhd.cum, hormesis$trt, se)
hormesisStatsCells

fwrite(hormesisStatsCells, file="cells_hormesisStats.csv")

### For Chl
hormesisStatsChl <- ddply(hormesis, "trt", summarise,
                            n = length(chl),
                            min = min(chl),
                            max = max(chl),
                            mean = mean(chl),
                            median = median(chl),
                            iqr = IQR(chl, na.rm=TRUE),
                            sd = sd(chl),
                            se = sd/ sqrt(n),
                            ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
hormesisStatsChl
hormesisStatsChl$prime.dhd <- tapply(hormesis$dhd.prime, hormesis$trt, mean)
hormesisStatsChl$cum.dhd <- tapply(hormesis$dhd.cum, hormesis$trt, mean)
hormesisStatsChl$cum.dhd.se <- tapply(hormesis$dhd.cum, hormesis$trt, se)
hormesisStatsChl

fwrite(hormesisStatsChl, file="chl_hormesisStats.csv")

##Plot the means with errorbars on both x and y

# First for Cells
p<- ggplot(hormesisStatsCells, aes(x=cum.dhd, y=median, colour = trt)) +
  stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
  geom_point(size=5) +
  geom_errorbar(aes(ymax=median+se, ymin=median-se), width=.5) +
  geom_errorbarh(aes(xmax=cum.dhd+cum.dhd.se, xmin=cum.dhd-cum.dhd.se), height=.5) +
  #ylim(.5e06, 3E06) +
  #xlim(0,12) +
  theme_bw() +
  theme(legend.position="none") #c(.07,.81), legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white"))
p.labs<- p + labs(title=" ", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y=expression(bold(paste('Symbiont Density (cells cm'^-2*") "))))
plot(p.labs + theme(axis.title= black.bold.text, axis.text=small.text, legend.text = legend.text))


### Now for Chl
p<- ggplot(hormesisStatsChl, aes(x=cum.dhd, y=median, colour = trt)) +
  stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
  geom_point(size=5) +
  geom_errorbar(aes(ymax=median+se, ymin=median-se), width=.5) +
  #ylim(.5e06, 3E06) +
  #xlim(0,12) +
  theme_bw() +
  theme(legend.position=c(.07,.18), legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white"))
p.labs<- p + labs(title=" ", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y= expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") "))))
plot(p.labs + theme(axis.title= black.bold.text, axis.text=small.text, legend.text = legend.text))

####################################################

########################
#POLYNOMIAL REGRESSIONS#
########################
# General equation: y = ax^2 - bx + c
# Cells = -44740x2 + 356237x + 2E+06
#    R² = 0.9229
# Chl = -0.037x2 + 0.2742x + 3.5987
#    R² = 0.9996

### y is response variable (i.e., sym density, chl)
### x is the predictor variable (i.e., cum.dhd)
### a-d are constants, the coefficients

# First, set x and y
x <- hormesis$dhd.prime
y <- hormesis$chl
#y <- hormesis$cells
plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 

# Perform the regression
model <- lm(y ~ I(x^2) + I(x))
# Examine regression summary details
model
print(summary(model))
# Examine the names of the summary output table
names(model)
names(summary(model))
AovTbl<-anova(model)
AovTbl
# Write out the Anova Table to a file
write.csv(AovTbl, file = "chl_v_dhdPriming.csv", na = " ")
#write.csv(AovTbl, file= "cells_v_dhdPriming.csv,na = " ")
