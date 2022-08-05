## USAGE: "Thermal Priming Costs Outweigh Benefits in the Staghorn Coral Acropora cervicornis (Lamarck 1816)"

# by Harmony Martell

# Stress Memory Experiment script
# loads libraries
# loads in SM data (all timepoints and values)
# searches for and identifies outliers via Cleveland dotplots
# checks assumptions using Shapiro-Wilk for normality and Levene's Tests for Equal Variances
# performs LME model hypothesis tests to identify significant differences in response vars due to fixed effects (trt within tp) while accounting for random effect of genet
# transforms data when necessary
# performs multiple comparison tests to identify different trts at each timepoint, when appropriate
# performs power analysis to determine effect sizes using Cohen's crtieria
# creates plots for symbiont density, total chl, chl per cell, total algal protein and protein per cell

#### load libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("multtest")
# library(multtest)
library(ggplot2)
library(psych)
library(lme4)
library(car)
library(plyr)
library(DescTools)
library(data.table)
library(lattice)
library(onewaytests)
library(PMCMRplus)
library(multcomp)
library(dplyr)
library(FSA)
library(LambertW) # https://stats.stackexchange.com/questions/85687/how-to-transform-leptokurtic-distribution-to-normality
#install.packages('Rcpp')


####################################

#### Set the working directory for all analyses & load the chapter 3 data
	setwd("~/Documents/Documents - Harmony’s MacBook Pro - 1/Dissertation/Chapter3/analyses/revisionJune2022/")
	data=read.csv("StressMemory_allSamples.csv", header=TRUE, stringsAsFactors=TRUE) # all values and timepoints
	
	# reorder factors to appear as desired
	data$genet=factor(data$genet,levels=c("C","D","E","F","G","H","I","J"))
	data$trt=factor(data$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	data$tp=factor(data$tp, levels=c("1","2","3"), labels=c("1","2","3"))

#### Examine the complete dataset
	names(data) # what variables are there
	str(data) # look at the data structure
	headTail(data) # look at the first few rows
	summary(data) # look for possible missing values or unbalanced design
	colSums(is.na(data))
	
	#### Parse data with into timepoints
	
	time1<-data[data$tp=="1",]; dim(time1)
	time1$trt=factor(time1$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
	time2<-data[data$tp=="2",]; dim(time2)
	time2$trt=factor(time2$trt,levels=c("C","LL","LH","HL","HH"), labels=c("C","LL","LH","HL","HH"))
	time3<-data[data$tp=="3",]; dim(time3)

	#### Identifying Outliers
	#### Function for multi-panel Cleveland dotplot. The input file must contain no categorical variables
	# code graciously provided by Richard P. Dunne 07 Sept 2020
	
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
	
	MyVar <- c("cells","chl","chlpcell","prot","protpcell") 
	
	# TIME 1
	Mydotplot(time1[,MyVar])
	par(mfrow=c(1,1))
	boxplot(time1$cells ~ time1$trt)
	boxplot(time1$chl ~ time1$trt)
	boxplot(time1$chlpcell ~ time1$trt)
	boxplot(time1$prot ~ time1$trt)
	boxplot(time1$protpcell ~ time1$trt)
	
	# TIME 2
	Mydotplot(time2[,MyVar])
	boxplot(time2$cells ~ time2$trt)
	boxplot(time2$chl ~ time2$trt)
	boxplot(time2$chlpcell ~ time2$trt)
	boxplot(time2$prot ~ time2$trt)
	boxplot(time2$protpcell ~ time2$trt)
	
	# TIME 3
	Mydotplot(time3[,MyVar])
	boxplot(time3$cells ~ time3$trt)
	boxplot(time3$chl ~ time3$trt)
	boxplot(time3$chlpcell ~ time3$trt)
	boxplot(time3$prot ~ time3$trt)
	boxplot(time3$protpcell ~ time3$trt)

	#### Assumption Checking
	#### Shapiro Wilk Test for Normality
	### Time 1
	nor.test(cells~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
	### Time 2
	nor.test(cells~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
	### Time 3
	nor.test(cells~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
	#### Levene's Test of Equal Variances
	### Time 1
	levene.cells1<-homog.test(cells~trt, data=time1, method=c("Levene"))
	levene.chl1<-homog.test(chl~trt, data=time1, method=c("Levene"))
	levene.chlpcell1<-homog.test(chlpcell~trt, data=time1, method=c("Levene"))
	levene.prot1<-homog.test(prot~trt, data=time1, method=c("Levene"))
	levene.protpcell1<-homog.test(protpcell~trt, data=time1, method=c("Levene"))
	
	### Time 2
	levene.cells2<-homog.test(cells~trt, data=time2, method=c("Levene"))
	levene.chl2<-homog.test(chl~trt, data=time2, method=c("Levene"))
	levene.chlpcell2<-homog.test(chlpcell~trt, data=time2, method=c("Levene"))
	levene.prot2<-homog.test(prot~trt, data=time2, method=c("Levene"))
	levene.protpcell2<-homog.test(protpcell~trt, data=time2, method=c("Levene"))
	
	### Time 3
	levene.cells3<-homog.test(cells~trt, data=time3, method=c("Levene"))
	levene.chl3<-homog.test(chl~trt, data=time3, method=c("Levene"))
	levene.chlpcell3<-homog.test(chlpcell~trt, data=time3, method=c("Levene"))
	levene.prot3<-homog.test(prot~trt, data=time3, method=c("Levene"))
	levene.protpcell3<-homog.test(protpcell~trt, data=time3, method=c("Levene"))
	
	#### Hypothesis Tests
	
	## TIME 1

### Cells - Time 1 (with LMEs)
	cells.tp1 <- lmer(cells ~ trt + (1|genet), data=time1, REML = TRUE)
	summary(cells.tp1)
	par(mfrow=c(1,1))
	boxplot(cells ~ trt, data=time1)
	summary(glht(cells.tp1, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(cells.tp1,test.statistic = "F"); j
	AovTbl<-anova(cells.tp1) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(cells.tp1)
	plot(cells.tp1)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(cells.tp1)) # histogram of residuals which should be normal
	plot(fitted(cells.tp1), residuals(cells.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.cells.time1<-dunnettTest(cells ~ trt, time1, alternative=c("less"))
	dt<-DunnettTest(cells ~ trt, time1) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.cells.time1$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	VarCorr(cells.tp1)
	
	
### Chl - Time 1
	chl.tp1 <- lmer(chl ~ trt + (1|genet), data=time1, REML = TRUE)
	summary(chl.tp1)
	boxplot(chl ~ trt, data=time1)
	summary(glht(chl.tp1, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(chl.tp1,test.statistic = "F")
	AovTbl<-anova(chl.tp1) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chl.tp1)
	plot(chl.tp1)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chl.tp1)) # histogram of residuals which should be normal
	plot(fitted(chl.tp1), residuals(chl.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.chl.time1<-dunnettTest(chl ~ trt, time1, alternative=c("less"))
	dt<-DunnettTest(chl ~ trt, time1) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.chl.time1$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	VarCorr(chl.tp1)
	
### Chlpcell Time 1 
	chlpcell.tp1 <- lmer(chlpcell ~ trt + (1|genet), data=time1, REML = TRUE)
	summary(chlpcell.tp1)
	boxplot(chlpcell ~ trt, data=time1)
	summary(glht(chlpcell.tp1, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(chlpcell.tp1,test.statistic = "F")
	AovTbl<-anova(chlpcell.tp1) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chlpcell.tp1) # Produces r-square, overall p-value, parameter estimates
	plot(chlpcell.tp1)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chlpcell.tp1)) # histogram of residuals which should be normal
	plot(fitted(chlpcell.tp1), residuals(chlpcell.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.chlpcell.time1<-dunnettTest(chlpcell ~ trt, time1)
	dt<-DunnettTest(chlpcell ~ trt, time1) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.chlpcell.time1$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	VarCorr(chlpcell.tp1)
	
	
### Protein - Time 1
	time1$logprot <- log(time1$prot) # data have been log transformed to correct right skew to meet assumptions
	time1$sqrtprot <- sqrt(time1$prot) # data have been sqrt transformed to correct right skew to meet assumptions
	time1$recipprot <- 1/time1$prot # data have been reciprocal transformed to correct right skew to meet assumptions
	prot.tp1 <- lmer(recipprot ~ trt + (1|genet), data=time1, REML = TRUE)
	summary(prot.tp1)
	par(mfrow=c(1,1))
	boxplot(prot ~ trt, data=time1)
	summary(glht(prot.tp1, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(prot.tp1,test.statistic = "F")
	AovTbl<-anova(prot.tp1) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(prot.tp1) # Produces r-square, overall p-value, parameter estimates
	plot(prot.tp1)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(prot.tp1)) # histogram of residuals which should be normal
	plot(fitted(prot.tp1), residuals(prot.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	## NS
	
	VarCorr(prot.tp1)
	
	### Protein per cell = Time 1
	protpcell.tp1 <- lmer(protpcell ~ trt + (1|genet), data=time1, REML = TRUE)
	summary(protpcell.tp1)
	boxplot(protpcell ~ trt, data=time1)
	summary(glht(protpcell.tp1, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(protpcell.tp1,test.statistic = "F")
	AovTbl<-anova(protpcell.tp1) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(protpcell.tp1) # Produces r-square, overall p-value, parameter estimates
	plot(protpcell.tp1)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(protpcell.tp1)) # histogram of residuals which should be normal
	plot(fitted(protpcell.tp1), residuals(protpcell.tp1)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is diff than the control
	dunnett.protpcell.time1<-dunnettTest(protpcell ~ trt, time1)
	dt<-DunnettTest(protpcell ~ trt, time1) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.protpcell.time1$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt

		### TIME 2 

	### Cells - Time 2
	cells.tp2 <- lmer(cells ~ trt + (1|genet), data=time2, REML = TRUE)
	summary(cells.tp2)
	boxplot(cells ~ trt, data=time2)
	summary(glht(cells.tp2, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(cells.tp2,test.statistic = "F")
	AovTbl<-anova(cells.tp2) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(cells.tp2) # Produces r-square, overall p-value, parameter estimates
	plot(cells.tp2)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(cells.tp2)) # histogram of residuals which should be normal
	plot(fitted(cells.tp2), residuals(cells.tp2)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.cells.time2<-dunnettTest(cells ~ trt, time2, alternative=c("less"))
	dt<-DunnettTest(cells ~ trt, time2) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.cells.time2$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	### Chl - Time 2
	chl.tp2 <- lmer(chl ~ trt + (1|genet), data=time2, REML = TRUE)
	summary(chl.tp2)
	boxplot(chl ~ trt, data=time2)
	summary(glht(chl.tp2, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(chl.tp2,test.statistic = "F")
	AovTbl<-anova(chl.tp2 ) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chl.tp2)
	plot(chl.tp2)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chl.tp2 )) # histogram of residuals which should be normal
	plot(fitted(chl.tp2 ), residuals(chl.tp2 )) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.chl.time2<-dunnettTest(chl ~ trt, time2, alternative=c("less"))
	dt<-DunnettTest(chl ~ trt, time2) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.chl.time2$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	### Chlpcell - Time 2
	time2<- na.omit(time2)
	par(mfrow=c(1,2))
	hist(time2$chlpcell); qqnorm(time2$chlpcell)
	mod.Lh<- MLE_LambertW(time2$chlpcell, distname = "normal", type = "s") 
	summary(mod.Lh)
	time2$probintchlpcell <- get_input(mod.Lh) # data have been probability integral transformed to correct Leptokurtic distribution to meet assumptions
	par(mfrow=c(1,2))
	hist(time2$probintchlpcell); qqnorm(time2$probintchlpcell)
	
	chlpcell.tp2 <- lmer(probintchlpcell ~ trt + (1|genet), data=time2, REML = TRUE)
	summary(chlpcell.tp2)
	boxplot(chlpcell ~ trt, data=time2)
	summary(glht(chlpcell.tp2, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(chlpcell.tp2,test.statistic = "F")
	AovTbl<-anova(chlpcell.tp2 ) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chlpcell.tp2)
	plot(chlpcell.tp2)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chlpcell.tp2 )) # histogram of residuals which should be normal
	plot(fitted(chlpcell.tp2 ), residuals(chlpcell.tp2 )) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.chlpcell.time2<-dunnettTest(chlpcell ~ trt, time2)
	dt<-DunnettTest(chlpcell ~ trt, time2) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.chlpcell.time2$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	### Protein - Time 2
	par(mfrow=c(1,2))
	hist(time2$prot); qqnorm(time2$prot)
	mod.Lh<- MLE_LambertW(time2$prot, distname = "normal", type = "s") 
	summary(mod.Lh)
	time2$probintprot <- get_input(mod.Lh) # data have been probability integral transformed to correct Leptokurtic distribution to meet assumptions
	par(mfrow=c(1,2))
	hist(time2$probintprot); qqnorm(time2$probintprot)
	time2$recipprot <- 1/(time2$prot+1) # data have been reciprocal transformed to correct right skew to meet assumptions
	par(mfrow=c(1,2))
	hist(time2$recipprot); qqnorm(time2$recipprot)
	
	prot.tp2 <- lmer(recipprot ~ trt + (1|genet), data=time2, REML = TRUE)
	summary(prot.tp2)
	boxplot(prot ~ trt, data=time2)
	summary(glht(prot.tp2, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(prot.tp2,test.statistic = "F")
	AovTbl<-anova(prot.tp2 ) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(prot.tp2) # Produces r-square, overall p-value, parameter estimates
	plot(prot.tp2)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(prot.tp2 )) # histogram of residuals which should be normal
	plot(fitted(prot.tp2 ), residuals(prot.tp2 )) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.prot.time2<-dunnettTest(prot ~ trt, time2)
	dt<-DunnettTest(prot ~ trt, time2) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.prot.time2$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	### Protein per cell - Time 2
	time2$sqrtprotpcell <- sqrt(time2$protpcell)
	par(mfrow=c(1,2))
	hist(time2$protpcell); hist(time2$sqrtprotpcell)
	
	protpcell.tp2 <- lmer(sqrtprotpcell ~ trt + (1|genet), data=time2, REML = TRUE)
	summary(protpcell.tp2)
	boxplot(protpcell ~ trt, data=time2)
	summary(glht(protpcell.tp2, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
	j <-Anova(protpcell.tp2,test.statistic = "F")
	AovTbl<-anova(protpcell.tp2 ) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(protpcell.tp2) # Produces r-square, overall p-value, parameter estimates
	plot(protpcell.tp2)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(protpcell.tp2 )) # histogram of residuals which should be normal
	plot(fitted(protpcell.tp2 ), residuals(protpcell.tp2 )) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	dunnett.protpcell.time2<-dunnettTest(protpcell ~ trt, time2)
	dt<-DunnettTest(protpcell ~ trt, time2) # ignore p-values, they are not adjusted!
	dt = data.frame(dt$`C`) # convert dt to a dataframe
	pval.adj<- p.adjust(dunnett.protpcell.time2$p.value,method="holm"); # do the p.adjustment using FWER with Holm Method
	dt$pval.adj <- pval.adj # add it to the dt df
	dt
	
	### TIME 3
	### Cells time 3
	cells.tp3 <- lmer(cells ~ trt + (1|genet), data=time3, REML = TRUE)
	summary(cells.tp3)
	boxplot(cells ~ trt, data=time3)
	j <-Anova(cells.tp3,test.statistic = "F"); j
	AovTbl<-anova(cells.tp3) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(cells.tp3) # Produces r-square, overall p-value, parameter estimates
	plot(cells.tp3)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(cells.tp3)) # histogram of residuals which should be normal
	plot(fitted(cells.tp3), residuals(cells.tp3)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: trt is less than the control
	summary(glht(cells.tp3, linfct=mcp(trt="Tukey")), test = adjusted("holm"))
	
	### Chl - Time 3
	par(mfrow=c(1,2))
	hist(time3$chl); qqnorm(time3$chl)
	mod.Lh<- MLE_LambertW(time3$chl, distname = "normal", type = "h") 
	summary(mod.Lh)
	time3$probintchl <- get_input(mod.Lh) # data have been probability integral transformed to correct Leptokurtic distribution to meet assumptions
	par(mfrow=c(1,2))
	hist(time3$probintchl); qqnorm(time3$probintchl)
	
	par(mfrow = c(1,1))
	chl.tp3 <- lmer(probintchl ~ trt + (1|genet), data=time3, REML = TRUE)
	summary(chl.tp3)
	boxplot(chl ~ trt, data=time3)
	j <-Anova(chl.tp3,test.statistic = "F")
	AovTbl<-anova(chl.tp3) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chl.tp3) # Produces r-square, overall p-value, parameter estimates
	plot(chl.tp3)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chl.tp3)) # histogram of residuals which should be normal
	plot(fitted(chl.tp3), residuals(chl.tp3)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: identify significantly different trts
	summary(glht(chl.tp3, linfct=mcp(trt="Tukey")), test = adjusted("holm"))
	
	### Chlpcell - Time 3
	par(mfrow = c(1,1))
	chlpcell.tp3 <- lmer(chlpcell ~ trt + (1|genet), data=time3, REML = TRUE)
	summary(chlpcell.tp3)
	boxplot(chlpcell ~ trt, data=time3)
  j <-Anova(chlpcell.tp3,test.statistic = "F"); j
	AovTbl<-anova(chlpcell.tp3) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(chlpcell.tp3) # Produces r-square, overall p-value, parameter estimates
	plot(chlpcell.tp3)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(chlpcell.tp3)) # histogram of residuals which should be normal
	plot(fitted(chlpcell.tp3), residuals(chlpcell.tp3)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# NS different
	
	### Prot - Time 3
	par(mfrow = c(1,1))
	prot.tp3 <- lmer(prot ~ trt + (1|genet), data=time3, REML = TRUE)
	summary(prot.tp3)
	boxplot(prot ~ trt, data=time3)
	j <-Anova(prot.tp3,test.statistic = "F"); j
	AovTbl<-anova(prot.tp3) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(prot.tp3) # Produces r-square, overall p-value, parameter estimates
	plot(prot.tp3)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(prot.tp3)) # histogram of residuals which should be normal
	plot(fitted(prot.tp3), residuals(prot.tp3)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: identify significantly different trts
	summary(glht(prot.tp3, linfct=mcp(trt="Tukey")), test = adjusted("holm"))
	
	### Protpcell - Time 3
	par(mfrow = c(1,1))
	protpcell.tp3 <- lmer(protpcell ~ trt + (1|genet), data=time3, REML = TRUE)
	summary(protpcell.tp3)
	boxplot(protpcell ~ trt, data=time3)
	j <-Anova(protpcell.tp3,test.statistic = "F"); j
	AovTbl<-anova(protpcell.tp3) # produces type I sum of squares
	adj<-p.adjust(j[,4],method="holm") #adjusting p-values for multiple testing
	AovTbl$Holm<-adj #this is adding adjusted pvalues
	AovTbl
	summary(protpcell.tp3) # Produces r-square, overall p-value, parameter estimates
	plot(protpcell.tp3)
	# check assumptions of the model
	par(mfrow=c(1,2))
	hist (residuals(protpcell.tp3)) # histogram of residuals which should be normal
	plot(fitted(protpcell.tp3), residuals(protpcell.tp3)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
	# multiple comparisons: identify significantly different trts
	summary(glht(protpcell.tp3, linfct=mcp(trt="Tukey")), test = adjusted("holm"))
	
##########################################################
	#### power analysis
	library(pwr)
	citation("pwr")
	pwr.anova.test(k = 5, n = 8, f = 0.25, sig.level = 0.05, power = NULL) # for general groupings at Timepoint 1 and 2
	pwr.anova.test(k = 6, n = 8, f = 0.25, sig.level = 0.05, power = NULL) # for general groupings at Timepoint 3
	
	j
	p.t.f2 <- pwr.f2.test(u = 5, v = 8, f2 = 0.25, sig.level = 0.05, power = NULL) # for lmes
	pwr.f2.test(u = 4, v = 33.964, f2 = 0.25, sig.level = 0.05, power = NULL) # for lmes
	pwr.f2.test(u = 5, v = 35, f2 = 0.25, sig.level = 0.05, power = NULL) # for lmes
	
	pwr.f2.test(u = 4, v = 36, f2 = NULL, sig.level = 0.05, power = 0.66) # effect sizes are greater with more power
	pwr.f2.test(u = 4, v = 33.964, f2 = NULL, sig.level = 0.05, power = 0.63) # for lmes
	pwr.f2.test(u = 5, v = 35, f2 = NULL, sig.level = 0.05, power = 0.61) # for lmes

	### Power was between 0.61-0.66 at moderate effect size (Cohen's f = 0.25) (considered a medium Cohen's effect)
	### With a power of 0.8, the effect sizes were increased to between 0.33-0.37 (considered a large Cohen's effects)
	### 
	cohen.ES(test = "f2", size = "small")
	cohen.ES(test = "f2", size = "medium")
	cohen.ES(test = "f2", size = "large")
	
	#### Perform summary statistics for boxplots

	cells<- ddply(data, c("day","trt"), summarise,
			n = length(cells),
			mean = mean(cells),
			median = median(cells),
			min = min(cells),
			max = max(cells),
			iqr = IQR(cells, na.rm=TRUE),
			sd = sd(cells),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	cells
fwrite(cells,"summaryCells_x_Trt_Time.csv")

	chl<- ddply(data, c("day","trt"), summarise,
			n = length(chl),
			mean = mean(chl),
			median = median(chl),
			min = min(chl),
			max = max(chl),
			iqr = IQR(chl, na.rm=TRUE),
			sd = sd(chl),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	chl
fwrite(chl,"summaryChl_x_Trt_Time.csv")

	chlpcell<- ddply(data, c("day","trt"),summarise,
			n = length(!is.na(chlpcell)),
			mean = mean(chlpcell, na.rm=TRUE),
			median = median(chlpcell, na.rm=TRUE),
			min = min(chlpcell, na.rm=TRUE),
			max = max(chlpcell, na.rm=TRUE),
			iqr = IQR(chlpcell, na.rm=TRUE),
			sd = sd(chlpcell, na.rm=TRUE),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	chlpcell
fwrite(chlpcell,"summaryChlpCell_x_Trt_Time.csv")

	prot<- ddply(data, c("day","trt"), summarise,
			n = length(prot),
			mean = mean(prot),
			median = median(prot),
			min = min(prot),
			max = max(prot),
			iqr = IQR(prot, na.rm=TRUE),
			sd = sd(prot),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	prot
fwrite(prot,"summaryProt_x_Trt_Time.csv")

	protpcell<- ddply(data, c("day","trt"), summarise,
			n = length(!is.na(protpcell)),
			mean = mean(protpcell, na.rm=TRUE),
			median = median(protpcell, na.rm=TRUE),
			min = min(protpcell, na.rm=TRUE),
			max = max(protpcell, na.rm=TRUE),
			iqr = IQR(protpcell, na.rm=TRUE),
			sd = sd(protpcell, na.rm=TRUE),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	protpcell
fwrite(protpcell,"summaryProtpCell_x_Trt_Time.csv")

##########################################################3
###PLOTTING

## set the controls 
controls<-read.csv("controlNaive_day10-11SM.csv", header=TRUE, stringsAsFactors = TRUE)
controls$trt=factor(controls$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
controls$new.trt.C=factor(controls$new.trt.C,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
controls$new.trt.N=factor(controls$new.trt.N,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

## Set the plot parameters	
  black.bold.text<-element_text(face="bold",color="black", size=14)
	black.italic.text<-element_text(face="bold.italic",color="black", size=18)
	italic.text<-element_text(face="italic",color="black", size=14)
	plain.text<-element_text(face="plain",color="black", size=14)

## Individual Boxplots by tp
		### Cells
	p<-ggplot(time1, aes(x=trt, y=cells/1000000)) + 
		geom_boxplot(width=.5, lwd=1) +
		ylim(0,5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('cells cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y= cells/1000000)) + 
		geom_boxplot(width=.5, lwd=1) +
		ylim(0,5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('cells cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=cells/1000000)) + 
		geom_boxplot(width=.5, lwd=1) +
		#ylim(0,5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('cells cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

		### Chl
	p<-ggplot(time1, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,3) +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~mu *'g Chl cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,8.5) +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~mu *'g Chl cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,7)
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~ mu *'g Chl cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Chl per Cell
	p<-ggplot(time1, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('pg Chl cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,3.5) +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('pg Chl cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('pg Chl cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	### Protein
	p<-ggplot(time1, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~mu *'g prot cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,0.7) +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~mu *'g prot cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		theme_bw() +
	  ylim(0,.7) +
		theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste( ~ mu *'g prot cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Protein per Cell
	p<-ggplot(time1, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  theme_bw() +
	  ylim(0,0.3) +
	  theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('ng prot cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(time2, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  theme_bw() +
	  theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('ng prot cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(time3, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  theme_bw() +
	  theme(legend.position="none")
	p.labs <- p + labs(x="", y=expression(bold(paste('ng prot cell'^-1*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Set plotting parameters
	black.italic.text<- element_text(face="bold.italic", color="black", size=24)
	black.bold.text<- element_text(face="bold", color="black", size=18)
	italic.text<- element_text(face="italic",color="black", size=18)
	plain.text<- element_text(face="plain",color="black", size=18)
	cbPalette <- c("#000000", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00") #color blindness friendly palette

	###### Control plots	
	
	cells <- read.csv("summaryCells_x_Trt_Time_modified.csv", header=TRUE, stringsAsFactors = TRUE)
	cells$trt=factor(cells$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	chls <- read.csv("summaryChl_x_Trt_Time_modified.csv", header=TRUE, stringsAsFactors = TRUE)
	chls$trt=factor(chls$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	chlpcells <- read.csv("summaryChlpcell_x_Trt_Time_modified.csv", header=TRUE, stringsAsFactors = TRUE)
	chlpcells$trt=factor(chlpcells$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	prots<- read.csv("summaryProt_x_Trt_Time_modified.csv", header=TRUE, stringsAsFactors = TRUE)
	prots$trt=factor(prots$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	protpcells<- read.csv("summaryProtpcell_x_Trt_Time_modified.csv", header=TRUE, stringsAsFactors = TRUE)
	protpcells$trt=factor(protpcells$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
	
	## plots of controls and naives on days 10-11
	p<-ggplot(controls, aes(x=day, y=cells/1000000)) +
	  geom_line(data=cells[c(7,8,13,14),], aes(y=median/1000000, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) +
	  scale_x_continuous(breaks=c(9,10,11,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste('cells cm'^-2*" "))))
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
	
	p<-ggplot(controls, aes(x=day, y=chl)) +
	  geom_line(data=chls[c(7,8,13,14),], aes(y=median, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) +
	  scale_x_continuous(breaks=c(9,10,11,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste( ~mu *'g Chl cm'^-2*" ")))) 
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
	
	p<-ggplot(controls, aes(x=day, y=prot)) +
	  geom_line(data=prots[c(7,8,13,14),], aes(y=median, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) +
	  scale_x_continuous(breaks=c(9,10,11,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste('ng Prot cm'^-2*" ")))) 
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
	
	p<-ggplot(controls, aes(x=day, y=chlpcell)) +
	  geom_line(data=chlpcells[c(7,8,13,14),], aes(y=median, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) +
	  scale_x_continuous(breaks=c(9,10,11,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  ylim(.8,2.5) +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste('pg Chl cell'^-1*" ")))) 
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
	
	p<-ggplot(controls, aes(x=day, y=protpcell)) +
	  geom_line(data=protpcells[c(7,8,13,14),], aes(y=median, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, new.trt.C), fill=new.trt.C)) +
	  scale_x_continuous(breaks=c(9,10,11,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste('ng Prot cell'^-1*" ")))) 
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
	
	### Cells across entire timeseries
	p<-ggplot(cells, aes(y=(mean/1000000), x=day, color=trt)) + 
	  geom_line(lwd=.5,position=position_dodge(width=.4)) +
	  geom_errorbar(aes(ymin=(mean/1000000)-(se/1000000), ymax=(mean/1000000)+(se/1000000)), width=0,position=position_dodge(width=.4))  +
	  geom_point(size=3,position=position_dodge(width=.4))  +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12), labels = c("0","2","4","6","8","10","12"), limits = c(0,12)) +
	  scale_y_continuous(breaks=c(1,2,3.0), labels=c("1.0","2.0","3.0")) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) + 
	  theme_bw() +
	  theme(legend.position=c(0.09,0.68), legend.title=element_blank())
	#theme(legend.position="none", legend.title=element_blank())
	p.labs <- p + labs(x="Time (days)", y=expression(bold(paste('cells cm'^-2*" "))))
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(chls, aes(y=mean, x=day, color=trt)) + 
	  geom_line(lwd=.5,position=position_dodge(width=.4)) +
	  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0,
	                position=position_dodge(.4)) +
	  geom_point(size=3,position=position_dodge(width=.4)) +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
	  scale_y_continuous(breaks=c(2,3,4,5), labels = c("2.0","3.0","4.0","5.0")) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) + 
	  theme_bw() +
	  theme(legend.position="none", legend.title=element_blank())
	p.labs <- p + labs(x="Time (days)",  y=expression(bold(paste( ~mu *'g Chl cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(prots, aes(y=mean, x=day, color=trt)) + 
	  geom_line(lwd=.5,position=position_dodge(width=.4)) +
	  geom_errorbar(aes(ymin=(mean)-(se), ymax=(mean)+(se)), width=0,
	                position=position_dodge(.4)) +
	  geom_point(size=3,position=position_dodge(width=.4)) +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
	  #scale_y_continuous(breaks=c(0.1,1.5,3)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) + 
	  theme_bw() +
	  theme(legend.position="none", legend.title=element_blank())
	p.labs <- p + labs(x="Time (days)",  y=expression(bold(paste( ~mu *'g Prot cm'^-2*" ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(chlpcells, aes(y=mean, x=day, color=trt)) + 
	  geom_line(lwd=.5,position=position_dodge(width=.4)) +
	  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0,
	                position=position_dodge(.4)) +
	  geom_point(size=3,position=position_dodge(width=.4)) +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
	  #scale_y_continuous(breaks=c(0.1,1.5,3)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) + 
	  theme_bw() +
	  theme(legend.position="none", legend.title=element_blank())
	p.labs <- p + labs(x="Time (days)",  y=expression(bold(paste('pg Chl cell'^-1*" "))))
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(protpcells, aes(y=mean, x=day, color=trt)) + 
	  geom_line(lwd=.5,position=position_dodge(width=.4)) +
	  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0,
	                position=position_dodge(.4)) +
	  geom_point(size=3,position=position_dodge(width=.4)) +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
	  scale_y_continuous(breaks=c(0,.1,.2)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) + 
	  theme_bw() +
	  theme(legend.position="none", legend.title=element_blank())
	p.labs <- p + labs(x="Time (days)",  y=expression(bold(paste('ng Prot cell'^-1*" "))))
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	
	## Boxplot of cells thru time
	p<-ggplot(data, aes(x=day, y=cells/1000000)) +
	  #geom_line(data=cells[c(1,2,7,8,13,14),], aes(y=median, x=day, color=trt),lwd=1) +
	  geom_boxplot(aes(group=interaction(day, trt), fill=trt)) + #,position=position_dodge(width=.8))  +
	  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
	  scale_fill_manual(values=cbPalette) +
	  scale_colour_manual(values=cbPalette) +
	  theme_bw() +
	  theme(legend.position = "none")
	p.labs <- p +  labs(x = "Time (days)", y=expression(bold(paste('cells cm'^-2*" ")))) 
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
