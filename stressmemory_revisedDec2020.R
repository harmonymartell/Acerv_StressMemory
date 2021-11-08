## USAGE: "Thermal Priming Costs Outweigh Benefits in the Staghorn Coral Acropora cervicornis (Lamarck 1816)

# Stress Memory Experiment script
# loads libraries
# loads in SM data (all timepoints and values)
# creates boxplots for symbiont density, total chl, chl per cell, total algal protein and protein per cell
# searches for and identifies outliers
# checks assumptions using Shapiro Wilk for normality and Levene's Tests for Equal Variances
# performs hypothesis tests (i.e., ANOVAs and non-parametric tests) to identify significant differences
# performs multiple comparison tests to identify different trts at each timepoint, when appropriate

#### load libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("multtest")
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

#### Perform summary statistics for boxplots

	cells<- ddply(data, c("tp","trt"), summarise,
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

	chl<- ddply(data, c("tp","trt"), summarise,
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

	chlpcell<- ddply(data, c("tp","trt"), summarise,
			n = length(chlpcell),
			mean = mean(chlpcell),
			median = median(chlpcell),
			min = min(chlpcell),
			max = max(chlpcell),
			iqr = IQR(chlpcell, na.rm=TRUE),
			sd = sd(chlpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	chlpcell
fwrite(chlpcell,"summaryChlpCell_x_Trt_Time.csv")

	prot<- ddply(data, c("tp","trt"), summarise,
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

	protpcell<- ddply(data, c("tp","trt"), summarise,
			n = length(protpcell),
			mean = mean(protpcell),
			median = median(protpcell),
			min = min(protpcell),
			max = max(protpcell),
			iqr = IQR(protpcell, na.rm=TRUE),
			sd = sd(protpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	protpcell
fwrite(protpcell,"summaryProtpCell_x_Trt_Time.csv")

	fvfm<- ddply(data, c("tp","trt"), summarise,
			n = length(fvfm),
			mean = mean(fvfm),
			median = median(fvfm),
			min = min(fvfm),
			max = max(fvfm),
			iqr = IQR(fvfm, na.rm=TRUE),
			sd = sd(fvfm),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
			)
	fvfm
fwrite(fvfm,"summaryFvFm_x_Trt_Time.csv")


#### Parse data with into timepoints for plotting

time1<-data[data$tp=="1",]; dim(time1)
time2<-data[data$tp=="2",]; dim(time2)
time3<-data[data$tp=="3",]; dim(time3)

#### Plot the complete dataset

	black.bold.text<-element_text(face="bold",color="black", size=14)
	black.italic.text<-element_text(face="bold.italic",color="black", size=18)
	italic.text<-element_text(face="italic",color="black", size=14)
	plain.text<-element_text(face="plain",color="black", size=14)

	### Cells
	p<-ggplot(time1, aes(x=trt, y=cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y= cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(alpha=0.3, size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Time 3 cell stats
	time3stats <- ddply(time3, "trt", summarise,
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
	time3stats

	p<-ggplot(time3stats, aes(x=trt, y= median)) + 
		#geom_boxplot(width=.5, lwd=1) +
		geom_bar(stat="identity") +
		geom_errorbar(aes(ymin=median-se,ymax=median+se), width=.2) +
		#geom_point(alpha=0.3, size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Chl
	p<-ggplot(time1, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Chl per Cell
	p<-ggplot(time1, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,28) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,28) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,28) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	### Protein
	p<-ggplot(time1, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~ mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	### Protein per Cell
	p<-ggplot(time1, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  #geom_point(size=3, shape=21) +
	  ylim(0,0.6) +
	  theme_bw() +
	  theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(time2, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  #geom_point(size=3, shape=21) +
	  ylim(0,0.6) +
	  theme_bw() +
	  theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	p<-ggplot(time3, aes(x=trt, y=protpcell)) + 
	  geom_boxplot(width=.5, lwd=1) +
	  #geom_point(size=3, shape=21) +
	  ylim(0,0.6) +
	  theme_bw() +
	  theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
	
	### Fv/Fm
	p<-ggplot(time1, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time2, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

	p<-ggplot(time3, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))

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

	MyVar <- c("cells","chl","chlpcell","prot","protpcell","fvfm") 

	# TIME 1
	Mydotplot(time1[,MyVar])
	boxplot(time1$cells ~ time1$trt)
	boxplot(time1$chl ~ time1$trt)
	boxplot(time1$chlpcell ~ time1$trt)
	boxplot(time1$prot ~ time1$trt)
	boxplot(time1$protpcell ~ time1$trt)
	boxplot(time1$fvfm ~ time1$trt)

	# TIME 2
	Mydotplot(time2[,MyVar])
  boxplot(time2$cells ~ time2$trt)
  boxplot(time2$chl ~ time2$trt)
  boxplot(time2$chlpcell ~ time2$trt)
  boxplot(time2$prot ~ time2$trt)
  boxplot(time2$protpcell ~ time2$trt)
  boxplot(time2$fvfm ~ time2$trt)
  
	# TIME 3
	Mydotplot(time3[,MyVar])
	boxplot(time3$cells ~ time3$trt)
	boxplot(time3$chl ~ time3$trt)
	boxplot(time3$chlpcell ~ time3$trt)
	boxplot(time3$prot ~ time3$trt)
	boxplot(time3$protpcell ~ time3$trt)
	boxplot(time3$fvfm ~ time3$trt)

#### Assumption Checking

#### Shapiro Wilk Test for Normality
	### Time 1
	nor.test(cells~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(fvfm~trt, data=time1, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
	### Time 2
	nor.test(cells~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(fvfm~trt, data=time2, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
	### Time 3
	nor.test(cells~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chl~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(chlpcell~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(prot~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(protpcell~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	nor.test(fvfm~trt, data=time3, method=c('SW'), alpha=0.05, plot=c("qqplot-histogram"))
	
#### Levene's Test of Equal Variances
	### Time 1
	levene.cells1<-homog.test(cells~trt, data=time1, method=c("Levene"))
	levene.chl1<-homog.test(chl~trt, data=time1, method=c("Levene"))
	levene.chlpcell1<-homog.test(chlpcell~trt, data=time1, method=c("Levene"))
	levene.prot1<-homog.test(prot~trt, data=time1, method=c("Levene"))
	levene.protpcell1<-homog.test(protpcell~trt, data=time1, method=c("Levene"))
	levene.fvfm1<-homog.test(fvfm~trt, data=time1, method=c("Levene"))
	
	### Time 2
	levene.cells2<-homog.test(cells~trt, data=time2, method=c("Levene"))
	levene.chl2<-homog.test(chl~trt, data=time2, method=c("Levene"))
	levene.chlpcell2<-homog.test(chlpcell~trt, data=time2, method=c("Levene"))
	levene.prot2<-homog.test(prot~trt, data=time2, method=c("Levene"))
	levene.protpcell2<-homog.test(protpcell~trt, data=time2, method=c("Levene"))
	levene.fvfm2<-homog.test(fvfm~trt, data=time2, method=c("Levene"))
	
	### Time 3
	levene.cells3<-homog.test(cells~trt, data=time3, method=c("Levene"))
	levene.chl3<-homog.test(chl~trt, data=time3, method=c("Levene"))
	levene.chlpcell3<-homog.test(chlpcell~trt, data=time3, method=c("Levene"))
	levene.prot3<-homog.test(prot~trt, data=time3, method=c("Levene"))
	levene.protpcell3<-homog.test(protpcell~trt, data=time3, method=c("Levene"))
	levene.fvfm3<-homog.test(fvfm~trt, data=time3, method=c("Levene"))

#### Hypothesis Tests

	## TIME 1
    ### Cells - Time 1 with ANOVA
      aov <- aov(cells ~ trt, data=time1)
      summary(aov)
      boxplot(cells ~ trt, data=time1)
      summary(glht(aov, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
      model = lm (cells ~ trt, data=time1)
      Anova (model,type="II")
      AovTbl<-anova(model) # produces type I sum of squares
      AovTbl
      adj<-p.adjust(AovTbl[,5],method="holm") #adjusting p-values for multiple testing
      AovTbl$Holm<-adj #this is adding adjusted pvalues
      AovTbl
      summary(model) # Produces r-square, overall p-value, parameter estimates
      plot(model)
    # check assumptions of the model
      par(mfrow=c(1,2))
      hist (residuals(model)) # histogram of residuals which should be normal
      plot(fitted(model), residuals(model)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # multiple comparisons: trt is less than the control
      dunnett.cells.time1<-dunnettTest(cells ~ trt, time1, alternative=c("less"))
      p.adj<- p.adjust(dunnett.cells.time1$p.value,method="holm");
      t(p.adj)
      dunnett.cells.time1$padj <- t(p.adj)
      dt<-DunnettTest(cells ~ trt, time1) # ignore p-values, they are not adjusted!
      dunnett.cells.time1<- cbind(dunnett.cells.time1,dt[,1:3])
    # adjust the p values for FWER using Holm Step Down Method
      dunnett.cells.time1$p.adj <-p.adjust(dunnett.cells.time1$p.value,method="holm") 

      
    ### Chl - Time 1 with max-T + HCF3 Test
    amod <- aov(chl~trt, data = time1)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt="Dunnett"),vcov = vcovHC) # two-sided is default
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt = "Dunnett"))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    
    ### Chlpcell Time 1 with ANOVA
    aov <- aov(chlpcell ~ trt, data=time1)
    summary(aov)
    summary(glht(aov, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
    model = lm (chlpcell ~ trt, data=time1)
    Anova (model,type="II")
    anova(model) # produces type I sum of squares
    summary(model) # Produces r-square, overall p-value, parameter estimates
    plot(model)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(model)) # histogram of residuals which should be normal
    plot(fitted(model), residuals(model)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # multiple comparisons: trt is less than the control
    dunnett.chlpcell.time1<-dunnettTest(chlpcell ~ trt, time1, alternative=c("less"))
    p.adj<- p.adjust(dunnett.chlpcell.time1$p.value,method="holm");
    dunnett.chlpcell.time1$padj <- p.adj
    dt<-DunnettTest(chlpcell ~ trt, time1) # ignore p-values, they are not adjusted!
    #dunnett.cells.time1<- cbind(dunnett.cells.time1,dt[,1:3])
    # adjust the p values for FWER using Holm Step Down Method
    dunnett.chlpcell.time1$p.adj <-p.adjust(dunnett.chlpcell.time1$p.value,method="holm") 
    
    
    ### Protein - Time 1 with ANOVA
    aov <- aov(prot ~ trt, data=time1)
    summary(aov)
    summary(glht(aov, linfct=mcp(trt="Dunnett")))	# Dunnett - comparing all treatments with the control. Adjusted p values for the FWER.
    model = lm (prot ~ trt, data=time1)
    Anova (model,type="II")
    anova(model) # produces type I sum of squares
    summary(model) # Produces r-square, overall p-value, parameter estimates
    plot(model)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(model)) # histogram of residuals which should be normal
    plot(fitted(model), residuals(model)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # multiple comparisons: trt is less than the control
    dunnett.prot.time1<-dunnettTest(prot ~ trt, time1, alternative=c("less"))
    p.adj<- p.adjust(dunnett.prot.time1$p.value,method="holm");
    dunnett.prot.time1$padj <- p.adj
    dt<-DunnettTest(prot ~ trt, time1) # ignore p-values, they are not adjusted!
    #dunnett.cells.time1<- cbind(dunnett.cells.time1,dt[,1:3])
    # adjust the p values for FWER using Holm Step Down Method
    dunnett.chlpcell.time1$p.adj <-p.adjust(dunnett.chlpcell.time1$p.value,method="holm") 
    
    
    
    
    ### Protpcell - Time 1 with max-T + HCF3 Test
    amod <- aov(protpcell~trt, data = time1)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt="Dunnett"),vcov = vcovHC) # two-sided is default
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt = "Dunnett"))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    
    ### FvFm - Time 1 with ANOVA
    aov <- aov(fvfm ~ trt, data=time1)
    summary(aov)
  # no difference in trts
    model = lm (fvfm ~ trt, data=time1)
    Anova (model,type="II")
    anova(model) # produces type I sum of squares
    summary(model) # Produces r-square, overall p-value, parameter estimates
    plot(model)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(model)) # histogram of residuals which should be normal
    plot(fitted(model), residuals(model)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
### TIME 2 - (no sig relationships)
    
    ### Cells - Time 2 with max-T + HC3 Test
    amod <- aov(cells~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### Chl - Time 2 with max-T + HC3 Test  
    amod <- aov(chl~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### Chlpcell - Time 2 with max-T + HC3 Test
    amod <- aov(chlpcell~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.

    
    ### PROTEIN - Time 2 with max-T + HCF3 Test
    amod <- aov(prot~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### PROTPCELL - Time 2 with max-T + HCF3 Test
    amod <- aov(protpcell~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### FVFM - Time 2 with max-T + HCF3 Test
    amod <- aov(fvfm~trt, data = time2)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.

    
### TIME 3
    ### remove controls from the dataset here to rerun the dunnett's test relative to the naive trt
    ### Cells - Time 3 with max-T + HC3 Test
    amod <- aov(cells~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
      # are the trt values significantly different from the control?
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt="Dunnett"),vcov = vcovHC, alternative=c("less")) # two-sided is default
    
      # are the trt values signficantly different from the naive?
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt = c("C - N = 0",
                                                         "LL - N = 0",
                                                         "LH - N = 0",
                                                         "HL - N = 0",
                                                         "HH - N = 0")),
                              vcov = vcovHC)
      
    
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt = c("C - N = 0",
                                                           "LL - N = 0",
                                                           "LH - N = 0",
                                                           "HL - N = 0",
                                                           "HH - N = 0")))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    boxplot(cells~trt, data=time3)
    
    
    ### Chl - Time 3 with max-T + HC3 Test  
    amod <- aov(chl~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt = c("C - N = 0",
                                                         "LL - N = 0",
                                                         "LH - N = 0",
                                                         "HL - N = 0",
                                                         "HH - N = 0")),
                              vcov = vcovHC)
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt = c("C - N = 0",
                                                           "LL - N = 0",
                                                           "LH - N = 0",
                                                           "HL - N = 0",
                                                           "HH - N = 0")))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    
    
    ### Chlpcell - Time 3 with max-T + HC3 Test
    amod <- aov(chlpcell~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### PROTEIN - Time 3 with max-T + HCF3 Test
    amod <- aov(prot~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt = c("C - N = 0",
                                                         "LL - N = 0",
                                                         "LH - N = 0",
                                                         "HL - N = 0",
                                                         "HH - N = 0")),vcov = vcovHC) # two-sided is default
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt= c("C - N = 0",
                                                       "LL - N = 0",
                                                       "LH - N = 0",
                                                       "HL - N = 0",
                                                       "HH - N = 0")))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    
    ### PROTPCELL - Time 3 with max-T + HCF3 Test
    amod <- aov(protpcell~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    
    
    ### FVFM - Time 3 with max-T + HCF3 Test
    amod <- aov(fvfm~trt, data = time3)
    summary(amod)
    amod
    plot(amod)
    # check assumptions of the model
    par(mfrow=c(1,2))
    hist (residuals(amod)) # histogram of residuals which should be normal
    plot(fitted(amod), residuals(amod)) # A plot of residuals vs. predicted values. The residuals should be unbiased and homoscedastic.
    # Dunnett's following max-T HC3
    amod_glht_dunnett <- glht(amod, linfct = mcp(trt= c("C - N = 0",
                                                         "LL - N = 0",
                                                         "LH - N = 0",
                                                         "HL - N = 0",
                                                         "HH - N = 0")),vcov = vcovHC) # two-sided is default
    summary(amod_glht_dunnett)
    confint(amod_glht_dunnett) # confidence intervals for those
    # plot of the confidence intervals - set up so both plots will appear side by side
    #  ON THE LEFT IS THE PLOT WITH THE CORRECTION FOR HETEROSCEDASTCITY - ON THE RIGHT THE PLOT WHICH ASSUMES HOMOSECEDASTCITY
    # mai sets the bottom, left, top, right margins respectively in inches
    par (mfrow=c(1,2), cex=1,cex.lab=1,cex.axis=0.8,cex.main=1, mai=c(0.8,0.8,0.8,0.4))
    plot(confint(amod_glht_dunnett))
    # compare this to a model where homoscedasticity is assumed (note the omission of vcov=vcovHC to achieve this)
    amod_glht_dunnett_1 <- glht(amod, linfct = mcp(trt = "Dunnett"))
    summary(amod_glht_dunnett_1)
    confint(amod_glht_dunnett_1)
    plot(confint(amod_glht_dunnett_1))
    

### Cells over time
    ### Some summary stats
    
    chlpcellSM<- ddply(data, c("day","trt"), summarise,
                       n = length(chlpcell),
                       mean = mean(chlpcell),
                       median = median(chlpcell),
                       iqr = IQR(chlpcell, na.rm=TRUE),
                       sd = sd(chlpcell),
                       se = sd/ sqrt(n),
                       ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
    )
    chlpcellSM

    protpcellSM<- ddply(data, c("day","trt"), summarise,
                        n = length(protpcell),
                        mean = mean(protpcell),
                        median = median(protpcell),
                        iqr = IQR(protpcell, na.rm=TRUE),
                        sd = sd(protpcell),
                        se = sd/ sqrt(n),
                        ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
    )
    protpcellSM

    totalchlSM<- ddply(data, c("day","trt"), summarise,
                       n = length(chl),
                       mean = mean(chl),
                       median = median(chl),
                       iqr = IQR(chl, na.rm=TRUE),
                       sd = sd(chl),
                       se = sd/ sqrt(n),
                       ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
    )
    totalchlSM

    totalprotSM<- ddply(data, c("day","trt"), summarise,
                        n = length(prot),
                        mean = mean(prot),
                        median = median(prot),
                        iqr = IQR(prot, na.rm=TRUE),
                        sd = sd(prot),
                        se = sd/ sqrt(n),
                        ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
    )
    totalprotSM

    cellsSM<- ddply(data, c("day","trt"), summarise,
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
    cellsSM
    
    fvfmSM<- ddply(data, c("day","trt"), summarise,
                   n = length(fvfm),
                   mean = mean(fvfm),
                   median = median(fvfm),
                   iqr = IQR(fvfm, na.rm=TRUE),
                   sd = sd(fvfm),
                   se = sd/ sqrt(n),
                   ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
    )
    fvfmSM
    
# first, add a row of control values for time zero: the median Control values
# dummy_data <- data.frame(day= c(0,0,0,0,0,0), trt= c('C','N','LL','LH','HL','HH'), n= c(8,8,8,8,8,8), mean= c(3504030,3504030,3504030,3504030,3504030,3504030) , median= c(3504030,3504030,3504030,3504030,3504030,3504030), iqr= c(730466.8,730466.8,730466.8,730466.8,730466.8,730466.8), sd=c(482316.5, 482316.5, 482316.5, 482316.5, 482316.5, 482316.5), se=c(170524.64, 170524.64, 170524.64, 170524.64, 170524.64, 170524.64), ci=c(403226.7, 403226.7, 403226.7, 403226.7, 403226.7, 403226.7))
# headTail(dummy)

# plot the data
		black.italic.text<- element_text(family="Arial", face="bold.italic", color="black", size=28)
		black.bold.text<- element_text(family="Arial", face="bold", color="black", size=24)
		italic.text<- element_text(family="Arial", face="italic",color="black", size=24)
		plain.text<- element_text(family="Arial", face="plain",color="black", size=18)
#pdf(file="cells_trt_time.pdf")
	p<-ggplot(cellsSM, aes(x=day, y=median, colour=trt)) +
		geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
		geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
		geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
		#geom_ribbon(aes(fill=trt,ymin=min, ymax=max), alpha=0.05,position=position_dodge(width=0.6)) +
		scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
		#xlim(2,14) +
		theme_bw() +
		theme(legend.position=c(0.08,0.18), legend.title=element_blank())
	p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
#dev.off()


########### Now load all the Fv/Fm data from all timepoints
## set the working directory for all analyses
setwd("~/Desktop/PUBS/StressMemoryPaper/additionalAnalyses/")

####################################
## load the adjusted Chl data
adj_chl<-read.csv("~/Desktop/PUBS/stressmemorypaper/additionalAnalyses/adj_totalchlSM.csv", header=TRUE, stringsAsFactors=TRUE)
adj_chl$trt=factor(adj_chl$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
## load the adjusted cell data
adj_cells<-read.csv("~/Desktop/PUBS/stressmemorypaper/additionalAnalyses/adj_cellsSM.csv", header=TRUE, stringsAsFactors=TRUE)
adj_cells$trt=factor(adj_cells$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

## load the fvfm data
dat=read.csv("~/Desktop/PUBS/StressMemoryPaper/additionalAnalyses/allfvfmtp.csv", header=TRUE, stringsAsFactors=TRUE)
## reorder factors to appear as desired
dat$genet=factor(dat$genet,levels=c("C","D","E","F","G","H","I","J"))
dat$trt=factor(dat$trt,levels=c("C","U","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
#dat$day=factor(dat$day, levels=c("1","2","3","6","10","11")) #day is a numerical value, only do this for timepoint

## examine the complete dataset
names(dat) # what variables are there
str(dat) # look at the data structure
headTail(dat) # look at the first few rows
tail(dat) # look at the last few rows
summary(dat) # look for possible missing values or unbalanced design

# some summary stats
library(plyr)

fvfmStats <- ddply(dat, c("day","trt"), summarise,
                   n = length(fvfm),
                   mean = mean(round(fvfm,3)),
                   median = median(round(fvfm,3)),
                   min = min(round(fvfm,3)),
                   max = max(round(fvfm,3)),
                   iqr = IQR(round(fvfm,3), na.rm=TRUE),
                   sd = sd(round(fvfm,3)),
                   se = sd/ sqrt(n),
                   ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
fvfmStats
fvfmStats %>% mutate_if(is.numeric, round, 3)
fwrite(fvfmStats,"summaryFvFm_x_Trt_Timeseries.csv")


## This is to identify the outliers via Mahalanobis distance
library(stats)
library(reshape)
library(psych)

data_wide<-reshape(dat, idvar = "trt_genet", timevar = "day", direction = "wide")
headTail(data_wide)
fvfm<- select(data_wide, fvfm.0, fvfm.1, fvfm.2, fvfm.5, fvfm.10, fvfm.11); dim(fvfm); summary(fvfm)
mahal_fvfm = mahalanobis(fvfm, colMeans(fvfm, na.rm=TRUE), cov(fvfm,use="pairwise.complete.obs"));
mahal_fvfm
## determine a cutoff score
cutoff=qchisq(1-.001, ncol(fvfm))
cutoff
ncol(fvfm) # df
summary(mahal_fvfm < cutoff)  ## NAs
noout_fvfm = subset(fvfm, mahal_fvfm < cutoff) # remove samples with mahal values < the cutoff (> are outlying)
#noout_fvfm ## outliers #NONE  were removed from all timepoints
#THIS DF HAS NO NAs BECAUSE THERE ARE NO OUTLIERS!



####### This plots median +-1sd of cells, total chl and Fv/Fm over time in all trts #######

black.italic.text<- element_text(family="Arial", face="bold.italic", color="black", size=28)
black.bold.text<- element_text(family="Arial", face="bold", color="black", size=24)
italic.text<- element_text(family="Arial", face="italic",color="black", size=24)
plain.text<- element_text(family="Arial", face="plain",color="black", size=18)

# cells
p<-ggplot(cellsSM, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# total chl
p<-ggplot(totalchlSM, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# adjusted total chl (only day 10 and 11)
p<-ggplot(adj_chl, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.05)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.05)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.05)) +
  scale_x_continuous(breaks=c(9, 10,11,12)) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# adjusted cells (only day 10 and 11)
p<-ggplot(adj_cells, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.05)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.05)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.05)) +
  scale_x_continuous(breaks=c(9, 10,11,12)) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# fvfm with all timepoints
p<-ggplot(fvfmStats, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bolditalic(paste('Photochemical Efficiency (F'[V]*'/F'[M]*") " )))  )
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
