# Working Directory
setwd("~/Desktop/Invasive.Signal/Invasive.Signal")

# Trait Data
trait.data=read.csv("trait.data.port.csv", header=T)

# Determine if everything is spelled correctly

table(trait.data$Dispersal.Syn.Simplify)
# 21 Unknown
table(trait.data$Duration.Categories)
# Remove Biennial: only 1
table(trait.data$Status)
# 110 native, 60 introduced
table(trait.data$Growth.Habitat.Categories)
# Remove Subshrub: only 1
# Remove Vine: Only 3
table(trait.data$Pollination.Syn.Simplify)
# Remove Hydrophily: Only 1
# 17 Unknown

# Make contingency tables

table(trait.data$Status, trait.data$Dispersal.Syn.Simplify)
table(trait.data$Status, trait.data$Pollination.Syn.Simplify)
table(trait.data$Status, trait.data$Duration.Categories)
table(trait.data$Status, trait.data$Growth.Habitat.Categories)

# Read in tables

disperse.table=read.csv("disperse.table.csv", header=T, row.names = 1)
duration.table=read.csv("duration.table.csv", header=T, row.names = 1)
growth.table=read.csv("growth.table.csv", header=T, row.names = 1)
pollination.table=read.csv("pollination.table.csv", header=T, row.names = 1)

# graphing
library("graphics")
disperse.dt=as.table(as.matrix(disperse.table))
mosaicplot(disperse.dt, shade = TRUE, main="")
library(vcd)
assoc(disperse.dt, shade=TRUE)

# Chi-square test
# Significance means: row and column variables are significantly associated
# To know the most contributing cell to the total Chisq, calculate chisq stat for each cell
# Use formula for so-called Pearson residuals (r) for each cell or standardized residuals
# Cells with highest absolute standardized residuals contribute most to total Chisq score
# Positive residuals specify attraction between row and column variables
# Negative residuals specify repulsion between corresponding row and column
# Contribution in % of a given cell to total Chisq

# Perform post hoc analysis based on residuals of Pearson's Chisq test for Count Data
install.packages("chisq.posthoc.test")
library(chisq.posthoc.test)

disperse.chisq=chisq.test(disperse.table) # p = 0.511, total Chisq = 3.2868, df = 4
# No relationship between dispersal syndrom and status
disperse.chisq$residuals 
# negative association between Native and Zoochory and Polychory, but positive association between Introduced and Zoochory and Polychory
disperse.contrib=100*disperse.chisq$residuals^2/disperse.chisq$statistic
disperse.posthoc=chisq.posthoc.test(disperse.table)
# No post hoc analysis was significant

library(corrplot)
corrplot(disperse.chisq$residuals, is.cor = FALSE)
corrplot(disperse.contrib, is.cor = FALSE)


duration.chisq=chisq.test(duration.table) # p = 0.006504, total Chisq = 10.071, df = 2
# Significant relationship between duration and status
duration.chisq$residuals 
# negative association between Native and Annual, but positive association between Introduced and Annual
# positive association between Native and Perennial, but negative association between Introduced and Perennial
duration.contrib=100*duration.chisq$residuals^2/duration.chisq$statistic
duration.posthoc=chisq.posthoc.test(duration.table)
# Post Hoc Significant: Native Perennial, Introduced Perennial

corrplot(duration.chisq$residuals, is.cor = FALSE)
corrplot(duration.contrib, is.cor = FALSE)

growth.chisq=chisq.test(growth.table) # p = 0.8989, total Chisq = 0.5891, df = 3
# No relationship between growth and status
growth.chisq$residuals 
growth.contrib=100*growth.chisq$residuals^2/growth.chisq$statistic
growth.posthoc=chisq.posthoc.test(growth.table)
# Post Hoc Significant: Native Annual, Native Introduced, Native Perennial, Introduced Perennial

corrplot(duration.chisq$residuals, is.cor = FALSE)
corrplot(duration.contrib, is.cor = FALSE)


# Rerun Growth removing Tree category
growth.table.2=growth.table[1:3,]
growth.chisq.2=chisq.test(growth.table.2) # p = 0.8056, total Chisq = 0.43244, df = 2
# No relationship between growth and status
growth.chisq.2$residuals 
growth.contrib.2=100*growth.chisq.2$residuals^2/growth.chisq.2$statistic
growth.posthoc.2=chisq.posthoc.test(growth.table.2)

pollination.chisq=chisq.test(pollination.table) # p = 0.0001, total Chisq = 19.867, df = 3
# Significant Relationship between Pollination Syndrome and Status
pollination.chisq$residuals 
pollination.contrib=100*pollination.chisq$residuals^2/pollination.chisq$statistic
pollination.posthoc=chisq.posthoc.test(pollination.table)
# Post Hoc Significant: Native Multiple,Introduced Multiple, Native Zoophily, Introduced Zoophily

# Rerun Pollination removing Selfing category
pollination.table.2=pollination.table[c(1,2,4),]
pollination.chisq.2=chisq.test(pollination.table.2) # p = 0.001403, total Chisq = 13.139, df = 2
pollination.chisq.2$residuals 
pollination.contrib.2=100*pollination.chisq.2$residuals^2/pollination.chisq.2$statistic
pollination.posthoc.2=chisq.posthoc.test(pollination.table.2)


corrplot(pollination.chisq$residuals, is.cor = FALSE)
corrplot(pollination.contrib, is.cor = FALSE)


# Plotting data
# Plots include all data
# Read in tables for all data

disperse.table.all=read.csv("disperse.table.all.csv", header=T, row.names = 1)
duration.table.all=read.csv("duration.table.all.csv", header=T, row.names = 1)
growth.table.all=read.csv("growth.table.all.csv", header=T, row.names = 1)
pollination.table.all=read.csv("pollination.table.all.csv", header=T, row.names = 1)

disperse.dt.all=as.table(as.matrix(disperse.table.all))
duration.dt.all=as.table(as.matrix(duration.table.all))
growth.dt.all=as.table(as.matrix(growth.table.all))
pollination.dt.all=as.table(as.matrix(pollination.table.all))


# Dispersal
barplot(t(disperse.dt.all), beside = TRUE, xlab="Dispersal Syndrome", axis.lty = 1,
        ylab="Number of Species", col=c("black", "gray"), ylim=c(0, 35), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)
  
# Duration
barplot(t(duration.dt.all), beside = TRUE, xlab="Duration", axis.lty = 1,
        ylab="Number of Species", col=c("black", "gray"), ylim=c(0, 70), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Growth Habitat
barplot(t(growth.dt.all), beside = TRUE, xlab="Growth Habit", axis.lty = 1,
        ylab="Number of Species", col=c("black", "gray"), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Pollination
barplot(t(pollination.dt.all), beside = TRUE, xlab="Pollination Syndrome", axis.lty = 1,
        ylab="Number of Species", col=c("black", "gray"), ylim=c(0, 60),cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Remake BarGraph with numbers of species in each category as proportions.

colSums(disperse.dt.all) #110/610 = 170
disperse.proport=disperse.dt.all
disperse.proport=disperse.proport/170

colSums(duration.dt.all) #110/610 = 170
duration.proport=duration.dt.all
duration.proport=duration.proport/170

colSums(growth.dt.all) #110/610 = 170
growth.proport=growth.dt.all
growth.proport=growth.proport/170

colSums(pollination.dt.all) #110/610 = 170
pollination.proport=pollination.dt.all
pollination.proport=pollination.proport/170

# Dispersal
barplot(t(disperse.proport), beside = TRUE, xlab="Dispersal Syndrome", axis.lty = 1,
        ylab="Proportion of Species", col=c("black", "gray"), ylim=c(0,0.20), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Duration
barplot(t(duration.proport), beside = TRUE, xlab="Duration", axis.lty = 1,
        ylab="Proportion of Species", col=c("black", "gray"), ylim=c(0,0.35), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Growth Habitat
barplot(t(growth.proport), beside = TRUE, xlab="Growth Habit", axis.lty = 1,
        ylab="Proportion of Species", col=c("black", "gray"), ylim=c(0,0.35), cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

# Pollination
barplot(t(pollination.proport), beside = TRUE, xlab="Pollination Syndrome", axis.lty = 1,
        ylab="Proportion of Species", col=c("black", "gray"),ylim=c(0,0.35),cex.names=1.5,
        cex.axis = 1.5, cex.lab=1.5)

