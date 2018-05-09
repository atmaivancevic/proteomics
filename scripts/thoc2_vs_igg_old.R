
######### METHOD1 ###########

####################################################################################

# Set working dir and import/filter data

setwd("/Users/atma/Desktop/Raman/thoc2_old")

data = read.table("matrix", sep = "\t", header = T)
data

# remove line numbers 
dimnames(data)[[1]] <- data[,1]
data = data[,-1]
data

# check column names
colnames(data)

# remove rows which contain 4 or more 0 (too many missing values)
data <- data[rowSums(data == 0) <= 4, ]
data

# add 0.001 to each entry
data = data + 0.001
data

########################################################################################
#Since the data are essentially constant
#R will not compute pvalues
#Work-around is to use Excel to compute the pvalues

# export table to do paried t-test in Excel
write.csv(data, file="matrix_prelim.csv")

##########################################################################################

# Then import the table back in, with pvalue as an additional column

setwd("/Users/atma/Desktop/Raman/thoc2_old")

data = read.table("matrix_with_pvalue", sep = "\t", header = T)
data

# remove line numbers 
dimnames(data)[[1]] <- data[,1]
data = data[,-1]
data

# check column names
colnames(data)

rawpvalue <- data[,7]
rawpvalue
class(rawpvalue)

hist(rawpvalue)

#ttestData <- function(df, grp1, grp2) {
#  x = df[grp1]
#  y = df[grp2]
#  x = as.numeric(x)
#  y = as.numeric(y)  
#  results = t.test(x, y, paired=TRUE)
#  results$p.value
#}
#rawpvalue = apply(data, 1, ttestData, grp1 = c(1:3), grp2 = c(4:6))

########################################################################################

# Calculate fold change by log transform then subtracting control from test

#transform our data into log2 base.
data = data[,1:6]
data
data = log2(data)
data

############### METHOD 1

# First, take the means for each gene and subtract

#calculate the mean of each gene per control group
iggfirst = apply(data[,1:3], 1, mean)

uspfirst = apply(data[, 4:6], 1, mean) 

log2foldchange_firstmethod <- uspfirst - iggfirst 
log2foldchange_firstmethod

hist(log2foldchange_firstmethod, xlab = "log2 Fold Change (Test minus Control)")

########################################################################################

# Plot volcano

results = cbind(log2foldchange_firstmethod, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)
results

#svg('THOC2_old.svg')

# Try to make labels not overlap with each other
# separate data into 3 subsets
value1 = subset(results, rawpvalue<0.05 & abs(log2foldchange_firstmethod)<8)
value1

value2 = subset(results, rawpvalue>0.05 & abs(log2foldchange_firstmethod)>8)
value2

value3 = subset(results,rawpvalue<0.05 & abs(log2foldchange_firstmethod)>8)
value3

#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)

ggplot(results, aes(log2foldchange_firstmethod, -log10(rawpvalue)), colour="grey") +
  # Set all dots color to grey
  geom_point(
    data=results, colour = "grey") + 
  # If pvalue<0.05, change dot color to green
  geom_point(
    data=results[which(results$rawpvalue<0.05),], colour = "green") + 
  # If log2FC >1, change dot color to red
  geom_point(
    data=results[which(abs(results$log2foldchange_firstmethod)>8),], colour = "red") +
  # If both, change dot color to blue
  geom_point(
    data=results[which(abs(results$log2foldchange_firstmethod)>8 & results$rawpvalue<0.05),], colour = "blue") +
  # Add text if pvalue < 0.05
  geom_text_repel(
    data = value1,
    mapping = aes(log2foldchange_firstmethod, -log10(rawpvalue), label = rownames(value1), colour = 'pvalue<0.05'), 
    size = 2,
  ) +
  # Add text if log2FC > 1
  geom_text_repel(
    data = value2,
    mapping = aes(log2foldchange_firstmethod, -log10(rawpvalue), label = rownames(value2), colour = '|log2FC|>8'),
    size = 2,
  ) +
  # Add text if both
  geom_text_repel(
    data = value3,
    mapping = aes(log2foldchange_firstmethod, -log10(rawpvalue), label = rownames(value3), colour = 'pvalue<0.05 & |log2FC|>8'),
    fontface = 'bold',
    size = 4,
  ) +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("Volcano plot: THOC2 vs IgG (old)") +
  labs(y="-log10(Pvalue)", x = "log2(FoldChange)") +
  xlim(-7,15)
  
  #dev.off()
