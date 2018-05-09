setwd("/Users/atma/Desktop/LJ_proteomics")

data = read.table("USP9X_vs_IgG.txt", sep = "\t", header = T)
data
dimnames(data)[[1]] <- data[,1]
data = data[,-1]
data


colnames(data)

ttestData <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y,paired=TRUE)
  results$p.value
}

rawpvalue = apply(data, 1, ttestData, grp1 = c(1:3), grp2 = c(4:6))
rawpvalue

#sort p-values from lowest to highest
rawpvalue_sorted = sort(rawpvalue)
rawpvalue_sorted[1:10]

# R provides 7 options for correcting p-values: 
# fdr, bonferroni, holm, hochberg, hommel, BH, or BY
# test all of these and compare to raw p-values

# use the list of raw p-values to get a corresponding list of FDR p-values
fdr_corrected_pval = p.adjust(rawpvalue, method="fdr")
bonferroni_corrected_pval = p.adjust(rawpvalue, method="bonferroni")
holm_corrected_pval = p.adjust(rawpvalue, method="holm")
hochberg_corrected_pval = p.adjust(rawpvalue, method="hochberg")
hommel_corrected_pval = p.adjust(rawpvalue, method="hommel")
BH_corrected_pval = p.adjust(rawpvalue, method="BH")
BY_corrected_pval = p.adjust(rawpvalue, method="BY")


#fdr_corrected_pval_sorted = fdr_corrected_pval[order(fdr_corrected_pval)]
#fdr_corrected_pval_sorted
#fdr_corrected_pval_sorted[1:10]





results = cbind(rawpvalue, fdr_corrected_pval, bonferroni_corrected_pval, holm_corrected_pval, hochberg_corrected_pval, hommel_corrected_pval, BH_corrected_pval, BY_corrected_pval)
results = as.data.frame(results)
results

results_sorted = results[order(results$rawpvalue), ]
results_sorted



hist(rawpvalue)

##transform our data into log2 base.
data = log2(data)

#calculate the mean of each gene per control group
control = apply(data[,1:3], 1, mean)

#calcuate the mean of each gene per test group
test = apply(data[, 4:6], 1, mean) 

#confirming that we have a vector of numbers
class(control) 

#confirming we have a vector of numbers
class(test)

log2foldchange <- control - test 
log2foldchange

hist(log2foldchange, xlab = "log2 Fold Change (Control vs Test)")

results = cbind(log2foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)
results

#library(ggplot2)
#volcano = ggplot(data = results, aes(x = log2foldchange, y = -1*log10(rawpvalue)))
#volcano + geom_point()

# draw a basic volcano plot

pdf('test.pdf')

#with(results,plot(log2foldchange, -log10(rawpvalue),pch=20,main="Volcano plot"))
with(results,plot(log2foldchange, -log10(rawpvalue),pch=20,main="USP9X_vs_IgG paired t-test",xlab="log2(foldchange)",xlim=c(-2.75,2)))

# add colored points: red if pvalue<0.05, orange if log2FC>1, green if both
with(subset(results, rawpvalue<0.1),points(log2foldchange,-log10(rawpvalue),pch=20,col="red"))

with(subset(results,abs(log2foldchange)>1),points(log2foldchange,-log10(rawpvalue),pch=20,col="orange"))

test = subset(results,abs(log2foldchange)>1)
test

with(subset(results,rawpvalue<0.1 & abs(log2foldchange)>1),points(log2foldchange,-log10(rawpvalue),pch=20,col="green"))

# label points with the textxy function from the calibrate package
#install.packages("calibrate")
library(calibrate)

# this will label just the ones that are significant with high fold change
with(subset(results,rawpvalue<0.1 & abs(log2foldchange)>1),textxy(log2foldchange,-log10(rawpvalue),labs=probename,cex=.5))

# this will also label ones with high fold change (no restriction on pvalue)
with(subset(results,abs(log2foldchange)>1),textxy(log2foldchange,-log10(rawpvalue),labs=probename,cex=.5))

# this will also label ones which are significant (no restriction on fold change)
with(subset(results,rawpvalue<0.1 & abs(log2foldchange)>0),textxy(log2foldchange,-log10(rawpvalue),labs=probename,cex=.5))

dev.off()
