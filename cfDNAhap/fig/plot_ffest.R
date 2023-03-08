list.of.packages <- c("ggplot2", "gridExtra", "grid", "lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

args = commandArgs(TRUE)
homodiffest = args[1]
heteropatest = args[2]
pref = args[3]
workdir = args[4]
ffsex = args[5]
ffsex = as.character(ffsex)

readin_ffestfh <- function(homodiffest, heteropatest, ffsex) {
    homodiffest = read.delim(homodiffest, header = TRUE, sep = "\t")
    heteropatest = read.delim(heteropatest, header = TRUE, sep = "\t")
    #ffdf = data.frame(CHR=character(), POS=numeric(), FFest=numeric(), Estimation=character(), stringsAsFactors=FALSE)
    homodiffdf = subset(homodiffest, X.CHROM != "X" & X.CHROM != "Y")
    homodiffdf = cbind(homodiffdf[,c(1,2,4)], rep("Diff.Homo.SNPs", dim(homodiffdf)[1]))
    heteropatdf = subset(heteropatest, X.CHROM != "X" & X.CHROM != "Y" & as.numeric(FFest) > 0.01)
    heteropatdf = cbind(heteropatdf[,c(1,2,4)], rep("HeteroPat.SNPs", dim(heteropatdf)[1]))
    if (ffsex == "male") {
        ychrdf = subset(heteropatest, X.CHROM == "Y")
        ychrdf = cbind(ychrdf[,c(1,2,4)], rep("ChrY.SNPs", dim(ychrdf)[1]))
        colnames(homodiffdf) = colnames(heteropatdf) = colnames(ychrdf) = c("CHR", "POS", "FFest", "Estimation")
        ffdf = rbind(homodiffdf, heteropatdf, ychrdf)
    } else {
        colnames(homodiffdf) = colnames(heteropatdf) = c("CHR", "POS", "FFest", "Estimation")
        ffdf = rbind(homodiffdf, heteropatdf)
    }
    return(ffdf)
}


get_ffest_median <- function(homodiffest, heteropatest, ffsex) {
    ffestfh = readin_ffestfh(homodiffest, heteropatest, ffsex)
    ffest.med = aggregate(FFest ~ Estimation, data = ffestfh, median)
    ffest.med$FFest = round(ffest.med$FFest, digits=4)
    return(ffest.med)
}


figdir = paste0(workdir, "6_fig/")
if (dir.exists(figdir) != TRUE) {
    dir.create(figdir)
}
ffestfh = readin_ffestfh(homodiffest, heteropatest, ffsex)
ffest.med = get_ffest_median(homodiffest, heteropatest, ffsex)
fo = paste0(figdir, pref, ".ffest.png")
png(file = fo, width = 560, height = 320, res = 100)
if (ffsex == "male") {
    ggplot(ffestfh, aes(x = Estimation, y = FFest, color = Estimation)) + geom_boxplot() + theme_minimal() + labs(x="SNPs type", y="Fetal Fraction Estimation") + scale_color_manual(values = c("#C60000", "#000B29", "orange")) + scale_y_continuous(limits=c(0, 0.4), breaks=seq(0, 0.4, 0.05)) + geom_text(data = ffest.med, aes(label = FFest, y = FFest + 0.015), size = 4)
} else {
    ggplot(ffestfh, aes(x = Estimation, y = FFest, color = Estimation)) + geom_boxplot() + theme_minimal() + labs(x="SNPs type", y="Fetal Fraction Estimation") + scale_color_manual(values = c("#C60000", "#000B29")) + scale_y_continuous(limits=c(0, 0.4), breaks=seq(0, 0.4, 0.05)) + geom_text(data = ffest.med, aes(label = FFest, y = FFest + 0.015), size = 4)
}
dev.off()

