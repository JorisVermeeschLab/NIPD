list.of.packages <- c("ggplot2", "gridExtra", "grid", "lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

args = commandArgs(TRUE)
heteropat = args[1]
heteromat = args[2]
heteroboth = args[3]
heteropataerc = args[4]
heteromataerc = args[5]
heterobothaerc = args[6]
pref = args[7]
workdir = args[8]

plot_heteropat_raf <- function(heteropat, heteropataerc) {
    #read in snp type aerc file, which is a subset of the snp type file; thus only to plot the distribution of parental heterozygous sites that passed filtering in cfDNA 
    heteropataerc = read.delim(heteropataerc, comment.char = "#", sep = "\t", header = FALSE)
    heteropataercdf = heteropataerc[, c(1,2)]
    colnames(heteropataercdf) = c("CHROM", "POS")
    heteropat = read.delim(heteropat, comment.char = "#", sep = "\t", header = FALSE)
    heteropatdf = heteropat[, c(1,2,7,10,16)]
    colnames(heteropatdf) = c("CHROM", "POS", "PATRAF", "PH1", "TYPE")
    heteropatdf = merge(heteropataercdf, heteropatdf, by=c("CHROM", "POS"))
    heteropatdf$PATRAF = as.numeric(as.character(heteropatdf$PATRAF))
    pat = ggplot(heteropatdf, aes(x=PATRAF)) + geom_histogram(color = "darkblue", na.rm = TRUE, binwidth = 0.01) + facet_grid(TYPE~PH1) + theme_minimal() + labs(x="RAF(heteropat)") + theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
    return(pat)
}

plot_heteromat_raf <- function(heteromat, heteromataerc) {
    heteromataerc = read.delim(heteromataerc, comment.char = "#", sep = "\t", header = FALSE)
    heteromataercdf = heteromataerc[, c(1,2)]
    colnames(heteromataercdf) = c("CHROM", "POS")
    heteromat = read.delim(heteromat, comment.char = "#", sep = "\t", header = FALSE)
    heteromatdf = heteromat[, c(1,2,8,12,16)]
    colnames(heteromatdf) = c("CHROM", "POS", "MATRAF", "MH1", "TYPE")
    heteromatdf = merge(heteromataercdf, heteromatdf, by=c("CHROM", "POS"))
    heteromatdf$MATRAF = as.numeric(as.character(heteromatdf$MATRAF))
    mat = ggplot(heteromatdf, aes(x=MATRAF)) + geom_histogram(color = "darkblue", na.rm = TRUE, binwidth = 0.01) + facet_grid(TYPE~MH1) + theme_minimal() + labs(x="RAF(heteromat)") + theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
    return(mat)
}

plot_heteroboth_raf <- function(heteroboth, heterobothaerc) {
    heterobothaerc = read.delim(heterobothaerc, comment.char = "#", sep = "\t", header = FALSE)
    heterobothaercdf = heterobothaerc[, c(1,2)]
    colnames(heterobothaercdf) = c("CHROM", "POS")
    heteroboth = read.delim(heteroboth, comment.char = "#", sep = "\t", header = FALSE)
    bothpatdf = cbind(heteroboth[, c(1,2,7,10)], rep("Pat", dim(heteroboth)[1]))
    bothmatdf = cbind(heteroboth[, c(1,2,8,12)], rep("Mat", dim(heteroboth)[1]))
    colnames(bothpatdf) = colnames(bothmatdf) = c("CHROM", "POS", "RAF", "H1", "TYPE")
    bothpatdf = merge(heterobothaercdf, bothpatdf, by=c("CHROM", "POS"))
    bothmatdf = merge(heterobothaercdf, bothmatdf, by=c("CHROM", "POS"))
    hetbothdf = rbind(bothpatdf, bothmatdf)
    hetbothdf$RAF = as.numeric(as.character(hetbothdf$RAF))
    hetboth = ggplot(hetbothdf, aes(x=RAF)) + geom_histogram(color = "darkblue", na.rm = TRUE, binwidth = 0.01) + facet_grid(TYPE~H1) + theme_minimal() + labs(x="RAF(heteroboth)") + theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
    return(hetboth)
}

plot_raf_main <- function(heteropat, heteromat, heteroboth, heteropataerc, heteromataerc, heterobothaerc, pref, workdir) {
    figdir = paste0(workdir, "6_fig/")
    if (dir.exists(figdir) != TRUE) {
	dir.create(figdir)
    }
    fo = paste0(figdir, pref, ".raf.pdf")
    pdf(file = fo, width = 7, height = 4)
    pat = plot_heteropat_raf(heteropat, heteropataerc)
    mat = plot_heteromat_raf(heteromat, heteromataerc)
    hetboth = plot_heteroboth_raf(heteroboth, heterobothaerc)
    grid.arrange(pat, mat, hetboth, ncol=3)
    dev.off()
}

plot_raf_main(heteropat, heteromat, heteroboth, heteropataerc, heteromataerc, heterobothaerc, pref, workdir)
