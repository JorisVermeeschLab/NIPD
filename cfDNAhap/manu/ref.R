.libPaths(c("/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/R/3.4.2",.libPaths()))

list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

args <- commandArgs(TRUE)
mathap <- args[1] #heteromat.aerc.CBSseg.resolved.tsv
scmathap <- args[2] #PGD074_C3test_Itp2_E11.hap.mat.resolve.2.block
aerc <- args[3]
samplename <- args[4]
figdir <- "/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/"


EvalRefAlleleFreqLocal <- function(mathap, aerc, samplename) {
	mathap <- read.delim(mathap, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	mathap <- mathap[,c(1,2,6,7)]

	aerc <- read.delim(aerc, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	aerc <- aerc[,c(1,2,6,7)]
	colnames(aerc) <- c("chr", "pos", "refcnt", "altcnt")
	
	dt <- merge(mathap, aerc, by=c("chr","pos"))
	
	RAF <- data.frame(chr=character(), position=character(), raf=numeric(), type=character(), resHL=character(), stringsAsFactors=FALSE)
	for (i in 1:nrow(dt)) {
		if (dt$resHL[i] == "MH1" && dt$type[i] == "M2") {
			pfline <- data.frame(dt$chr[i], dt$pos[i], dt$refcnt[i]/(dt$refcnt[i]+dt$altcnt[i]), dt$type[i], dt$resHL[i])

			RAF <- rbind(RAF, pfline)
		} else if (dt$resHL[i] == "MH2" && dt$type[i] == "M1") {
			pfline <- data.frame(dt$chr[i], dt$pos[i], dt$refcnt[i]/(dt$refcnt[i]+dt$altcnt[i]), dt$type[i], dt$resHL[i])

                        RAF <- rbind(RAF, pfline)
		}
	}
        colnames(RAF) <- c("chr", "position", "raf", "type", "resHL")
	RAF <- subset(RAF, chr != "X" & chr != "Y")
	RAF$chr <- paste0("chr", RAF$chr)
	fo1 <- paste0(figdir, samplename, ".cfdna.raf.bias.tsv")
        write.table(RAF, fo1, row.names=FALSE, quote=FALSE, sep="\t")
	return(RAF)
}

EvalRefAlleleFreqRef <- function(mathap, scmathap, aerc, samplename) {
    mathap <- read.delim(mathap, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    mathap <- mathap[,c(1,2,6,7)]

    aerc <- read.delim(aerc, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    aerc <- aerc[,c(1,2,6,7)]
    colnames(aerc) <- c("chr", "pos", "refcnt", "altcnt")

    dt <- merge(mathap, aerc, by=c("chr","pos"))
    dt <- dt[, c("chr", "pos", "refcnt", "altcnt", "type")]
	
    scmathap <- read.delim(scmathap, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    scmathap$resHL <- as.character(scmathap$resHL)

    procdt <- data.frame(chr=character(), pos=character(), refcnt=numeric(), altcnt=numeric(), type=character(), stringsAsFactors=FALSE)
    tempdt <- data.frame(chr=character(), pos=character(), refcnt=numeric(), altcnt=numeric(), type=character(), stringsAsFactors=FALSE)
    
    for (i in 1:nrow(scmathap)) {
        if (scmathap$resHL[i]=="1") {
	    tempdt <- dt[(dt$chr==scmathap$chr[i] & dt$pos>=scmathap$start[i] & dt$pos<=scmathap$end[i] & dt$type == "M2"),]
	    #print(tempdt)
	} else if (scmathap$resHL[i]=="2") {
	    tempdt <- dt[(dt$chr==scmathap$chr[i] & dt$pos>=scmathap$start[i] & dt$pos<=scmathap$end[i] & dt$type == "M1"),]
	}
	
	procdt <- rbind(procdt, tempdt)    
    }
    raf <- procdt$refcnt/(procdt$refcnt + procdt$altcnt)
    fdt <- cbind(procdt, raf)
    fdt <- fdt[,c("chr", "pos", "raf", "type")]
    fdt <- subset(fdt, chr != "X" & chr != "Y")
    fdt$chr <- paste0("chr", fdt$chr)
    fo1 <- paste0(figdir, samplename, ".cfdna.raf.bias.tsv")
    write.table(fdt, fo1, row.names=FALSE, quote=FALSE, sep="\t")
    return(fdt)
}

#RAF <- EvalRefAlleleFreqLocal(mathap, aerc, samplename)
RAF <- EvalRefAlleleFreqRef(mathap, scmathap, aerc, samplename)
chrorder <- paste0("chr", as.character(seq(1,22)))
RAF$chr <- factor(RAF$chr, levels = chrorder)
fo2 <- paste0(figdir, samplename, ".cfdna.raf.plot.chr.pdf")
pdf(fo2, width=8, height=8)
ggplot(RAF, aes(x=raf, color=type)) + geom_density(aes(y=..scaled..)) + facet_wrap(~chr) + theme_light() + labs(x="Reference Allele Frequency(cfDNA)", y="Density") + theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

fo3 <- paste0(figdir, samplename, ".cfdna.raf.plot.pdf")
pdf(fo3, width=6, height=6)
ggplot(RAF, aes(x=raf, color=type)) + geom_density(aes(y=..scaled..)) + theme_light() + labs(x="Reference Allele Frequency(cfDNA)", y="Density") + theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#######combined

s1 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family1_181.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
s2 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family2_186.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
s3 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family3_085.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
s4 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family8_078.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
s5 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family9_074.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
s6 <- read.delim("/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/supp/Family10_001.cfdna.raf.bias.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
Sample <- rep("Family1_181", nrow(s1))
s1 <- cbind(s1, Sample)
Sample <- rep("Family2_186", nrow(s2))
s2 <- cbind(s2, Sample)
Sample <- rep("Family3_085", nrow(s3))
s3 <- cbind(s3, Sample)
Sample <- rep("Family8_078", nrow(s4))
s4 <- cbind(s4, Sample)
Sample <- rep("Family9_074", nrow(s5))
s5 <- cbind(s5, Sample)
Sample <- rep("Family10_001", nrow(s6))
s6 <- cbind(s6, Sample)

dt <- rbind(s1,s2,s3,s4,s5,s6)
fo4 <- paste0(figdir, "cfdna.raf.crosample.pdf")
pdf(fo4, width=8, height=6)
ggplot(dt, aes(x=raf, color=Sample)) + geom_density(aes(y=..scaled..)) + theme_light() + labs(x="Reference Allele Frequency(cfDNA)", y="Density") + theme(legend.position=c(0.85,0.75), axis.text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1))+scale_x_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
dev.off()
