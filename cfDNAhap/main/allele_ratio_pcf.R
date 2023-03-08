####using piecewise constant fitting to do infer haplotype segmentation
####pcf script adapted from Nilsen et al. BMC Genomics 2012, 13:591 Copynumber: Efficient algorithms for single and multi-track copy number segmentation

args = commandArgs(TRUE)
snptype <- args[1]
pref <- args[2] 
ratio <- args[3]
winsize <- args[4]
kmin <- args[5]
gamma1 <- args[6]
gamma2 <- args[7]
penalty <- args[8]
cffhapdir <- args[9]
gender <- args[10]

script.dir <- dirname(sys.frame(1)$ofile)
fastPCFR <- file.path(script.dir, "fastPCF.R")
source(fastPCFR)

Run_PCF = function(x, winsize, kmin, gamma, penalty) {
    sdev = getMad(x, winsize)
    res = selectFastPcf(x, kmin, gamma*sdev, penalty)
    yhat = res$yhat
    return(yhat)
}


Paternal_FetalAF_Fitting = function(patfh, ratio, winsize, kmin, gamma1, penalty) {
    winsize = as.numeric(winsize)
    kmin = as.numeric(kmin)
    gamma1 = as.numeric(gamma1)
    fh = read.delim(patfh, header=TRUE, sep="\t")
    fhP1 = subset(fh, TYPE=="P1")
    fhP2 = subset(fh, TYPE=="P2")
    chrlist = as.character(unique(fh$X.CHROM))
    tdf = data.frame(chr=character(), pos=numeric(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    for (i in 1:length(chrlist)) {
	fhP1chr = subset(fhP1, X.CHROM == chrlist[i])
	fhP2chr = subset(fhP2, X.CHROM == chrlist[i])
        #convert ratio in P1 category to negative values for visualization
	fhP1chr[, ratio] = 0 - fhP1chr[, ratio]
        xP1 = fhP1chr[, ratio]
	xP2 = fhP2chr[, ratio]
	if (length(xP1) != 0 && length(xP2) != 0) {
            yhatP1 = Run_PCF(xP1, winsize, kmin, gamma1, penalty)
            yhatP2 = Run_PCF(xP2, winsize, kmin, gamma1, penalty)
            calP1 = cbind(fhP1chr[,c(1,2)], xP1, yhatP1, rep("P1", length(yhatP1)))
            calP2 = cbind(fhP2chr[,c(1,2)], xP2, yhatP2, rep("P2", length(yhatP2)))
	    colnames(calP1) = colnames(calP2) = c("chr","pos",ratio,"yhat","type")
	    calchr = rbind(calP1, calP2)
	    #Y chromosome are all P2 category by definition
	} else if (length(xP1) == 0 && length(xP2) != 0) {
	    yhatP2 = Run_PCF(xP2, winsize, kmin, gamma1, penalty)
	    calP2 = cbind(fhP2chr[,c(1,2)], xP2, yhatP2, rep("P2", length(yhatP2)))
            colnames(calP2) = c("chr","pos",ratio,"yhat","type")
	    calchr = calP2
	}
	tdf = rbind(tdf, calchr)
    }
    fo = gsub(".tsv", ".bp.tsv", patfh)
    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
}


Maternal_FetalAF_Fitting = function(matfh, ratio, winsize, kmin, gamma2, penalty, gender) {
    winsize = as.numeric(winsize)
    kmin = as.numeric(kmin)
    gamma2 = as.numeric(gamma2)
    fh = read.delim(matfh, header=TRUE, sep="\t")
    fhM1 = subset(fh, TYPE=="M1")
    fhM2 = subset(fh, TYPE=="M2")
    chrlist = as.character(unique(fh$X.CHROM))
    tdf = data.frame(chr=character(), pos=numeric(), RAF=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    for (i in 1:length(chrlist)) {
        if (gender == "male" & chrlist[i] == "X") {
	    fhchr = subset(fh, X.CHROM == chrlist[i])
            x = 2*fhchr[, ratio]
            yhat = Run_PCF(x, winsize, kmin, gamma1, penalty)
            calchr = cbind(fhchr[,c(1,2)], x, yhat, rep("M1", length(yhat)))
            colnames(calchr) = c("chr","pos",ratio,"yhat","type")
	} else {
	    fhM1chr = subset(fhM1, X.CHROM == chrlist[i])
	    fhM2chr = subset(fhM2, X.CHROM == chrlist[i])
	    xM1 = fhM1chr[, ratio]
	    xM2 = fhM2chr[, ratio]
	    yhatM1 = Run_PCF(xM1, winsize, kmin, gamma2, penalty)
	    yhatM2 = Run_PCF(xM2, winsize, kmin, gamma2, penalty)
	    calM1 = cbind(fhM1chr[,c(1,2)], xM1, yhatM1, rep("M1", length(yhatM1)))
	    calM2 = cbind(fhM2chr[,c(1,2)], xM2, yhatM2, rep("M2", length(yhatM2)))
	    colnames(calM1) = colnames(calM2) = c("chr","pos",ratio,"yhat","type")
	    calchr = rbind(calM1, calM2)
        }
	tdf = rbind(tdf, calchr)
    }
    fo = gsub(".tsv", ".bp.tsv", matfh)
    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
}



pcf_main = function(snptype, pref, ratio, winsize, kmin, gamma1, gamma2, penalty, cffhapdir, gender) {
    if (snptype == 'heteropat') {
        infh = paste0(cffhapdir, pref, '.heteropat.aerc.tsv')
        Paternal_FetalAF_Fitting(infh, ratio, winsize, kmin, gamma1, penalty)
    } else if (snptype == 'heteromat') {
        infh = paste0(cffhapdir, pref, '.heteromat.aerc.tsv')
        Maternal_FetalAF_Fitting(infh, ratio, winsize, kmin, gamma2, penalty, gender)
    } else if (snptype == 'heteromatboth') {
        infh = paste0(cffhapdir, pref, '.heteromatboth.aerc.tsv')
        Maternal_FetalAF_Fitting(infh, ratio, winsize, kmin, gamma2, penalty, gender)
    }
}

pcf_main(snptype, pref, ratio, winsize, kmin, gamma1, gamma2, penalty, cffhapdir, gender)



#specify the ratio to use - MH1_ratio or overhap_ratio
#Maternal_FetalAF_Fitting = function(matfh, ratio, winsize, kmin, gamma1, gamma2, penalty, workdir) {
#    fh = read.delim(matfh, header=TRUE, sep="\t")
#    fhM1 = subset(fh, Type=="M1")
#    fhM2 = subset(fh, Type=="M2")
#    chrlist = c(as.character(seq(1, 22, 1)), "X")
#    tdf = data.frame(chr=character(), pos=numeric(), RAF=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
#    for (i in 1:length(chrlist)) {
#        fhM1chr = subset(fhM1, X.chrom == chrlist[i])
#        fhM2chr = subset(fhM2, X.chrom == chrlist[i])
#        xM1 = fhM1chr[, ratio]
#        xM2 = fhM2chr[, ratio]
#        yhatM1 = Run_PCF(xM1, winsize, kmin, gamma1, penalty)
#        yhatM2 = Run_PCF(xM2, winsize, kmin, gamma2, penalty)
#        calM1 = cbind(fhM1chr[,c(1,2)], xM1, yhatM1, rep("M1", length(yhatM1)))
#        calM2 = cbind(fhM2chr[,c(1,2)], xM2, yhatM2, rep("M2", length(yhatM2)))
#        colnames(calM1) = colnames(calM2) = c("chr","pos",ratio,"yhat","type")
#        calchr = rbind(calM1, calM2)
#        tdf = rbind(tdf, calchr)
#    }
#    fo = paste0(workdir, ratio, ".tsv")
#    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
#}
#
#Maternal_FetalAF_Fitting2 = function(matfh, ratio, winsize, kmin, gamma1, gamma2, penalty, workdir) {
#    fh = read.delim(matfh, header=TRUE, sep="\t")
#    fhM1 = subset(fh, Type=="M1")
#    fhM2 = subset(fh, Type=="M2")
#    chrlist = c(as.character(seq(1, 22, 1)), "X")
#    tdf = data.frame(chr=character(), pos=numeric(), RAF=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
#    for (i in 1:length(chrlist)) {
#        if (chrlist[i] != "X") {
#            fhM1chr = subset(fhM1, X.chrom == chrlist[i])
#            fhM2chr = subset(fhM2, X.chrom == chrlist[i])
#            xM1 = fhM1chr[, ratio]
#            xM2 = fhM2chr[, ratio]
#            yhatM1 = Run_PCF(xM1, winsize, kmin, gamma1, penalty)
#            yhatM2 = Run_PCF(xM2, winsize, kmin, gamma2, penalty)
#            calM1 = cbind(fhM1chr[,c(1,2)], xM1, yhatM1, rep("M1", length(yhatM1)))
#            calM2 = cbind(fhM2chr[,c(1,2)], xM2, yhatM2, rep("M2", length(yhatM2)))
#            colnames(calM1) = colnames(calM2) = c("chr","pos",ratio,"yhat","type")
#            calchr = rbind(calM1, calM2)
#            tdf = rbind(tdf, calchr)
#        } else if (chrlist[i] == "X") {
#            fhchr = subset(fh, X.chrom == chrlist[i])
#            x = fhchr[, ratio]
#            yhat = Run_PCF(x, winsize, kmin, gamma1, penalty)
#            calchr = cbind(fhchr[,c(1,2)], x, yhat, rep("M1", length(yhat)))
#            colnames(calchr) = c("chr","pos",ratio,"yhat","type")
#            tdf = rbind(tdf, calchr)
#        }
#    }
#    fo = paste0(workdir, ratio, ".v2.tsv")
#    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
#}
#
#Paternal_FetalAF_Fitting = function(patfh, ratio, winsize, kmin, gamma3, penalty, workdir) {
#    fh = read.delim(patfh, header=TRUE, sep="\t")
#    fhP1 = subset(fh, Type=="P1")
#    fhP2 = subset(fh, Type=="P2")
#    chrlist = c(as.character(seq(1, 22, 1)), "Y")
#    tdf = data.frame(chr=character(), pos=numeric(), RAF=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
#    for (i in 1:length(chrlist)) {
#	fhP1chr = subset(fhP1, X.chrom == chrlist[i])
#	fhP2chr = subset(fhP2, X.chrom == chrlist[i])
#        #convert ratio in P2 category to negative values for visualization
#	fhP2chr[, ratio] = 0 - fhP2chr[, ratio]
#        xP1 = fhP1chr[, ratio]
#	xP2 = fhP2chr[, ratio]
#	if (length(xP1) != 0 && length(xP2) != 0) {
#	    yhatP1 = Run_PCF(xP1, winsize, kmin, gamma3, penalty)
#	    yhatP2 = Run_PCF(xP2, winsize, kmin, gamma3, penalty)
#	    calP1 = cbind(fhP1chr[,c(1,2)], xP1, yhatP1, rep("P1", length(yhatP1)))
#	    calP2 = cbind(fhP2chr[,c(1,2)], xP2, yhatP2, rep("P2", length(yhatP2)))
#	    colnames(calP1) = colnames(calP2) = c("chr","pos",ratio,"yhat","type")
#	    calchr = rbind(calP1, calP2)
#	} else if (length(xP1) != 0 && length(xP2) == 0) {
#	    yhatP1 = Run_PCF(xP1, winsize, kmin, gamma3, penalty)
#	    calP1 = cbind(fhP1chr[,c(1,2)], xP1, yhatP1, rep("P1", length(yhatP1)))
#	    colnames(calP1) = c("chr","pos",ratio,"yhat","type")
#	    calchr = calP1
#	}
#	tdf = rbind(tdf, calchr)
#    }
#    fo = paste0(workdir, ratio, ".tsv")
#    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
#}
