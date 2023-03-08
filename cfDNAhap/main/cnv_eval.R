library(tidyr)
script.dir <- dirname(sys.frame(1)$ofile)
fastPCFR <- file.path(script.dir, "fastPCF.R")
source(fastPCFR)

GetTarGC <- function(targetgcfh) {
    gcfh = read.delim(targetgcfh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    gcfh = gcfh[,c(1,2,3,6,13,14)]
    colnames(gcfh) = c("CHROM","START","END","GC","TARSIZE","INDEX")
    return(gcfh)
}

Run_PCF <- function(x, winsize, kmin, gamma, penalty) {
    sdev = getMad(x, winsize)
    res = selectFastPcf(x, kmin, gamma*sdev, penalty)
    yhat = res$yhat
    return(yhat)
}

#Cal_Frag_CopyRatio <- function(fragcov.sis, sample.name, targetgcfh, winsize, kmin, gamma, penalty, cnvdir, metric="totalcov", sex) {
#    #fragcov.sis - the output file from GATK coverage analysis with suffix (fragcov).sample_interval_summary
#    fragcov = read.delim(fragcov.sis, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#    #get Target, total_coverage, granular_median
#    fragcov = fragcov[,c(1,2,7)]
#    fragnew = fragcov %>% separate(col=Target, into=c("chrom","start"), sep=":") %>% separate(col=start, into=c("start", "end"), sep="-")
#    colnames(fragnew) = c("chrom", "start", "end", "totalcov", "median")
#    fragnew$median <- gsub(">500", 500, fragnew$median)
#    fragnew$totalcov <- gsub(">500", 500, fragnew$totalcov)
#    fragnew$totalcov <- as.numeric(fragnew$totalcov)
#    fragnew$start = as.numeric(fragnew$start)-1
#    fragnew$chrom = as.character(fragnew$chrom)
    
#    #output raw normalized target coverage for sample comparison
#    fragnew$median <- as.numeric(fragnew$median)
#    autofragnew <- subset(fragnew, chrom != "X" | chrom != "Y")
#    tarintervalmed <- median(autofragnew$median, na.rm=TRUE)
#    fragraw <- fragnew[,c("chrom", "start", "end", "median")]
#    fragraw$median <- fragraw$median/tarintervalmed
	
#    fo2 = paste0(cnvdir, sample.name, ".", metric ,".median.norm.rawcov.bed")
#    write.table(fragraw, file = fo2, sep = "\t", row.names = FALSE, quote = FALSE)


#    if (sex == "male") {
#        fragnew$totalcov[fragnew$chrom=="X"] <- 2*as.numeric(fragnew$totalcov[fragnew$chrom=="X"])
#        fragnew$median[fragnew$chrom=="X"]  <- 2*as.numeric(fragnew$median[fragnew$chrom=="X"])
#    }

#    gcfh = Read_TargetGCFile(targetgcfh)
#    fragmerge = merge(fragnew, gcfh, by=c("chrom", "start", "end"))
#    #exclude Y chromosome
#    fragmerge = subset(fragmerge, fragmerge$chrom != "Y")
#    if (metric=="totalcov") {
	#normalize total coverage per target by target size, totalcov of the target is used by default, seems to work better than median cov of the target
#	fragmerge$tcov.size = fragmerge$totalcov/fragmerge$tgsize
#    } else if (metric=="median") {
#	fragmerge$tcov.size = as.numeric(fragmerge$median)
#    } else {
#	print("The metric is unspecified. Please specify metric as totalcov or median")
#    }
    ##lowerquant = quantile(fragmerge$tcov.size, 0.05)
    ##higherquant = quantile(fragmerge$tcov.size, 0.985)
    #remove outliers - bottom and top 1% quantiles
    ##fragmerge = subset(fragmerge, fragmerge$tcov.size >= lowerquant & fragmerge$tcov.size <= higherquant)
    #LOESS is unbounded, when the coverage is close to 1, the predicted value can go to negative; work on log scale of the coverage, to avoid zero coverage, add 0.1 to coverage
    #fragmerge$tcov.size = log2(fragmerge$tcov.size+0.1)
    #loess GC correction
    #break up by GC level, and take median in each GC level group
#    fragloess = subset(fragmerge, fragmerge$chrom != "X")
#    targcmedian = tapply(fragloess$tcov.size, fragloess$gc, median)
#    cor.factor.gc = predict(loess(targcmedian ~ as.numeric(names(targcmedian))), fragmerge$gc)
#    fragmerge$tcov.gc = fragmerge$tcov.size * median(fragmerge$tcov.size)/cor.factor.gc
#    fragmerge$gccor.cov = fragmerge$tcov.gc / median(fragmerge$tcov.gc, na.rm=TRUE)
#    #apply log2 on positive values
#    fragmerge = subset(fragmerge, fragmerge$gccor.cov > 0)
#    fragmerge$logcov = log2(fragmerge$gccor.cov)
#    chrlist = as.character(unique(fragmerge$chrom))
#    tdf = data.frame(chrom=character(), start=numeric(), end=numeric(), logcov=numeric(), logcov.pcf=character(), stringsAsFactors=FALSE)
#    for (i in 1:length(chrlist)) {
#        fhchr = subset(fragmerge, chrom == chrlist[i])
#	lowerquant = quantile(fhchr$logcov, 0.015)
#	higherquant = quantile(fhchr$logcov, 0.985)
#	fhchr = subset(fhchr, fhchr$logcov>=lowerquant & fhchr$logcov<=higherquant)
#	fhchr$logcov.pcf = Run_PCF(fhchr$logcov, winsize, kmin, gamma, penalty)
#	calfh = fhchr[,c("chrom", "start", "end", "logcov", "logcov.pcf")]
#	tdf = rbind(tdf, calfh) 
#    }
#    fo = paste0(cnvdir, sample.name, ".", metric ,".gccor.noref.logR.tsv")
#    write.table(tdf, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)	
#}

GetSampleTarCount <- function(fragcov.sis, samplename, targetgcfh, cnvdir) {

    #samplefh <- paste0(sampledir, "/4_cov/", samplename, "/", samplename, ".sorted.merged.md.filtered.fragcov.sample_interval_summary")
    gcfh <- GetTarGC(targetgcfh)
    fragcov <- read.delim(fragcov.sis, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    #get Target, avg_coverage (= total.fragment.coverage/target.size)
    fragcov <- fragcov[,c(1,2,3)]
    fragnew <- fragcov %>% separate(col=Target, into=c("chrom","start"), sep=":") %>% separate(col=start, into=c("start", "end"), sep="-")
    colnames(fragnew) <- c("CHROM", "START", "END", "TOTAL", "MEAN")
    fragnew$START <- as.numeric(fragnew$START)-1
    fragnew$CHROM <- as.character(fragnew$CHROM)

    #generate chromosome-wise stats
    chrstats <- data.frame(CHROM=character(), NORM.COV=numeric())
    chrdt <- merge(fragnew, gcfh, by=c("CHROM", "START", "END"))
    allmedian <- median(chrdt$MEAN)
    chrlist <- unique(chrdt$CHROM)
    for (i in 1:length(chrlist)) {
        fhchrdt <- subset(chrdt, CHROM==chrlist[i])
	tarsize <- sum(fhchrdt$TARSIZE)
	tarcov <- sum(fhchrdt$TOTAL)
	chravgcov <- tarcov/tarsize
	norm.chr <- chravgcov/allmedian
	chrstats <- rbind(chrstats, data.frame(chrlist[i], norm.chr))

    }
    colnames(chrstats) <- c("CHROM", "NORM.COV")
    fochr <- paste0(cnvdir, "/", samplename, ".chr.stats.tsv")
    write.table(chrstats, fochr, sep="\t", row.names=FALSE, quote=FALSE)

    #apply GC correction

    fragmerge <- merge(fragnew, gcfh, by=c("CHROM", "START", "END"))
   
    fragloess <- subset(fragmerge, fragmerge$CHROM != "X" & fragmerge$CHROM != "Y" & fragmerge$MEAN != 0)
    #median coverage per GC level
    targcmedian <- tapply(fragloess$MEAN, fragloess$GC, median)

    gc.cor.factor <- predict(loess(targcmedian ~ as.numeric(names(targcmedian))), fragmerge$GC)
    fragmerge$GCC.MEAN <- fragmerge$MEAN * median(fragmerge$MEAN)/gc.cor.factor
    fragmerge$GCC.MEAN[!is.finite(fragmerge$GCC.MEAN)] <- 0

    #Normalization
    autofragmerge <- subset(fragmerge, CHROM != "X" & CHROM != "Y")
    auto.raw.median <- median(as.numeric(autofragmerge$MEAN), na.rm=TRUE)
    auto.gcc.median <- median(as.numeric(autofragmerge$GCC.MEAN), na.rm=TRUE)
    fragmerge$MEAN.RAW.NORM <- fragmerge$MEAN/auto.raw.median
    fragmerge$MEAN.GCC.NORM <- fragmerge$GCC.MEAN/auto.gcc.median
    fragmerge$MEAN.GCC.NORM[fragmerge$MEAN.GCC.NORM < 0] <- 0

    fragout <- fragmerge[,c("INDEX", "CHROM", "START", "END", "MEAN", "MEAN.RAW.NORM", "MEAN.GCC.NORM")]
    #write.table(fragout, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)
    fragout <- fragout[order(fragout$INDEX),]
    return(fragout)
}


Cal_Frag_CopyRatio <- function(fragcov.sis, samplename, targetgcfh, winsize, kmin, gamma, penalty, cnvdir) {
    fo <- paste0(cnvdir, "/", samplename, ".target.coverage.noref.tsv")
    
    fragout <- GetSampleTarCount(fragcov.sis, samplename, targetgcfh, cnvdir)

    chrlist <- as.character(unique(fragout$CHROM))
    tdf <- data.frame(INDEX=numeric(), CHROM=character(), START=numeric(), END=numeric(), MEAN=numeric(), MEAN.RAW.NORM=numeric(), MEAN.GCC.NORM=numeric(), MEAN.GCC.NORM.PCF=numeric(), stringsAsFactors=FALSE)
    for (i in 1:length(chrlist)) {
        fhchr <- subset(fragout, CHROM == chrlist[i])
        fhchr$MEAN.GCC.NORM.PCF = Run_PCF(fhchr$MEAN.GCC.NORM, winsize, kmin, gamma, penalty)
        tdf <- rbind(tdf, fhchr)
    }
    tdf <- tdf[order(tdf$INDEX),]
    write.table(tdf, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)

}



Cal_HomoDiffest_Seg <- function(HomoDiffest, winsize, kmin, gamma1, penalty, cnvdir) {
    homodiffest <- read.delim(HomoDiffest, header = TRUE, sep = '\t', stringsAsFactors=FALSE)
    chrlist <- as.character(unique(homodiffest$X.CHROM))
    tdf <- data.frame(CHROM=character(), POS=numeric(), PH1=character(), FFest=numeric(), FFest.pcf=numeric())
    for (i in 1:length(chrlist)) {
	fhchr <- subset(homodiffest, X.CHROM == chrlist[i])
	fhchr$FFest.pcf <- Run_PCF(fhchr$FFest, winsize, kmin, gamma1, penalty)
	tdf <- rbind(tdf, fhchr)
    }
    fo <- paste0(cnvdir, "homodiff.ffest.tsv")
    write.table(tdf, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)
}


CNV_eval_main <- function(fragcov.sis.list, targetgcfh, HomoDiffest, winsize, kmin, gamma, gamma1, penalty, workdir) {
    winsize <- as.numeric(winsize)
    kmin <- as.numeric(kmin)
    gamma <- as.numeric(gamma)
    gamma1 <- as.numeric(gamma1)
    cnvdir <- paste0(workdir, "5_cnv/")
    if (dir.exists(cnvdir) != TRUE) {
	dir.create(cnvdir)
    }

    fhlist <- read.delim(fragcov.sis.list, header = TRUE)
        for (i in 1:dim(fhlist)[1]) {
	    fragcov.sis <- as.character(fhlist[i,1])
	    samplename <- as.character(fhlist[i,2])
            sex <- as.character(fhlist[i,3])
	    Cal_Frag_CopyRatio(fragcov.sis, samplename, targetgcfh, winsize, kmin, gamma, penalty, cnvdir)
	}

    Cal_HomoDiffest_Seg(HomoDiffest, winsize, kmin, gamma1, penalty, cnvdir)
}
