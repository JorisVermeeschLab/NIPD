.libPaths(c("/scratch/leuven/323/vsc32370/sw/R/3.4.2",.libPaths()))

list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args = commandArgs(TRUE)
patfh = args[1]
matfh = args[2]
logcovfh = args[3]
homodiffest = args[4]
patlogcovfh = args[5]
matlogcovfh = args[6]
samplename = args[7]
figdir = args[8]
ffest_scale = args[9]
scriptdir = args[10]
ffsex = args[11]

Get_hg19_chrsize <- function(scriptdir) {
    fh <- paste0(scriptdir, "/",  "hg19.chrom.size.bed")
    fh <- read.delim(fh, header=FALSE, sep="\t")
    colnames(fh) <- c("chromosome", "start", "end")
    return(fh)
}

Get_SNPratio <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 <- fh$pos-1
        fh <- cbind(fh[,1], pos1, fh[,2:6])
        colnames(fh) <- c("chromosome", "start", "end", "siteinf", "FAR", "yhat", "type")
        fh$chromosome <- paste0("chr", fh$chromosome)
    } else {
        fh = data.frame(chromosome=character(), start=numeric(), end=numeric(), siteinf=character(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    }

    return(fh)
}

Get_logCov <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    colnames(fh) <- c("chromosome", "start", "end", "logcov", "logcov.pcf")
    fh$chromosome <- paste0("chr", fh$chromosome)
    return(fh)
}

Get_HomoDiffest <- function(homodiffest) {
    fh <- read.delim(homodiffest, header=TRUE, sep="\t")
    pos1 <- fh$POS-1
    fh <- cbind(fh[,1], pos1, fh[,2:5])
    colnames(fh) <- c("chromosome", "start", "end", "PH1", "FFest", "FFest.pcf")
    fh$chromosome <- paste0("chr", fh$chromosome)
    return(fh)
}

#Get_parentalFH <- function(infh) {
#    fh = read.delim(infh, header=FALSE, sep="\t", comment.char = "#", stringsAsFactors=FALSE)
#    subfh = fh[,c(1,2,7,8,9,16)] 
#    colnames(subfh) = c("chromosome","end","PATRAF","MATRAF","SIBRAF","type")
#    start = subfh$end - 1
#    pfh = cbind(subfh[,1], start, subfh[,2:6])
#    colnames(pfh) = c("chromosome","start","end","PATRAF","MATRAF","SIBRAF","type")
#	pfh$chromosome <- paste0("chr", pfh$chromosome)
#    return(pfh)
#}


Plot_gtrellis_singlechr <- function(patfh, matfh, logcovfh, homodiffest, patlogcovfh, matlogcovfh, samplename, figdir, ffest_scale, scriptdir, ffsex) {
    #chrom = as.character(chrom)
    chroms <- seq(1,22,1)
    chroms <- paste0("chr", chroms)
    ffest_scale <- as.numeric(ffest_scale)
    ffest_scale_raw <- ffest_scale*5
    ffest_scale_fit <- ffest_scale*2
    patfh <- Get_SNPratio(patfh)
    if (nrow(patfh) == 0) {
        patfh <- rbind(patfh, data.frame(chromosome=rep(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"), 2), start=rep(1, 46), end=rep(2, 46), siteinf=c(rep("PH1", 23), rep("PH2", 23)), FAR=rep(0, 46), yhat=rep(0, 46), type=c(rep("P1", 23), rep("P2", 23))))
    }

    matfh <- Get_SNPratio(matfh)
    if (nrow(matfh) == 0) {
        matfh <- rbind(matfh, data.frame(chromosome=rep(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"), 2), start=rep(1, 46), end=rep(2, 46), siteinf=c(rep("MH1", 23), rep("MH2", 23)), FAR=rep(0, 46), yhat=rep(0, 46), type=c(rep("M1", 23), rep("M2", 23))))
    }
    
    logcovfh <- Get_logCov(logcovfh)
    homodiffest <- Get_HomoDiffest(homodiffest)
    patlogcovfh <- Get_logCov(patlogcovfh)
    matlogcovfh <- Get_logCov(matlogcovfh)
    chrsizeinfo <- Get_hg19_chrsize(scriptdir)
    val <- c(-ffest_scale, 0, ffest_scale)
    valtri <- c(-ffest_scale*1.5, -ffest_scale*0.5, ffest_scale*0.5, ffest_scale*1.5)

    outdir <- paste0(figdir, "/chr/")
    if (dir.exists(outdir) != TRUE) {
        dir.create(outdir)
    }
	
    fo <- paste0(outdir, samplename, ".gtrellis_chr.dev.pdf")
    pdf(file=fo, width=7.5, height=4.5)
    for (ic in chroms) {
	chrsisub <- subset(chrsizeinfo, chromosome==ic)
	#chrsisub[3,] <- chrsisub[2,] <- chrsisub[1,]
        #chrsisubnorm <- chrsisub[rep(seq_len(nrow(chrsisub)), each=3),]
        chrsisubnorm <- rbind(chrsisub, chrsisub, chrsisub)
        chrsisubnorm <- cbind(chrsisubnorm, val)
	chrsisubtri <- rbind(chrsisub, chrsisub, chrsisub, chrsisub)
	chrsisubtri <- cbind(chrsisubtri, valtri)
        patfhsub <- subset(patfh, chromosome==ic)
	patfhdens <- genomicDensity(patfhsub[,1:4], 1e6)
	patfhP1 <- subset(patfhsub, type=="P1")
	patfhP2 <- subset(patfhsub, type=="P2")
	
	matfhsub <- subset(matfh, chromosome==ic)
	matfhdens <- genomicDensity(matfhsub[,1:4], 1e6)
	matfhM1 <- subset(matfhsub, type=="M1")
	matfhM2 <- subset(matfhsub, type=="M2")
	
	logcovfhsub <- subset(logcovfh, chromosome==ic)
	patlogcovsub <- subset(patlogcovfh, chromosome==ic)
	matlogcovsub <- subset(matlogcovfh, chromosome==ic)
	        #covupper <- max(as.numeric(logcovfh$logcov))
	        #covlower <- min(as.numeric(logcovfh$logcov))
	        #print(covupper)
	homodiffestsub <- subset(homodiffest, chromosome==ic)
		#col_fun = circlize::colorRamp2(seq(0, covupper, length = 11), rev(brewer.pal(11, "RdYlBu")))
	gtrellis_layout(n_track = 10, category=ic, track_axis = TRUE, border=FALSE, padding = unit(c(2, 2, 2, 2), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 0.00015, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 0.00015, 0, ffest_scale_fit, -3, 3, -1, 1, -1, 1), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null")), track_ylab = c("PatRaw", "PatFit", "PatSNP", "MatRaw", "MatFit", "MatSNP", "FFest", "cfDNA.log2(Cov)", "Pat.log2(Cov)", "Mat.log2(Cov)"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=9)
	add_points_track(patfhP1, patfhP1$FAR, pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(patfhP1$siteinf=="PH2", "#c82027", "#344D90")))
	add_points_track(patfhP2, patfhP2$FAR, track = current_track(), pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(patfhP2$siteinf=="PH2", "#c82027", "#344D90"))) 
        add_segments_track(chrsisubnorm, chrsisubnorm$val, gp=gpar(col="#d9d9d9", lty=2, lwd=0.6))
	add_segments_track(chrsisubtri, chrsisubtri$valtri, track=current_track(), gp=gpar(col="#d4b9da", lty=3, lwd=0.6))
        add_points_track(patfhP1, patfhP1$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#c82027"))
	add_points_track(patfhP2, patfhP2$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#344D90"))
	add_lines_track(patfhdens, patfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
	add_points_track(matfhM1, matfhM1$FAR, pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(matfhM1$siteinf=="MH1", "#efb509", "#3f681c")))
	add_points_track(matfhM2, matfhM2$FAR, track = current_track(), pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(matfhM2$siteinf=="MH1", "#efb509", "#3f681c")))
	add_segments_track(chrsisubnorm, chrsisubnorm$val, gp=gpar(col="#d9d9d9", lty=2, lwd=0.6))	
        add_segments_track(chrsisubtri, chrsisubtri$valtri, track=current_track(), gp=gpar(col="#d4b9da", lty=3, lwd=0.6))
	add_points_track(matfhM1, matfhM1$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#efb509"))
	add_points_track(matfhM2, matfhM2$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#3f681c"))
	add_lines_track(matfhdens, matfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
	add_points_track(homodiffestsub, homodiffestsub$FFest, pch=1, size=unit(0.5, "mm"), gp=gpar(col="#d9d9d9"))
	add_points_track(homodiffestsub, homodiffestsub$FFest.pcf, track = current_track(), pch=16, size=unit(1, "mm"), gp=gpar(col="#595959"))
		#add_lines_track(logcovfhsub, logcovfhsub$logcov, gp=gpar(col = "#a6cee3"))
	add_rect_track(logcovfhsub, h1=logcovfhsub$logcov, h2=0, gp=gpar(fill = "#92c5de", col = NA))
	add_points_track(logcovfhsub, logcovfhsub$logcov.pcf, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
	add_points_track(patlogcovsub, patlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
	add_points_track(matlogcovsub, matlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
    }
	
    #plot sex chromosomes
    chrsisub <- subset(chrsizeinfo, chromosome=="chrX")
    chrsisubnorm <- rbind(chrsisub, chrsisub, chrsisub)
    chrsisubnorm <- cbind(chrsisubnorm, val)
    chrsisubtri <- rbind(chrsisub, chrsisub, chrsisub, chrsisub)
    chrsisubtri <- cbind(chrsisubtri, valtri)
    matfhsub <- subset(matfh, chromosome=="chrX")
    matfhdens <- genomicDensity(matfhsub[,1:4], 1e6)
    logcovfhsub <- subset(logcovfh, chromosome=="chrX")
    patlogcovsub <- subset(patlogcovfh, chromosome=="chrX")
    matlogcovsub <- subset(matlogcovfh, chromosome=="chrX")
    homodiffestsub <- subset(homodiffest, chromosome=="chrX")

    if (ffsex=="male") {
        matfhM1 <- subset(matfh, type=="M1" & chromosome=="chrX")
        gtrellis_layout(n_track = 7, category="chrX", track_axis = TRUE, border=FALSE, padding = unit(c(2, 2, 2, 2), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, range(matfhdens[,4]), 0, ffest_scale_fit, -3, 3, -1, 1, -1, 1), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null")), track_ylab = c("MatRaw", "MatFit", "MatSNP", "FFest", "cfDNA.log2(Cov)", "Pat.log2(Cov)", "Mat.log2(Cov)"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=9)
        add_points_track(matfhM1, matfhM1$FAR, pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(matfhM1$siteinf=="MH1", "#efb509", "#3f681c")))
	add_segments_track(chrsisubnorm, chrsisubnorm$val, gp=gpar(col="#d9d9d9", lty=2, lwd=0.6))
        add_segments_track(chrsisubtri, chrsisubtri$valtri, track=current_track(), gp=gpar(col="#d4b9da", lty=3, lwd=0.6))
	add_points_track(matfhM1, matfhM1$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#efb509"))
	add_lines_track(matfhdens, matfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
        add_points_track(homodiffestsub, homodiffestsub$FFest, pch=1, size=unit(0.5, "mm"), gp=gpar(col="#d9d9d9"))
        add_points_track(homodiffestsub, homodiffestsub$FFest.pcf, track = current_track(), pch=16, size=unit(1, "mm"), gp=gpar(col="#595959"))
        add_rect_track(logcovfhsub, h1=logcovfhsub$logcov, h2=0, gp=gpar(fill = "#92c5de", col = NA))
        add_points_track(logcovfhsub, logcovfhsub$logcov.pcf, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
        add_points_track(patlogcovsub, patlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
        add_points_track(matlogcovsub, matlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
    } else {
        matfhM1 <- subset(matfh, type=="M1" & chromosome=="chrX")
        matfhM2 <- subset(matfh, type=="M2" & chromosome=="chrX")
        gtrellis_layout(n_track = 7, category="chrX", track_axis = TRUE, border=FALSE, padding = unit(c(2, 2, 2, 2), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, range(matfhdens[,4]), 0, ffest_scale_fit, -3, 3, -1, 1, -1, 1), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null")), track_ylab = c("MatRaw", "MatFit", "MatSNP", "FFest", "cfDNA.log2(Cov)", "Pat.log2(Cov)", "Mat.log2(Cov)"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=9)
        add_points_track(matfhM1, matfhM1$FAR, pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(matfhM1$siteinf=="MH1", "#efb509", "#3f681c")))
        add_points_track(matfhM2, matfhM2$FAR, track = current_track(), pch=1, size=unit(0.8, "mm"), gp=gpar(col=ifelse(matfhM2$siteinf=="MH1", "#efb509", "#3f681c")))
        add_segments_track(chrsisubnorm, chrsisubnorm$val, gp=gpar(col="#d9d9d9", lty=2, lwd=0.6))
        add_segments_track(chrsisubtri, chrsisubtri$valtri, track=current_track(), gp=gpar(col="#d4b9da", lty=3, lwd=0.6))
        add_points_track(matfhM1, matfhM1$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#efb509"))
        add_points_track(matfhM2, matfhM2$yhat, track = current_track(), pch=16, size=unit(0.8, "mm"), gp=gpar(col="#3f681c"))
        add_lines_track(matfhdens, matfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
        add_points_track(homodiffestsub, homodiffestsub$FFest, pch=1, size=unit(0.5, "mm"), gp=gpar(col="#d9d9d9"))
        add_points_track(homodiffestsub, homodiffestsub$FFest.pcf, track = current_track(), pch=16, size=unit(1, "mm"), gp=gpar(col="#595959"))
        add_rect_track(logcovfhsub, h1=logcovfhsub$logcov, h2=0, gp=gpar(fill = "#92c5de", col = NA))
        add_points_track(logcovfhsub, logcovfhsub$logcov.pcf, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
        add_points_track(patlogcovsub, patlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
        add_points_track(matlogcovsub, matlogcovsub$logcov.pcf, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
    }	
    dev.off()
}

Plot_gtrellis_singlechr(patfh, matfh, logcovfh, homodiffest, patlogcovfh, matlogcovfh, samplename, figdir, ffest_scale, scriptdir, ffsex)
