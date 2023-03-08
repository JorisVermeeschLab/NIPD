.libPaths(c("/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/R/3.4.2",.libPaths()))

list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args <- commandArgs(TRUE)
patfh <- args[1]
matfh <- args[2]
patres <- args[3]
matres <- args[4]
samplename <- args[5]
ffest_scale <- args[6] #input percent - if FF 10%, then input 10
chrom <- args[7] #chr1
figdir <- "/staging/leuven/stg_00019/research/Huiwen/project/1_Haplotype/manu/fig2/"

Get_SNPratio <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 <- fh$pos-1
        fh <- cbind(fh[,1], pos1, fh[,2:6])
        colnames(fh) <- c("chromosome", "start", "end", "siteinf", "FAR", "yhat", "type")
        fh$chromosome <- paste0("chr", fh$chromosome)
	fh$FAR <- 100*as.numeric(fh$FAR)
	fh$yhat <- 100*as.numeric(fh$yhat)
    } else {
        fh = data.frame(chromosome=character(), start=numeric(), end=numeric(), siteinf=character(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    }

    return(fh)
}

Get_SegRes <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    colnames(fh) <- c("chromosome", "start", "end", "resHL")
    fh$chromosome <- paste0("chr", fh$chromosome)
    fh <- subset(fh, resHL != "inconclusive")
    return(fh)
}

singlechrsib <- function(patfh, matfh, patres, matres, samplename, ffest_scale, chrom) {
    ffest_scale <- as.numeric(ffest_scale)
    ffest_scale_raw <- ffest_scale*5
    ffest_scale_fit <- ffest_scale*2

    patfh <- Get_SNPratio(patfh)
    matfh <- Get_SNPratio(matfh)

    patres <- Get_SegRes(patres)
    matres <- Get_SegRes(matres)

    fo <- paste0(figdir, samplename, ".", chrom, ".manu.pdf")
    pdf(file=fo, width=7, height=4.5)

    patfhsub <- subset(patfh, chromosome==chrom)
    matfhsub <- subset(matfh, chromosome==chrom)
    patressub <- subset(patres, chromosome==chrom)
    matressub <- subset(matres, chromosome==chrom)
    
    gtrellis_layout(n_track = 6, category=chrom, track_axis = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE), border=TRUE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 1, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 1), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(0.3, "null"), unit(2.5, "null"), unit(2, "null"), unit(0.3, "null")), track_ylab = c("Pat.FAR (%)", "Pat.FAR.seg (%)", "", "Mat.FAR (%)", "Mat.FAR.seg (%)", ""), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 5, lab_fontsize = 6, name_fontsize = 8, title=samplename, title_fontsize = 10, xlab = NULL)

    add_track(patfhsub, panel_fun = function(gr) {
        x <- patfhsub$start
	y <- patfhsub$FAR
        grid.points(x, y, pch = 1, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$siteinf=="PH2", "#c82027", "#344D90")))
    })

    add_track(patfhsub, panel_fun = function(gr) {
        x <- patfhsub$start
        y <- patfhsub$yhat
	grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
        grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.6))
        grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
        grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$type=="P1", "#c82027", "#344D90")))
    })

    add_rect_track(patressub, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(patressub$resHL=="PH2", "cornflowerblue", "blue")))

    add_track(matfhsub, panel_fun = function(gr) {
        x <- matfhsub$start
	y <- matfhsub$FAR
	grid.points(x, y, pch = 1, size = unit(0.8, "mm"), gp = gpar(col=ifelse(matfhsub$siteinf=="MH1", "#c82027", "#344D90")))
    })
    
    add_track(matfhsub, panel_fun = function(gr) {
        x <- matfhsub$start
	y <- matfhsub$yhat
	grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.6))
        grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(matfhsub$type=="M1", "#c82027", "#344D90")))
    })

    add_rect_track(matressub, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(matressub$resHL=="MH1", "red", "pink")))

    dev.off()
}

singlechrsib(patfh, matfh, patres, matres, samplename, ffest_scale, chrom)

singlechrpatgrand <- function(patfh, patres, samplename, ffest_scale, chrom) {
    ffest_scale <- as.numeric(ffest_scale)
    ffest_scale_raw <- ffest_scale*5
    ffest_scale_fit <- ffest_scale*2

    patfh <- Get_SNPratio(patfh)

    patres <- Get_SegRes(patres)

    fo <- paste0(figdir, samplename, ".", chrom, ".manu.pdf")
    pdf(file=fo, width=7, height=3)

    patfhsub <- subset(patfh, chromosome==chrom)
    patressub <- subset(patres, chromosome==chrom)
   
    gtrellis_layout(n_track = 3, category=chrom, track_axis = c(TRUE, TRUE, FALSE), border=TRUE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 1), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(0.3, "null")), track_ylab = c("Pat.FAR (%)", "Pat.FAR.seg (%)", ""), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 5, lab_fontsize = 6, name_fontsize = 8, title=samplename, title_fontsize = 10, xlab = NULL)

    add_track(patfhsub, panel_fun = function(gr) {
        x <- patfhsub$start
        y <- patfhsub$FAR
        grid.points(x, y, pch = 1, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$siteinf=="PH2", "#c82027", "#344D90")))
    })

    add_track(patfhsub, panel_fun = function(gr) {
        x <- patfhsub$start
        y <- patfhsub$yhat
        grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
        grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.6))
        grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
        grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$type=="P1", "#c82027", "#344D90")))
    })

    add_rect_track(patressub, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(patressub$resHL=="PH2", "cornflowerblue", "blue")))

    dev.off()
}

#singlechrpatgrand(patfh, patres, samplename, ffest_scale, chrom)
