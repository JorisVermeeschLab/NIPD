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
homodiffest = args[3]
cflogcov = args[4]
samplename = args[5]
workdir = args[6]
format = args[7]
ffest_scale = args[8]


Get_SNPratio = function(fh) {
    fh = read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 = fh$pos-1
        fh = cbind(fh[,1], pos1, fh[,2:6])
        colnames(fh) = c("chromosome", "start", "end", "siteinf", "FAR", "yhat", "type")
    } else {
        fh = data.frame(chromosome=character(), start=numeric(), end=numeric(), siteinf=character(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    }

    return(fh)
}

Get_HomoDiffest = function(homodiffest) {
    fh = read.delim(homodiffest, header=TRUE, sep="\t")
    pos1 = fh$POS-1
    fh = cbind(fh[,1], pos1, fh[,2:5])
    colnames(fh) = c("chromosome", "start", "end", "PH1", "FFest", "FFest.pcf")
    return(fh)
}

Get_cfDNAlogcovNOREF <- function(cflogcov) {
    fh <- read.delim(cflogcov, header=TRUE, sep="\t")
    colnames(fh) <- c("chromosome", "start", "end", "logcov", "logcov.pcf")
    return(fh)
}

Get_cfDNAlogcovREF <- function(cflogcov) {
    fh <- read.delim(cflogcov, header=TRUE, sep="\t")
    fh <- fh[,c("CHROM","START","END","Z","LOG2R","Seg.Z","Seg.LOG2R")]
    colnames(fh) <- c("chromosome","start","end","Z","LOG2R","SegZ","SegLOG2R")
    return(fh)
}

Plot_gtrellis_all_2col <- function(patfh, matfh, homodiffest, samplename, figdir, ffest_scale) {
    ffest_scale_raw <- as.numeric(ffest_scale)*4
    ffest_scale_fit <- as.numeric(ffest_scale)*2
    patfh <- Get_SNPratio(patfh)

    if (nrow(patfh) != 0) {
        patfh$chromosome <- paste0("chr", patfh$chromosome)
    } else {
        patfh <- rbind(patfh, data.frame(chromosome="chr1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="P1"), data.frame(chromosome="chr1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="P2"))
    }
    patfhP1 <- subset(patfh, type=="P1")
    patfhP2 <- subset(patfh, type=="P2")

    matfh <- Get_SNPratio(matfh)
    if (nrow(matfh) != 0) {
        matfh$chromosome <- paste0("chr", matfh$chromosome)
    } else {
        matfh <- rbind(matfh, data.frame(chromosome="chr1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="M1"), data.frame(chromosome="chr1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="M2"))
    }
    matfhM1 <- subset(matfh, type=="M1")
    matfhM2 <- subset(matfh, type=="M2")
    
    homodiffest <- Get_HomoDiffest(homodiffest)
    homodiffest$chromosome = paste0("chr", homodiffest$chromosome)
    
    fo = paste0(figdir, samplename, ".gtrellis_all.dev.pdf")
    pdf(file=fo, width=8.27, height=11.69)
    gtrellis_layout(n_track = 5, ncol = 2, byrow = FALSE, track_axis = TRUE, border=FALSE, padding = unit(c(1, 1, 1, 1), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, ffest_scale_fit), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(2.5, "null"), unit(2, "null"), unit(1, "null")), track_ylab = c("PatR", "PatF", "MatR", "MatF", "FFest"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=11)
	add_points_track(patfhP1, patfhP1$FAR, pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(patfhP1$siteinf=="PH2", "#c82027", "#344D90")))
	add_points_track(patfhP2, patfhP2$FAR, track = current_track(), pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(patfhP2$siteinf=="PH2", "#c82027", "#344D90")))
	add_points_track(patfhP1, patfhP1$yhat, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#c82027"))
	add_points_track(patfhP2, patfhP2$yhat, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#344D90"))
	add_points_track(matfhM1, matfhM1$FAR, pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(matfhM1$siteinf=="MH1", "#c82027", "#344D90")))
	add_points_track(matfhM2, matfhM2$FAR, track = current_track(), pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(matfhM2$siteinf=="MH1", "#c82027", "#344D90")))
	add_points_track(matfhM1, matfhM1$yhat, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#c82027"))
	add_points_track(matfhM2, matfhM2$yhat, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#344D90"))
	add_points_track(homodiffest, homodiffest$FFest, pch=1, size=unit(0.05, "mm"), gp=gpar(col="#d9d9d9"))
	add_points_track(homodiffest, homodiffest$FFest.pcf, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
     dev.off()				
}

Plot_gtrellis_all_1row = function(patfh, matfh, homodiffest, cflogcov, samplename, figdir, format="png", ffest_scale) {
    ffest_scale_raw = as.numeric(ffest_scale)*4
    ffest_scale_fit = as.numeric(ffest_scale)*2
    patfh = Get_SNPratio(patfh)
    if (nrow(patfh) == 0) {
        patfh <- rbind(patfh, data.frame(chromosome="1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="P1"), data.frame(chromosome="1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="P2"))
    }
    patfhP1 = subset(patfh, type=="P1")
    patfhP2 = subset(patfh, type=="P2")

    matfh = Get_SNPratio(matfh)
    if (nrow(matfh) == 0) {
        matfh <- rbind(matfh, data.frame(chromosome="1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="M1"), data.frame(chromosome="1", start=1, end=2, siteinf="NA", FAR=0, yhat=0, type="M2"))
    }
    matfhM1 = subset(matfh, type=="M1")
    matfhM2 = subset(matfh, type=="M2")

    homodiffest = Get_HomoDiffest(homodiffest)
    cflogR = Get_cfDNAlogcovREF(cflogcov)

    if (format=="png") {
        fo = paste0(figdir, samplename, ".gtrellis_all_horiz.", format)
        png(file=fo, width=1400, height=600, res=160)
    } else {
        fo = paste0(figdir, samplename, ".gtrellis_all_horiz.pdf")
        pdf(file=fo, width=11.7, height=3)
    }

    gtrellis_layout(n_track = 6, track_axis = TRUE, border = FALSE, remove_chr_prefix = TRUE, padding = unit(c(1, 1, 1, 1), "mm"), track_ylim = c(-ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, ffest_scale_fit, -3, 3), track_height = unit.c(unit(2.5, "null"), unit(2, "null"), unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(2.5, "null")), track_ylab = c("PatR", "PatF", "MatR", "MatF", "FFest", "cfDNA(logCov)"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=11)
	add_points_track(patfhP1, patfhP1$FAR, pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(patfhP1$siteinf=="PH2", "#c82027", "#344D90")))
	add_points_track(patfhP2, patfhP2$FAR, track = current_track(), pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(patfhP2$siteinf=="PH2", "#c82027", "#344D90")))
	add_points_track(patfhP1, patfhP1$yhat, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#c82027"))
	add_points_track(patfhP2, patfhP2$yhat, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#344D90"))
	add_points_track(matfhM1, matfhM1$FAR, pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(matfhM1$siteinf=="MH1", "#c82027", "#344D90")))
	add_points_track(matfhM2, matfhM2$FAR, track = current_track(), pch=1, size=unit(0.05, "mm"), gp=gpar(col=ifelse(matfhM2$siteinf=="MH1", "#c82027", "#344D90")))
	add_points_track(matfhM1, matfhM1$yhat, pch=16, size=unit(0.5, "mm"), gp=gpar(col="#c82027"))
	add_points_track(matfhM2, matfhM2$yhat, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#344D90"))
	add_points_track(homodiffest, homodiffest$FFest, pch=1, size=unit(0.05, "mm"), gp=gpar(col="#d9d9d9"))
	add_points_track(homodiffest, homodiffest$FFest.pcf, track = current_track(), pch=16, size=unit(0.5, "mm"), gp=gpar(col="#595959"))
        add_points_track(cflogR, cflogR$LOG2R, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#d9d9d9"))
	add_points_track(cflogR, cflogR$SegLOG2R, track=current_track(), pch=16, size=unit(0.1, "mm"), gp=gpar(col="#b30000"))
    dev.off()				
}

Plot_gtrellis_main <- function(patfh, matfh, homodiffest, cflogcov, samplename, workdir, format, ffest_scale) {
    figdir = paste0(workdir, "6_fig/")
    if (dir.exists(figdir) != TRUE) {
	dir.create(figdir)
    }
    Plot_gtrellis_all_2col(patfh, matfh, homodiffest, samplename, figdir, ffest_scale)
    Plot_gtrellis_all_1row(patfh, matfh, homodiffest, cflogcov, samplename, figdir, format, ffest_scale)
}

Plot_gtrellis_main(patfh, matfh, homodiffest, cflogcov, samplename, workdir, format, ffest_scale)

#Get_Copyratio = function(cr) {
#	cr = read.delim(cr, header=TRUE, sep="\t")
#	colnames(cr) = c("chromosome","start","end","value1", "value2")
#	return(cr)
#}
