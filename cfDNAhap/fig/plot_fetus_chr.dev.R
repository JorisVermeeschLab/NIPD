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
logcovfh <- args[3]
homodiffest <- args[4]
patlogcovfh <- args[5]
matlogcovfh <- args[6]
samplename <- args[7]
workdir <- args[8]
ffest_scale <- args[9]
scriptdir <- args[10]
ffsex <- args[11]

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
	fh <- data.frame(chromosome=character(), start=numeric(), end=numeric(), siteinf=character(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
     }
    return(fh)
}

Get_logCov <- function(fh) {
    if (fh != "NULL") {
        fh <- read.delim(fh, header=TRUE, sep="\t")
        colnames(fh) <- c("INDEX", "CHROM", "START", "END", "MEAN.GCC.NORM", "Z", "LOG2R", "Seg.Z", "Seg.LOG2R")
        fh <- fh[,c("CHROM", "START", "END", "Z", "LOG2R", "Seg.Z", "Seg.LOG2R")]
        fh$CHROM <- paste0("chr", fh$CHROM)
    } else {
        fh <- data.frame(CHROM=character(), START=numeric(), END=numeric(), Z=numeric(), LOG2R=numeric(), Seg.Z=numeric(), Seg.LOG2R=numeric())
    }
    return(fh)
}

Get_HomoDiffest <- function(homodiffest) {
    fh <- read.delim(homodiffest, header=TRUE, sep="\t")
    pos1 <- fh$POS-1
    fh <- cbind(fh[,1], pos1, fh[,2:5])
    colnames(fh) <- c("chromosome", "start", "end", "PH1", "FFest", "FFest.pcf")
    fh$FFest <- 100*as.numeric(fh$FFest)
    fh$FFest.pcf <- 100*as.numeric(fh$FFest.pcf)
    fh$chromosome <- paste0("chr", fh$chromosome)
    return(fh)
}


Plot_gtrellis_singlechr <- function(patfh, matfh, logcovfh, homodiffest, patlogcovfh, matlogcovfh, samplename, workdir, ffest_scale, scriptdir, ffsex) {
    chroms <- seq(1,22,1)
    chroms <- paste0("chr", chroms)
    chroms <- c(chroms, "chrX")

    ffest_scale <- 100*as.numeric(ffest_scale)
    ffest_scale_raw <- ffest_scale*5
    ffest_scale_fit <- ffest_scale*2

    patfh <- Get_SNPratio(patfh)
    matfh <- Get_SNPratio(matfh)
    logcovfh <- Get_logCov(logcovfh)
    patlogcovfh <- Get_logCov(patlogcovfh)
    matlogcovfh <- Get_logCov(matlogcovfh)
    homodiffest <- Get_HomoDiffest(homodiffest)

    outdir <- paste0(workdir, "/6_fig/")
    if (dir.exists(outdir) != TRUE) {
        dir.create(outdir)
    }
    tstp <- format(Sys.time(), "%Y%m%d%H%M")
    fo <- paste0(outdir, samplename, ".", tstp, ".gtrellis_chr.dev.pdf")
    pdf(file=fo, width=7.5, height=4.5)

    for (ic in chroms) {
        patfhsub <- subset(patfh, chromosome==ic)

        matfhsub <- subset(matfh, chromosome==ic)

        logcovfhsub <- subset(logcovfh, CHROM==ic)
	patlogcovsub <- subset(patlogcovfh, CHROM==ic)
	matlogcovsub <- subset(matlogcovfh, CHROM==ic)
	homodiffestsub <- subset(homodiffest, chromosome==ic)

        gtrellis_layout(n_track=13, category=ic, xaxis = FALSE, xaxis_bin = 25000000, track_axis = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), padding = unit(c(2, 2, 2, 2), "mm"), track_ylim = c(0, 1, 0, 1, 0, 0.5, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 0.00015, -ffest_scale_raw, ffest_scale_raw, -ffest_scale_fit, ffest_scale_fit, 0, 0.00015, 0, ffest_scale_fit, -1.5, 1.5, -1.5, 1.5, -1.5, 1.5), track_height = unit.c(unit(0.5, "null"), unit(0.3, "null"), unit(0.5, "null"), unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(2.5, "null"), unit(2, "null"), unit(1, "null"), unit(1, "null"), unit(1, "null"), unit(1, "null"), unit(1, "null")), track_ylab = c("", "" , "", "PatRaw", "PatFit", "SNPden", "MatRaw", "MatFit", "SNPden", "FFest", "cfDNA.log2R", "Pat.log2R", "Mat.log2R"), axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=NULL, xlab=NULL, title_fontsize=9)
        #tr1
        add_track(panel_fun = function(gr) {
	    grid.text(c(ic), gp=gpar(fontsize=8))
        })

	#tr2
        add_ideogram_track(track=2)

	#tr3
        add_track(panel_fun = function(gr) {
	    grid.text("25MB", x=unit(25000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
            grid.text("50MB", x=unit(50000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("75MB", x=unit(75000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("100MB", x=unit(100000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("125MB", x=unit(125000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("150MB", x=unit(150000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("175MB", x=unit(175000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("200MB", x=unit(200000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
	    grid.text("225MB", x=unit(225000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=4))
        })

	#tr4
        if (nrow(patfhsub)!=0) {
	    add_track(patfhsub, panel_fun = function(gr) {
	        x <- patfhsub$start
	        y <- patfhsub$FAR
		grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$siteinf=="PH2", "#c82027", "#344D90")))
	    })
	} else {
	    add_track(panel_fun = function(gr) {

            })
	}

	#tr5
        if (nrow(patfhsub)!=0) {
	    add_track(patfhsub, panel_fun = function(gr) {
                x <- patfhsub$start
	        y <- patfhsub$yhat
                grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
                grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$type=="P1", "#c82027", "#344D90")))
	    })
	} else {
	    add_track(panel_fun = function(gr) {
	        grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	        grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            })
	}

	#tr6
	if (nrow(patfhsub)!=0) {
	    patfhdens <- genomicDensity(patfhsub[,1:4], 1e6)
	    add_lines_track(patfhdens, patfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
	} else {
	    add_track(panel_fun = function(gr){})
	}

        #maternal
	#tr7
	if (nrow(matfhsub)!=0) {
	    add_track(matfhsub, panel_fun = function(gr) {
	        x <- matfhsub$start
	        y <- matfhsub$FAR
	        grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(matfhsub$siteinf=="MH1", "#c82027", "#344D90")))
	    })
	} else {
	    add_track(panel_fun = function(gr) {})
	}

	#tr8
	if (nrow(matfhsub)!=0) {
	    add_track(matfhsub, panel_fun = function(gr) {
		x <- matfhsub$start
		y <- matfhsub$yhat
		grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(matfhsub$type=="M1", "#c82027", "#344D90")))
	    })
	} else {
	    add_track(panel_fun = function(gr) {
		grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
		grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	    })
	}

	#tr9
	if (nrow(matfhsub)!=0) {
	    matfhdens <- genomicDensity(matfhsub[,1:4], 1e6)
            add_lines_track(matfhdens, matfhdens[, 4], area = TRUE, gp = gpar(fill = "#999999", col = NA))
	} else {
            add_track(panel_fun = function(gr){})
	}

	#tr10
	if (nrow(homodiffestsub)!=0) {
	    add_points_track(homodiffestsub, homodiffestsub$FFest, pch=1, size=unit(0.5, "mm"), gp=gpar(col="#d9d9d9"))
	    add_points_track(homodiffestsub, homodiffestsub$FFest.pcf, track = current_track(), pch=16, size=unit(1, "mm"), gp=gpar(col="#595959"))
	} else {
	    add_track(panel_fun = function(gr){})
	}

	#tr11
	if (nrow(logcovfhsub)!=0) {
	    add_track(logcovfhsub, panel_fun = function(gr) {
	        x <- logcovfhsub$START
		y1 <- logcovfhsub$LOG2R
		y2 <- logcovfhsub$Seg.LOG2R
		grid.points(x, y1, pch=1, size=unit(0.2, "mm"), gp=gpar(col="#595959"))
		grid.points(x, y2, pch=16, size=unit(0.5, "mm"), gp=gpar(col="red"))
	    })
	} else {
	    add_track(panel_fun=function(gr) {})
	}

        #tr12
	if (nrow(patlogcovsub)!=0) {
	    add_track(patlogcovsub, panel_fun = function(gr) {
	        x <- patlogcovsub$START
		y1 <- patlogcovsub$LOG2R
		y2 <- patlogcovsub$Seg.LOG2R
		grid.points(x, y1, pch=1, size=unit(0.2,"mm"), gp=gpar(col="#595959"))
		grid.points(x, y2, pch=16, size=unit(0.5, "mm"), gp=gpar(col="red"))
	    })
	} else {
	    add_track(panel_fun=function(gr) {})
	}
	
	#tr13
	if (nrow(matlogcovsub)!=0) {
	    add_track(matlogcovsub, panel_fun = function(gr) {
	        x <- matlogcovsub$START
		y1 <- matlogcovsub$LOG2R
		y2 <- matlogcovsub$Seg.LOG2R
		grid.points(x, y1, pch=1, size=unit(0.2,"mm"), gp=gpar(col="#595959"))
		grid.points(x, y2, pch=16, size=unit(0.5, "mm"), gp=gpar(col="red"))
	    })
	} else {
	    add_track(panel_fun=function(gr) {})
	}

    }
    dev.off()

}


Plot_gtrellis_singlechr(patfh, matfh, logcovfh, homodiffest, patlogcovfh, matlogcovfh, samplename, workdir, ffest_scale, scriptdir, ffsex)
