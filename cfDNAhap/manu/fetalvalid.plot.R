.libPaths(c("/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/R/3.4.2",.libPaths()))

list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args = commandArgs(TRUE)
fetalfh = args[1]
cfpat20 = args[2]
cfmat20 = args[3]
cfpat10 = args[4]
cfmat10 = args[5]
cfpat5 = args[6]
cfmat5 = args[7]
workdir = args[8]
samplename = args[9]

readinfetalfh <- function(fetalfh) {
	fh <- read.delim(fetalfh, header=FALSE, sep="\t")
	colnames(fh) <- c("chromosome", "pos", "rsid", "ref", "alt", "FPH1", "FPH2", "FMH1", "FMH2", "FT1", "FT2", "pathap", "mathap")
	start0 <- fh$pos-1
	ft <- cbind(fh[,1], start0, fh[,c(2,12,13)])
	colnames(ft) <- c("chromosome", "start", "end", "pathap", "mathap")
	ft$chromosome <- paste0("chr", ft$chromosome)
	return(ft)
}

fetalpatimpute <- function(fetalfh) {
	fh <- readinfetalfh(fetalfh)
	ftpathap <- fh[,c(1,2,3,4)]
	imputepathap <- data.frame(chromosome=character(), start=numeric(), end=numeric(), pathap=character(), stringsAsFactors=FALSE)
	
	uchrom <- unique(ftpathap$chromosome)
	for (ic in uchrom) {
		dt.chr <- subset(ftpathap, chromosome == ic)
		for (i in 1:nrow(dt.chr)) {
			if (i == 1) {
				imputeval <- dt.chr$pathap
			} else if (dt.chr$pathap[i] != "PH1" & dt.chr$pathap[i] != "PH2") {
				dt.chr$pathap[i] == imputeval
			}
			imputeval <- dt.chr$pathap
		}
		imputepathap <- rbind(imputepathap, dt.chr)
	}
	imputepathap <- subset(imputepathap, pathap == "PH1" | pathap == "PH2")
	
	write.table(imputepathap, "/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/4_Jun2018/7_trisomy/GC063925_TTCACGCA_S3.fetal.genome.impute.pat.hap", row.names=FALSE,quote=FALSE, sep="\t")
	return(imputepathap)
}

fetalmatimpute <- function(fetalfh) {
	fh <- readinfetalfh(fetalfh)
	ftmathap <- fh[,c(1,2,3,5)]
	imputemathap <- data.frame(chromosome=character(), start=numeric(), end=numeric(), mathap=character(), stringsAsFactors=FALSE)
	
	uchrom <- unique(ftmathap$chromosome)
	for (ic in uchrom) {
		dt.chr <- subset(ftmathap, chromosome == ic)
		for (i in 1:nrow(dt.chr)) {
			if (i == 1) {
				imputeval <- dt.chr$mathap
			} else if (dt.chr$mathap[i] != "PH1" & dt.chr$mathap[i] != "PH2") {
				dt.chr$mathap[i] == imputeval
			}
			imputeval <- dt.chr$mathap
		}
		imputemathap <- rbind(imputemathap, dt.chr)
	}
	imputemathap <- subset(imputemathap, mathap == "MH1" | mathap == "MH2")
        write.table(imputemathap, "/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/4_Jun2018/7_trisomy/GC063925_TTCACGCA_S3.fetal.genome.impute.mat.hap", row.names=FALSE, quote=FALSE, sep="\t")
	return(imputemathap)
}


Get_SegRes <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    colnames(fh) <- c("chromosome", "start", "end", "resHL")
    fh$chromosome <- paste0("chr", fh$chromosome)
    fh <- subset(fh, resHL != "inconclusive")
    return(fh)
}

fetalvalidpatplot <- function(fetalfh, cfpat20, cfmat20, cfpat10, cfmat10, cfpat5, cfmat5, workdir, samplename) {
    ftpathap <- fetalpatimpute(fetalfh)
	pat20 <- Get_SegRes(cfpat20)
	#mat20 <- Get_SegRes(cfmat20)
	pat10 <- Get_SegRes(cfpat10)
	#mat10 <- Get_SegRes(cfmat10)
	pat5 <- Get_SegRes(cfpat5)
	#mat5 <- Get_SegRes(cfmat5)
    
    plotchr <- paste0("chr", seq(1,22))
    fo <- paste0(workdir, "/", samplename, ".fetalvalid.pat.pdf")
    pdf(file=fo, width=4.2, height=10)
    gtrellis_layout(n_track = 4, ncol = 1, category=plotchr, track_axis = FALSE, xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE, asist_ticks = FALSE, add_ideogram_track = TRUE, ideogram_track_height = unit(2, "mm"), xlab = NULL)
    add_rect_track(ftpathap, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(ftpathap$pathap=="PH2", "cornflowerblue", "blue")))
    add_rect_track(pat20, h1 = 0.5, h2 = 0, gp = gpar(col=NA, fill = ifelse(pat20$resHL=="PH2", "cornflowerblue", "blue")))
    add_rect_track(pat10, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(pat10$resHL=="PH2", "cornflowerblue", "blue")))
    add_rect_track(pat5, h1=0.5, h2=0,  gp = gpar(col=NA, fill = ifelse(pat5$resHL=="PH2", "cornflowerblue", "blue")))
    add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
        chr = get_cell_meta_data("name")
        grid.text(chr, x = 0, y = 0, just = c("left", "bottom"), gp=gpar(fontsize=8))
    })
    dev.off()
}

fetalvalidmatplot <- function(fetalfh, cfpat20, cfmat20, cfpat10, cfmat10, cfpat5, cfmat5, workdir, samplename) {
    ftmathap <- fetalmatimpute(fetalfh)

	#pat20 <- Get_SegRes(cfpat20)
	mat20 <- Get_SegRes(cfmat20)
	#pat10 <- Get_SegRes(cfpat10)
	mat10 <- Get_SegRes(cfmat10)
	#pat5 <- Get_SegRes(cfpat5)
	mat5 <- Get_SegRes(cfmat5)
    
    plotchr <- paste0("chr", seq(1,20))
	plotchr <- c(plotchr, "chr22", "chrX")
    fo <- paste0(workdir, "/", samplename, ".fetalvalid.mat.pdf")
    pdf(file=fo, width=4.2, height=10)
    gtrellis_layout(n_track = 4, ncol = 1, category=plotchr, track_axis = FALSE, xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE, asist_ticks = FALSE, add_ideogram_track = TRUE, ideogram_track_height = unit(2, "mm"),xlab = NULL)
    add_rect_track(ftmathap, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(ftmathap$mathap=="MH1", "red", "pink")))
    add_rect_track(mat20, h1 = 0.5, h2 = 0, gp = gpar(col=NA, fill = ifelse(mat20$resHL=="MH1", "red", "pink")))
    add_rect_track(mat10, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(mat10$resHL=="MH1", "red", "pink")))
    add_rect_track(mat5, h1=0.5, h2=0,  gp = gpar(col=NA, fill = ifelse(mat5$resHL=="MH1", "red", "pink")))
    add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
        chr = get_cell_meta_data("name")
        grid.text(chr, x = 0, y = 0, just = c("left", "bottom"),gp=gpar(fontsize=8))
    })
    dev.off()
}

fetalvalidpatplot(fetalfh, cfpat20, cfmat20, cfpat10, cfmat10, cfpat5, cfmat5, workdir, samplename)
fetalvalidmatplot(fetalfh, cfpat20, cfmat20, cfpat10, cfmat10, cfpat5, cfmat5, workdir, samplename)
