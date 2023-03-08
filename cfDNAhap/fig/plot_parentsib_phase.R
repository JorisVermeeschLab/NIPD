list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)


args = commandArgs(TRUE)
hetpatfh = args[1]
hetmatfh = args[2]
hetbothfh = args[3]
homodifffh = args[4]
hetpatbpfh = args[5]
hetmatbpfh = args[6]
patlogR = args[7]
matlogR = args[8]
siblogR = args[9]
samplename = args[10]
figdir = args[11]


ReadInFH <- function(infh) {
    fh = read.delim(infh, header=FALSE, sep="\t", comment.char = "#", stringsAsFactors=FALSE)
    subfh = fh[,c(1,2,7,8,9,16)] 
    colnames(subfh) = c("chromosome","end","PATRAF","MATRAF","SIBRAF","type")
    start = subfh$end - 1
    pfh = cbind(subfh[,1], start, subfh[,2:6])
    colnames(pfh) = c("chromosome","start","end","PATRAF","MATRAF","SIBRAF","type")
    return(pfh)
}

ReadINSIBbp <- function(infh) {
    fh = read.delim(infh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    start = fh$pos - 1
    pfh = cbind(fh[,1], start, fh[,2:5])
    colnames(pfh) = c("chromosome","start","end","SIBRAF","yhat","type")
    return(pfh)
}

ReadINlogR <- function(infh) {
    fh = read.delim(infh, header=TRUE, sep="\t")
    fh$end = fh$start + 1
    colnames(fh) = c("chromosome","start","end","logcov", "logcovseg")
    return(fh)
}


Plot_parentsib_2col <- function(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, samplename, figdir) {
    hetpat = ReadInFH(hetpatfh)
    hetmat = ReadInFH(hetmatfh)
    hetboth = ReadInFH(hetbothfh)
    homodiff = ReadInFH(homodifffh)
    hetpat$chromosome = paste0("chr", hetpat$chromosome)
    hetmat$chromosome = paste0("chr", hetmat$chromosome)
    hetboth$chromosome = paste0("chr", hetboth$chromosome)
    homodiff$chromosome = paste0("chr", homodiff$chromosome)
    hetpatbp = ReadINSIBbp(hetpatbpfh)
    hetmatbp = ReadINSIBbp(hetmatbpfh)
    hetpatbp$chromosome = paste0("chr", hetpatbp$chromosome)
    hetmatbp$chromosome = paste0("chr", hetmatbp$chromosome)
    hetpatbpP1 = subset(hetpatbp, type=="P1")
    hetpatbpP2 = subset(hetpatbp, type=="P2")
    hetmatbpM1 = subset(hetmatbp, type=="M1")
    hetmatbpM2 = subset(hetmatbp, type=="M2")
    fo = paste0(figdir, '/', samplename, ".parents_sibling_phasing.pdf")
    titlename = paste0(samplename, ".parents_phasing")
    pdf(file=fo, width=8.27, height=11.69)
    gtrellis_layout(n_track = 5, ncol = 2, byrow = FALSE, track_axis = TRUE, border=FALSE, padding = unit(c(1, 1, 1, 1), "mm"), track_ylim = c(-0.1, 1.1, -0.1, 1.1, -0.1, 1.1, -0.1, 1.1, -0.1, 1.1), track_height = unit.c(unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null")), track_ylab = c("PatRAF", "MatRAF", "SibRAF", "FRAF", "MRAF"), add_ideogram_track = TRUE, ideogram_track_height = unit(1, "mm"), add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=titlename, title_fontsize=11)
    add_points_track(hetpatbpP1, hetpatbpP1$SIBRAF, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#f6e0e0"))
    add_points_track(hetpatbpP2, hetpatbpP2$SIBRAF, track = current_track(), pch=1, size=unit(0.01, "mm"), gp=gpar(col="#e2e0f6"))
    add_points_track(hetpatbpP1, hetpatbpP1$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#c82027"))
    add_points_track(hetpatbpP2, hetpatbpP2$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#344D90"))
    add_points_track(hetmatbpM1, hetmatbpM1$SIBRAF, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#fef8e5"))
    add_points_track(hetmatbpM2, hetmatbpM2$SIBRAF, track = current_track(), pch=1, size=unit(0.01, "mm"), gp=gpar(col="#f2faeb"))
    add_points_track(hetmatbpM1, hetmatbpM1$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#efb509"))
    add_points_track(hetmatbpM2, hetmatbpM2$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#3f681c"))
    add_points_track(hetpat, hetpat$SIBRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetpat, hetpat$PATRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetpat, hetpat$MATRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    dev.off()				
}

Plot_parentsib_1row <- function(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, patlogR, matlogR, siblogR, samplename, figdir) {
    hetpat = ReadInFH(hetpatfh)
    hetmat = ReadInFH(hetmatfh)
    hetboth = ReadInFH(hetbothfh)
    homodiff = ReadInFH(homodifffh)
    hetpatbp = ReadINSIBbp(hetpatbpfh)
    hetmatbp = ReadINSIBbp(hetmatbpfh)
    patlogr = ReadINlogR(patlogR)
    matlogr = ReadINlogR(matlogR)
    siblogr = ReadINlogR(siblogR)
    hetpatbpP1 = subset(hetpatbp, type=="P1")
    hetpatbpP2 = subset(hetpatbp, type=="P2")
    hetmatbpM1 = subset(hetmatbp, type=="M1")
    hetmatbpM2 = subset(hetmatbp, type=="M2")
    fo = paste0(figdir, '/', samplename, ".parents_sibling_phasing.full.png")
    titlename = paste0(samplename, ".parents_phasing")
    png(file=fo, width=1600, height=650, res=150)
    gtrellis_layout(n_track = 8, track_axis = TRUE, remove_chr_prefix = TRUE, border=FALSE, padding = unit(c(1, 1, 1, 1), "mm"), track_ylim = c(-0.1, 1.1, -0.1, 1.1, -0.1, 1.1, -2, 2, -0.1, 1.1, -2, 2, -0.1, 1.1, -2, 2), track_height = unit.c(unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "null")), track_ylab = c("PatRAF", "MatRAF", "SibRAF", "SibLogR", "FRAF", "FLogR", "MRAF", "MLogR"), add_ideogram_track = TRUE, ideogram_track_height = unit(1, "mm"), add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=titlename, title_fontsize=9)
    add_points_track(hetpatbpP1, hetpatbpP1$SIBRAF, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#f6e0e0"))
    add_points_track(hetpatbpP2, hetpatbpP2$SIBRAF, track = current_track(), pch=1, size=unit(0.01, "mm"), gp=gpar(col="#e2e0f6"))
    add_points_track(hetpatbpP1, hetpatbpP1$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#c82027"))
    add_points_track(hetpatbpP2, hetpatbpP2$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#344D90"))
    add_points_track(hetmatbpM1, hetmatbpM1$SIBRAF, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#fef8e5"))
    add_points_track(hetmatbpM2, hetmatbpM2$SIBRAF, track = current_track(), pch=1, size=unit(0.01, "mm"), gp=gpar(col="#f2faeb"))
    add_points_track(hetmatbpM1, hetmatbpM1$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#efb509"))
    add_points_track(hetmatbpM2, hetmatbpM2$yhat, track = current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#3f681c"))
    add_points_track(hetpat, hetpat$SIBRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$SIBRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(siblogr, siblogr$logcov, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#adadad"))
    add_points_track(siblogr, siblogr$logcovseg, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="red"))
    add_points_track(hetpat, hetpat$PATRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$PATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(patlogr, patlogr$logcov, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#adadad"))
    add_points_track(patlogr, patlogr$logcovseg, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="red"))
    add_points_track(hetpat, hetpat$MATRAF, pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetmat, hetmat$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(hetboth, hetboth$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(homodiff, homodiff$MATRAF, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="#595959"))
    add_points_track(matlogr, matlogr$logcov, pch=1, size=unit(0.01, "mm"), gp=gpar(col="#adadad"))
    add_points_track(matlogr, matlogr$logcovseg, track=current_track(), pch=1, size=unit(0.03, "mm"), gp=gpar(col="red"))
    dev.off()				
}


Plot_parentsib_main <- function(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, patlogR, matlogR, siblogR, samplename, figdir) {
    Plot_parentsib_2col(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, samplename, figdir)
    Plot_parentsib_1row(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, patlogR, matlogR, siblogR, samplename, figdir)
}

Plot_parentsib_main(hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatbpfh, hetmatbpfh, patlogR, matlogR, siblogR, samplename, figdir)
