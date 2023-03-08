args <- commandArgs(TRUE)
CBScore <- args[1]
infh <- args[2] 
sampleid <- args[3]
ffest <- args[4]
ffsex <- args[5]

source(CBScore)

FetalHapSeg(infh, sampleid, ffest, ffsex)
