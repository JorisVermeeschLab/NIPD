script.dir <- dirname(sys.frame(1)$ofile)
createrefsp <- file.path(script.dir, "createRef.R")
source(createrefsp)

args <- commandArgs(TRUE)
samplelist <- args[1]
targetgcfh <- args[2]
outdir <- args[3]
setname <- args[4]

GetCombinedCov(samplelist, targetgcfh, outdir, setname)
