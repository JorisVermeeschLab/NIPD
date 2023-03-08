args <- commandArgs(TRUE)
SegResolve <- args[1]
CBSseg <- args[2]
CBSsegsum <- args[3]
ffsex <- args[4]

source(SegResolve)

Resolve_Inheritance_CBSseg(CBSseg, CBSsegsum, ffsex)
