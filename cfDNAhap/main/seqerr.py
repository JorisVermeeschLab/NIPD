#!/usr/bin/env python

'''
Author: Huiwen Che
Created: Mar 09 2018
Use: 
estimate sequencing and mapping error
'''

from __future__ import division
import sys
import os
from Helper import mkdir, check_dirslash, basename
from Utility_stats import cal_basic_stats
import numpy as np

def errest_homovar(homosame, aerc, mindepth, maxdepth, workdir):
    pref = basename(homosame)
    workdir = check_dirslash(workdir)
    parentshapdir = workdir + '1_parentshap/'
    errstats = '%s%s.homovar.err.summary' % (parentshapdir, pref)
    
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)

    homovar = dict()
    (allrefcnt, alltotalcnt, raf) = (0, 0, 0)
    (raffull, rafnonzero) = ([], [])
    
    with open(homosame, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                val = 0
                homovar[refID] = val
                
    with open(errstats, 'w+') as fo:
        fo.write('#CHROM\tPOS\tRSID\tREF\tALT\trefCOUNT\taltCOUNT\ttotalCOUNT\tRAF\n')
        with open(aerc, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    checkID = '%s:%s' % (contig, position)
                    if checkID in homovar:
                        #examine whether reference alleles present in the snp site, which indicates potential sequencing errors and mapping errors
                        if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                            allrefcnt += int(refCount)
                            alltotalcnt += int(totalCount)
                            raf = int(refCount)/int(totalCount)
                            raffull.append(raf)
                            if raf > 0:
                                rafnonzero.append(raf)
                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, raf)
                                fo.write(pfline)
                            
            raffullstat = cal_basic_stats(raffull)
            rafnonzerostat = cal_basic_stats(rafnonzero)
            overallraf = allrefcnt/alltotalcnt
            fo.write('#########\n#Summary Stats:\n')
            fo.write('Overall error rate\t%s\n' % overallraf)
            fo.write('#ALL reference allele frequency\n#Mean\tS.D.\tMedian\tMax\tMin\n')
            print >>fo, '\t'.join(map(str, raffullstat))
            fo.write('#Non zero reference allele frequency\n#Mean\tS.D.\tMedian\tMax\tMin\n')
            print >>fo, '\t'.join(map(str, rafnonzerostat))


if __name__ == '__main__':
    errest_homovar(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
