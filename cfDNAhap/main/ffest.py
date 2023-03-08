#!/usr/bin/env python

'''
use: calculate fetal fraction using different methods
'''

from __future__ import division
import sys
import os
from Helper import mkdir, check_dirslash, basename
from Utility_stats import cal_basic_stats
import numpy as np


def FFest_HomoDiffSNP(HomoDiff, cfAERC, ffdir, mindepth, maxdepth):
    ffdir = check_dirslash(ffdir)
    pref = basename(HomoDiff)
    ffhdout = '%s%s.ffest.out' % (ffdir, pref)
    snpsite = dict()
    (ffl, chrffl, tmpchrffl) = ([], [], [])
    pre_chr = '0'
    (autofetaltotal, autoalltotal, fetalsite, ff)  = (0, 0, 0, 0)
    
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    
    with open(HomoDiff, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                snpsite[refID] = PH1
                
    with open(ffhdout, 'w+') as fo:
        fo.write('#CHROM\tPOS\tPH1\tFFest\tDEPTH\n')
             
        with open(cfAERC, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    if (contig != 'X' and contig != 'Y'):
                        if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                            checkID = '%s:%s' % (contig, position)
                            if checkID in snpsite:
                                if contig != pre_chr and pre_chr != '0':
                                    chrfflstat = cal_basic_stats(tmpchrffl)
                                    chrfflstat.insert(0, pre_chr)
                                    chrffl.append(chrfflstat)
                                    tmpchrffl = []
                                if snpsite[checkID] == '0':
                                    fetalsite = int(refCount)
                                    autofetaltotal += int(refCount)
                                elif snpsite[checkID] == '1':
                                    fetalsite = int(altCount)
                                    autofetaltotal += int(altCount)
                                ff = 2*fetalsite/int(totalCount)
                                ffl.append(ff)
                                autoalltotal += int(totalCount)
                                tmpchrffl.append(ff)
                                fo.write('%s\t%s\t%s\t%s\t%s\n' % (contig, position, snpsite[checkID], ff, totalCount))
                                pre_chr = contig
        
        chrfflstat = cal_basic_stats(tmpchrffl)
        chrfflstat.insert(0, pre_chr)
        chrffl.append(chrfflstat)                            
        totalFFest = 2*autofetaltotal/autoalltotal
        fflstat = cal_basic_stats(ffl)
        fflstat.insert(0, totalFFest)
        return fflstat, chrffl
        

def FFest_HeteroPatSNP(HeteroPat, cfAERC, ffdir, mindepth, maxdepth, mediancov):
    ffdir = check_dirslash(ffdir)
    pref = basename(HeteroPat)
    ffhdout = '%s%s.ffest.out' % (ffdir, pref)
    snpsite = dict()
    (ffl, chrffl, tmpchrffl) = ([], [], [])
    pre_chr = '0'
    (autofetaltotal, autoalltotal, fetalsite, ff)  = (0, 0, 0, 0)
    ffthresh = 0.03
    
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    mediancov = int(mediancov)
    
    with open(HeteroPat, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                snpsite[refID] = MH1
                
    with open(ffhdout, 'w+') as fo:
        fo.write('#CHROM\tPOS\tPH1\tFFest\tDEPTH\n')
             
        with open(cfAERC, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    checkID = '%s:%s' % (contig, position)
                    if checkID in snpsite:
                        if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                            if contig != pre_chr and pre_chr != '0':
                                chrfflstat = cal_basic_stats(tmpchrffl)
                                chrfflstat.insert(0, pre_chr)
                                chrffl.append(chrfflstat)
                                tmpchrffl = []
                            if snpsite[checkID] == '0':
                                fetalsite = int(altCount)
                            elif snpsite[checkID] == '1':
                                fetalsite = int(refCount)
                            ff = 2*fetalsite/int(totalCount)
                            pre_chr = contig
                            if contig != 'X' and contig != 'Y' and ff >= ffthresh:
                                ffl.append(ff)
                                autoalltotal += int(totalCount)
                                if snpsite[checkID] == '0':
                                    autofetaltotal += int(altCount)
                                elif snpsite[checkID] == '1':
                                    autofetaltotal += int(refCount)
                                tmpchrffl.append(ff)
                                fo.write('%s\t%s\t%s\t%s\t%s\n' % (contig, position, snpsite[checkID], ff, totalCount))
                        elif int(totalCount) < maxdepth and contig == 'Y':
                            if contig != pre_chr and pre_chr != '0':
                                chrfflstat = cal_basic_stats(tmpchrffl)
                                chrfflstat.insert(0, pre_chr)
                                chrffl.append(chrfflstat)
                                tmpchrffl = []
                            fetalsite = int(altCount)
                            ff = 2*fetalsite/mediancov
                            tmpchrffl.append(ff)
                            fo.write('%s\t%s\t%s\t%s\t%s\n' % (contig, position, snpsite[checkID], ff, totalCount))
                            pre_chr = contig
        
        chrfflstat = cal_basic_stats(tmpchrffl)
        chrfflstat.insert(0, pre_chr)
        chrffl.append(chrfflstat)                            
        totalFFest = 2*autofetaltotal/autoalltotal
        fflstat = cal_basic_stats(ffl)
        fflstat.insert(0, totalFFest)
        return fflstat, chrffl


def FFest_chrY(cfAERC, ffdir, mindepth, maxdepth, mediancov):
    ffdir = check_dirslash(ffdir)
    pref = 'chrY'
    ffhdout = '%s%s.ffest.out' % (ffdir, pref)
    chrffl = []

    mindepth = int(mindepth)
    maxdepth = int(maxdepth)

    with open(ffhdout, 'w+') as fo:
        fo.write("#CHROM\tPOS\tPH1\tFFest\tDEPTH\n")
        with open(cfAERC, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    if contig == 'Y':
                        ff = 2*int(altCount)/int(mediancov)
                        chrffl.append(ff)
                        fo.write('%s\t%s\t%s\t%s\t%s\n' % (contig, position, '.', ff, totalCount))

    chrfflstat = cal_basic_stats(chrffl)
    return chrfflstat
        
    
def FFest_main(HomoDiff, HeteroPat, cfAERC, workdir, mindepth, maxdepth, mediancov):
    workdir = check_dirslash(workdir)
    ffdir = workdir + '2_FFest'
    mkdir(ffdir)
    pref = basename(HomoDiff).rsplit('.', 1)[0]
    sumstats = ffdir + '/' + pref + '.ffest.summary'
    
    hdfflstat, hdchrffl = FFest_HomoDiffSNP(HomoDiff, cfAERC, ffdir, mindepth, maxdepth)
    dim0 = len(hdchrffl)
    chrYfflstat = FFest_chrY(cfAERC, ffdir, mindepth, maxdepth, mediancov)
    #disable heteropat SNP FF estimation
    #hpfflstat, hpchrffl = FFest_HeteroPatSNP(HeteroPat, cfAERC, ffdir, mindepth, maxdepth, mediancov)
    #dim1 = len(hpchrffl)
    with open(sumstats, 'w+') as fo:
        fo.write('##Fetal fraction estimation from homodiff SNPs on autosomes\n#overallFFest\tMean\tS.D.\tMedian\tMax\tMin\n')
        print >>fo, '\t'.join(map(str, hdfflstat))
        fo.write('#chromosome FFest.\n')
        for i in range(dim0):
            print >>fo, '\t'.join(map(str, hdchrffl[i]))
        #fo.write('##Fetal fraction estimation from heteropat SNPs on autosomes\n#overallFFest\tMean\tS.D.\tMedian\tMax\tMin\n')
        #print >>fo, '\t'.join(map(str, hpfflstat))
        #fo.write('#chromosome FFest.\n')
        #for i in range(dim1):
        #    print >>fo, '\t'.join(map(str, hpchrffl[i]))
        fo.write('##Fetal fraction estimation from chrY\nMean\tS.D.\tMedian\tMax\tMin\n')
        print >>fo, '\t'.join(map(str, chrYfflstat))
    
    
if __name__ == '__main__':
    FFest_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
