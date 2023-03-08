#!/usr/bin/env python

'''
use: extract fetal-specific reads length and shared reads length
'''

from __future__ import division
import sys
import os
import Helper
from Utility_stats import cal_basic_stats
import numpy as np
from subprocess import call


(fastqdir, workdir, samtools, bwa, java, picard, fastqc, gatk, Rscript, reference, target, targetflank50, targetflank100, targetSNPref, dbsnpref, pedfile) = Helper.get_config()

def make_homodiffbed(HomoDiff, piledir):
    pref = Helper.basename(HomoDiff)
    homodiffbed = '%s%s.bed' % (piledir, pref)
    
    #create homodiff reference sites bed file
    with open(homodiffbed, 'w+') as fo0:
        with open(HomoDiff, 'rU') as fh0:
            for line in fh0:
                if not line.strip().startswith('#'):
                    L0 = line.strip().split('\t')
                    (CHROM, POS, RSID) = (L0[0], L0[1], L0[2])
                    fo0.write('%s\t%s\t%s\t%s\n' % (CHROM, int(POS)-1, POS, RSID))
                    
    return(homodiffbed)
    
   
def fflen_main(worksubdir, HomoDiff, cfbam, rawvcf, script_extreads, script_plotlen):
    worksubdir = Helper.check_dirslash(worksubdir)
    lendir = worksubdir + '3_cffLEN'
    Helper.mkdir(lendir)
    lendir = Helper.check_dirslash(lendir)
    piledir = lendir + 'gpile'
    Helper.mkdir(piledir)
    piledir = Helper.check_dirslash(piledir)
    
    homodiffbed = make_homodiffbed(HomoDiff, piledir)
    
    pref = Helper.basename(HomoDiff)
    gpilescript = '%sgpile.sh' % piledir
    extract_cffreads = script_extreads
    plotdist = script_plotlen
    pileout = '%s%s.gpile.out' % (piledir, pref)
    fetalreadsbam = '%s%s.fetalreads.bam' % (piledir, pref)
    sharereadsbam = '%s%s.sharedreads.bam' % (piledir, pref)
    
    pref1 = Helper.basename(pileout)
    fetalreads = '%s%s.fetalreads.txt' % (piledir, pref1)
    sharedreads = '%s%s.sharedreads.txt' % (piledir, pref1)
    fetallen = '%s%s.fetalreads.length.out' % (lendir, pref1)
    sharedlen = '%s%s.sharedreads.length.out' % (lendir, pref1)
                 
    with open(gpilescript, 'w+') as fo0:
        fo0.write('#!/bin/bash\n\n')
        gplie_task = '%s -Xmx4G -jar %s -T Pileup -R %s -I %s -L %s -o %s --metadata %s --showVerbose\n\n' % (java, gatk, reference, cfbam, homodiffbed, pileout, rawvcf)
        fo0.write(gplie_task)
        readsname_task = 'python %s %s %s\n\n' % (extract_cffreads, pileout, piledir)
        fo0.write(readsname_task)
        extfetalreads_task = '%s -Xmx4G -jar %s FilterSamReads I=%s O=%s READ_LIST_FILE=%s FILTER=includeReadList\n' % (java, picard, cfbam, fetalreadsbam, fetalreads)
        extsharereads_task = '%s -Xmx4G -jar %s FilterSamReads I=%s O=%s READ_LIST_FILE=%s FILTER=includeReadList\n' % (java, picard, cfbam, sharereadsbam, sharedreads)
        fo0.write(extfetalreads_task)
        fo0.write(extsharereads_task)
        extfetallen_task = '%s view %s | awk -F\'\\t\' \'{printf("%%s\\t%%s\\n"), $1, sqrt($9*$9)}\' | sort -u > %s\n' % (samtools, fetalreadsbam, fetallen)
        extsharelen_task = '%s view %s | awk -F\'\\t\' \'{printf("%%s\\t%%s\\n"), $1, sqrt($9*$9)}\' | sort -u > %s\n\n' % (samtools, sharereadsbam, sharedlen)
        fo0.write(extfetallen_task)
        fo0.write(extsharelen_task)

        plotdist_task = '%s %s %s %s %s %s\n' % (Rscript, plotdist, fetallen, sharedlen, 0, 500)
        fo0.write(plotdist_task)
    

if __name__ == '__main__':
    fflen_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
