#!/usr/bin/env python

'''
Author: Huiwen Che
Created: Mar 05 2018
Merging scripts from parsingAERC_heteropat.py, parsingAERC_heteromat.py, parsingAERC_hetboth.py; modifying input and output handling according to modification in previous steps
Use: 
parsing cfDNA pileup (gatk aerc count file) according to different SNPs types
'''

from __future__ import division
import sys
import os
from Helper import mkdir, check_dirslash, basename
from Utility_stats import cal_basic_stats
import numpy as np

                                
def cffhap_main(aerc, pref, mindepth, maxdepth, mediancov, correctionft, workdir, script_calAR, script_findbp):
    workdir = check_dirslash(workdir)
    cffhapdir = workdir + '4_cffhap'
    mkdir(cffhapdir)
    cffhapdir = check_dirslash(cffhapdir)
    cffhapscript = cffhapdir + 'cffhap.sh'
    parentshapdir = workdir + '1_parentshap'
    parentshapdir = check_dirslash(parentshapdir)
    #PCF parameters
    ratio = 'FAR'
    winsize = 20
    kmin = 5
    gamma1 = 4
    gamma2 = 4
    penalty = 'T'
    
    with open(cffhapscript, 'w+') as fo:
        fo.write('#!/bin/bash\n\n')
        tasktype = 'heteropat'
        fo.write('##processing heteropat type2 SNPs\n')
        heteropat_ar_task = 'python %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir)
        fo.write(heteropat_ar_task)
        heteropat_findbreak_task = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s %s %s %s %s\n\n' % (script_findbp, tasktype, pref, ratio, winsize, kmin, gamma1, gamma2, penalty, cffhapdir)
        fo.write(heteropat_findbreak_task)

        tasktype = 'patinherit'
        fo.write('##inferring paternal allele inheritance based on type2 SNPs breakpoint detection\n')
        pat_inheritance_task = 'python %s %s %s %s %s %s %s %s %s %s\n\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir)
        fo.write(pat_inheritance_task)

        tasktype = 'heteromat'
        fo.write('##processing heteromat type3 SNPs\n')
        heteromat_ar_task = 'python %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir)
        fo.write(heteromat_ar_task)
        heteromat_findbreak_task = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s %s %s %s %s\n\n' % (script_findbp, tasktype, pref, ratio, winsize, kmin, gamma1, gamma2, penalty, cffhapdir)
        fo.write(heteromat_findbreak_task)

        tasktype = 'heteromatboth'
        fo.write('##processing heteroboth type4 SNPs that have been categorized using paternal haplotyping information\n')
        heteroboth_ar_task = 'python %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir)
        fo.write(heteroboth_ar_task)
        heteromatboth_findbreak_task = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s %s %s %s %s\n\n' % (script_findbp, tasktype, pref, ratio, winsize, kmin, gamma1, gamma2, penalty, cffhapdir)
        fo.write(heteromatboth_findbreak_task)

    
if __name__ == '__main__':
    cffhap_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
