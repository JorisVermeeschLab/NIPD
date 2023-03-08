#!/usr/bin/env python

'''
Author: Huiwen Che
Created: May 15 2018
Use: 
parsing cfDNA pileup (gatk aerc count file) according to different SNPs types; perform segmentation using CBS
'''

from __future__ import division
import sys
import os
from Helper import mkdir, check_dirslash, basename
from Utility_stats import cal_basic_stats
import numpy as np

                                
def cffhap_main_CBS(aerc, pref, sampleid, mindepth, maxdepth, mediancov, correctionft, workdir, scriptdir, ffest, ffsex):
    workdir = check_dirslash(workdir)
    cffhapdir = workdir + '4_cffhap'
    mkdir(cffhapdir)
    cffhapdir = check_dirslash(cffhapdir)
    cffhapscript = cffhapdir + 'cffhap.sh'
    parentshapdir = workdir + '1_parentshap'
    parentshapdir = check_dirslash(parentshapdir)
    scriptdir = check_dirslash(scriptdir)
    script_calAR = scriptdir + 'main/parsing_cfAC.CBS.py'
    script_callFetalseg = scriptdir + 'seg/callFetalseg.R'
    script_CBS = scriptdir + 'seg/CBS.R'
    script_callresHL = scriptdir + 'seg/callresHL.R'
    script_SegResolve = scriptdir + 'seg/SegResolve.R'
    script_mergeT34 = scriptdir + 'seg/mergeType34.R'
    script_switchres = scriptdir + 'diagnose/switchResolution.py'
    
    with open(cffhapscript, 'w+') as fo:
        fo.write('#!/bin/bash\n\n')
        tasktype = 'heteropat'
        fo.write('#Heteropat SNPs\n')
        fo.write('##processing heteropat type2 SNPs\n')
        heteropat_ar_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir, ffest)
        fo.write(heteropat_ar_task)
        
        fo.write('##CBS segmentation\n')
        heteropat_ar = '%s%s.heteropat.aerc.tsv' % (cffhapdir, pref)
        heteropat_findbreak_task = 'Rscript %s %s %s %s %s %s\n' % (script_callFetalseg, script_CBS, heteropat_ar, sampleid, ffest, ffsex)
        fo.write(heteropat_findbreak_task)
        
        fo.write('##Resolve CBS segmentation\n')
        heteropat_cbs = '%s%s.heteropat.aerc.CBSseg.tsv' % (cffhapdir, pref)
        heteropat_cbssum = '%s%s.heteropat.aerc.CBSseg.sum.tsv' % (cffhapdir, pref)
        resolve_break_task = 'Rscript %s %s %s %s %s\n\n' % (script_callresHL, script_SegResolve, heteropat_cbs, heteropat_cbssum, ffsex)
        fo.write(resolve_break_task)
        

        tasktype = 'heteromat1'
        fo.write('#Heteromat SNPs\n')
        fo.write('##processing heteromat type3 SNPs\n')
        heteromat_ar_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir, ffest)
        fo.write(heteromat_ar_task)
        
        fo.write('##CBS segmentation\n')
        heteromat_ar = '%s%s.heteromat.aerc.tsv' % (cffhapdir, pref)
        heteromat_findbreak_task = 'Rscript %s %s %s %s %s %s\n' % (script_callFetalseg, script_CBS, heteromat_ar, sampleid, ffest, ffsex)
        fo.write(heteromat_findbreak_task)
        
        fo.write('##Resolve CBS segmentation\n')
        heteromat_cbs = '%s%s.heteromat.aerc.CBSseg.tsv' % (cffhapdir, pref)
        heteromat_cbssum = '%s%s.heteromat.aerc.CBSseg.sum.tsv' % (cffhapdir, pref)
        resolve_break_task = 'Rscript %s %s %s %s %s\n\n' % (script_callresHL, script_SegResolve, heteromat_cbs, heteromat_cbssum, ffsex)
        fo.write(resolve_break_task)

        tasktype = 'heteromatboth'
        fo.write('#Heteromatboth SNPs\n')
        fo.write('##merge heteromat and heteroboth snps\n')
        heteropat_segresolve = '%s%s.heteropat.aerc.CBSseg.resolved.sum.tsv' % (cffhapdir, pref)
        heteromatraw = '%s%s.heteromat.tsv' % (parentshapdir, pref)
        heterobothraw = '%s%s.heteroboth.tsv' % (parentshapdir, pref)
        merge_heteromat_heteroboth_task = 'Rscript %s %s %s %s %s\n' % (script_mergeT34, script_SegResolve, heteropat_segresolve, heteromatraw, heterobothraw)
        fo.write(merge_heteromat_heteroboth_task)
                
        fo.write('##processing heteroboth type4 SNPs that have been categorized using paternal haplotyping information\n')
        heteroboth_ar_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s %s %s %s %s %s %s %s %s %s\n' % (script_calAR, aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir, ffest)
        fo.write(heteroboth_ar_task)
        
        fo.write('##CBS segmentation\n')
        heteromatboth_ar = '%s%s.heteromatboth.aerc.tsv' % (cffhapdir, pref)
        heteromatboth_findbreak_task = 'Rscript %s %s %s %s %s %s\n' % (script_callFetalseg, script_CBS, heteromatboth_ar, sampleid, ffest, ffsex)
        fo.write(heteromatboth_findbreak_task)
        
        fo.write('##Resolve CBS segmentation\n')
        heteromatboth_cbs = '%s%s.heteromatboth.aerc.CBSseg.tsv' % (cffhapdir, pref)
        heteromatboth_cbssum = '%s%s.heteromatboth.aerc.CBSseg.sum.tsv' % (cffhapdir, pref)
        resolve_break_task = 'Rscript %s %s %s %s %s\n\n' % (script_callresHL, script_SegResolve, heteromatboth_cbs, heteromatboth_cbssum, ffsex)
        fo.write(resolve_break_task)

        fo.write('##Switchover resolutions\n')
        switchreshetpat_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s\n' % (script_switchres, heteropat_segresolve)
        heteromat_segresolve = '%s%s.heteromat.aerc.CBSseg.resolved.sum.tsv' % (cffhapdir, pref)
        switchreshetmat_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s\n' % (script_switchres, heteromat_segresolve)
        heteromatboth_segresolve = '%s%s.heteromatboth.aerc.CBSseg.resolved.sum.tsv' % (cffhapdir, pref)
        switchreshetmatboth_task = '/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/anaconda2/bin/python %s %s\n' % (script_switchres, heteromatboth_segresolve)
        fo.write(switchreshetpat_task)
        fo.write(switchreshetmat_task)
        fo.write(switchreshetmatboth_task)
        

    
if __name__ == '__main__':
    cffhap_main_CBS(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
