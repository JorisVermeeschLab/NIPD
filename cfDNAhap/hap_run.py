#!/usr/bin/env python

import sys
import os
from subprocess import call
from main.Helper import check_dirslash

def getped(pedfh):
    try:
        with open(pedfh) as fh:
            lines = fh.readlines()
            infoline = lines[2]
    except IOError:
        print 'Fail to open pedfile'
        
    return (infoline)
                


def preproc(hapdir, workdir, scriptdir, familyID, cfID, ffsex, ffmedcov):
    #workdir = check_dirslash(workdir)
    #hapdir = workdir + '5_variant'
    #hapdir = check_dirslash(hapdir)
    #scriptdir = check_dirslash(scriptdir)

    configfh = workdir + 'ConfigFile'
    shcmd1 = 'cp %s %s' % (configfh, hapdir)
    call(shcmd1, shell=True)

    pedfh = '%s%s.ped' % (workdir, familyID)
    infoline = getped(pedfh)
    (familyID, sibID, patID, matID, sibsex, phenotype) = infoline.strip().split()

    if (sibsex == '1' or sibsex == 1):
        sibsex = 'male'
    else:
        sibsex = 'female'

    fragsisfo = hapdir + 'fragcov.sis.list'
    with open(fragsisfo, 'w+') as fo1:
        fo1.write('filedir\tsamplename\tsex\n')
        patcovsum = '%s4_cov/%s/%s.sorted.merged.md.filtered.fragcov.sample_interval_summary\t%s\tmale\n' % (workdir, patID, patID, patID)
        matcovsum = '%s4_cov/%s/%s.sorted.merged.md.filtered.fragcov.sample_interval_summary\t%s\tfemale\n' % (workdir, matID, matID, matID)
        sibcovsum = '%s4_cov/%s/%s.sorted.merged.md.filtered.fragcov.sample_interval_summary\t%s\t%s\n' % (workdir, sibID, sibID, sibID, sibsex)
        #female is used for cfDNA coverage calculation; may change
        cfcovsum = '%s4_cov/%s/%s.sorted.merged.md.filtered.fragcov.sample_interval_summary\t%s\tfemale\n' % (workdir, cfID, cfID, cfID)
        fo1.write(patcovsum)
        fo1.write(matcovsum)
        fo1.write(sibcovsum)
        fo1.write(cfcovsum)
    
    ffmedcov = int(ffmedcov)
    covupper = 500
    covlower = 10
    cfcovupper = 500
    cfcovlower = 30
    
    cfaerc = '%s6_AERC/%s/%s.sorted.merged.md.filtered.bam.frag.aerc.2x.csv' % (workdir, cfID, cfID)
    script_parph = scriptdir + 'main/parentalvcf_phasing.py'
    script_seqerr = scriptdir + 'main/seqerr.py'
    script_sibraf = scriptdir + 'main/sib_raf_pcf.R'
    script_ffest = scriptdir + 'main/ffest.py'
    script_cflen = scriptdir + 'main/cflen.py'
    script_excfreads = scriptdir + 'main/extract_cffreads.py'
    script_plotlen = scriptdir + 'fig/plot_cff_lendist.R'

    
    prehap_fo = hapdir + 'PreHapAnalysisRun.sh'
    with open(prehap_fo, 'w+') as fo:
        fo.write('#!/bin/bash\n\n')

        fo.write('#Parental data phasing\n')
        parvcf = '%s0_rawvcf/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.vcf' % (hapdir, familyID)
        task1 = '/home/hche0/miniconda2/bin/python %s %s %s %s %s %s %s\n\n' % (script_parph, parvcf, pedfh, covlower, covupper, hapdir, ffsex)
        fo.write(task1)

        fo.write('#error rate estimation\n')
        homosamefh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.homosame.tsv' % (hapdir, familyID)
        task2 = '/home/hche0/miniconda2/bin/python %s %s %s %s %s %s\n\n' % (script_seqerr, homosamefh, cfaerc, cfcovlower, cfcovupper, hapdir)
        fo.write(task2)

        fo.write('#Parental-sibling data examination\n')
        hetpatfh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteropat.tsv' % (hapdir, familyID)
        hetmatfh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteromat.tsv' % (hapdir, familyID)
        task3 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s\n\n' % (script_sibraf, hetpatfh, hetmatfh, sibsex)
        fo.write(task3)

        fo.write('#Fetal fraction estimation\n')
        homodifffh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.homodiff.tsv' % (hapdir, familyID)
        task4 = '/home/hche0/miniconda2/bin/python %s %s %s %s %s %s %s %s\n\n' % (script_ffest, homodifffh, hetpatfh, cfaerc, hapdir, cfcovlower, cfcovupper, ffmedcov)
        fo.write(task4)

        fo.write('#extract homodiff fetal vs. shared fragments lengths\n')
        cfbam = '%s3_align/%s/%s.sorted.merged.md.filtered.bam' % (workdir, cfID, cfID)
        rawvcf = '%s0_rawvcf/%s_parental_sib.joint.raw.snps.indels.vcf' % (hapdir, familyID)
        task5 = '/home/hche0/miniconda2/bin/python %s %s %s %s %s %s %s\n\n' % (script_cflen, hapdir, homodifffh, cfbam, rawvcf, script_excfreads, script_plotlen)
        fo.write(task5)

    return(prehap_fo)


def runPreCompute(hapdir, workdir, scriptdir, familyID, cfID, ffsex, ffmedcov):

    prehap_fo = preproc(hapdir, workdir, scriptdir, familyID, cfID, ffsex, ffmedcov)
    shcmd2 = 'sh %s' % (prehap_fo)
    call(shcmd2, shell=True)

    ffestdir = hapdir + '2_FFest'
    for fh in os.listdir(ffestdir):
        if fh.endswith('.ffest.summary'):
            ffestfh = ffestdir + '/' + fh

    if os.path.isfile(ffestfh):
        with open(ffestfh, 'rU') as fh1:
            lines = fh1.readlines()
            tarline = lines[2]
            (overallff, meanff, sdff, medff, maxff, minff) = tarline.strip().split('\t')
            ffest = '{0:.5f}'.format(float(medff))
                
    else:
        raise ValueError('ffest.summary file not found')

    return(ffest)


def mainhap(hapdir, scriptdir, familyID, cfID, ffsex, ffmedcov):
    hapdir = check_dirslash(hapdir)
    workdir = hapdir.rsplit('/', 2)[0]
    workdir = check_dirslash(workdir)
    #hapdir = workdir + '5_variant'
    #hapdir = check_dirslash(hapdir)
    scriptdir = check_dirslash(scriptdir)

    pedfh = '%s%s.ped' % (workdir, familyID)
    infoline = getped(pedfh)
    (familyID, sibID, patID, matID, sibsex, phenotype) = infoline.strip().split()

    ffest = runPreCompute(hapdir, workdir, scriptdir, familyID, cfID, ffsex, ffmedcov)

    script_cbshap = scriptdir + 'main/cffhap.CBS.py'
    script_cnveval = scriptdir + 'main/cnv_eval.R'
    script_plotgtr = scriptdir + 'fig/plot_gtrellis.R'
    script_plotchr = scriptdir + 'fig/plot_fetus_chr.R'
    script_plotparsib = scriptdir + 'fig/plot_parentsib_phase.R'

    ffmedcov = int(ffmedcov)
    covupper = 500
    covlower = 10
    cfcovupper = 500
    cfcovlower = 30
    pref = '%s_parental_sib.joint.raw.snps.indels.phasedbytrans.%s_%sx' % (familyID, covlower, covupper)
    #reference allele bias correction; need to improve the correction method                                                                                                                        
    rafcorfac = 0.525
    
    cfaerc = '%s6_AERC/%s/%s.sorted.merged.md.filtered.bam.frag.aerc.2x.csv' % (workdir, cfID, cfID)
    fragsisfo = hapdir + 'fragcov.sis.list'

    mainhap_fo = hapdir + 'HapAnalysisRun.sh'
    with open(mainhap_fo, 'w+') as fo:
        fo.write('#!/bin/bash\n\n')

        fo.write('#main haplotyping\n')
        task6 = '/home/hche0/miniconda2/bin/python %s %s %s %s %s %s %s %s %s %s %s %s\n\n' % (script_cbshap, cfaerc, pref, familyID, cfcovlower, cfcovupper, ffmedcov, rafcorfac, hapdir, scriptdir, ffest, ffsex)
        fo.write(task6)

        fo.write('#cnv\n')
        homodiffestfh = '%s2_FFest/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.%s_%sx.homodiff.ffest.out' % (hapdir, familyID, covlower, covupper)
        task7 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s /uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/0_data/siChild_SNPsCoveredChrRmv.mergeinterval.gc_profile.bed %s 20 5 30 4 T %s totalcov\n\n' % (script_cnveval, fragsisfo, homodiffestfh, hapdir)
        fo.write(task7)

        fo.write('#plot\n')
        patcbs = '%s4_cffhap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.%s_%sx.heteropat.aerc.CBSseg.tsv' % (hapdir, familyID, covlower, covupper)
        matcbs = '%s4_cffhap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.%s_%sx.heteromat.aerc.CBSseg.tsv' % (hapdir, familyID, covlower, covupper)
        matbothcbs = '%s4_cffhap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.%s_%sx.heteromatboth.aerc.CBSseg.tsv' % (hapdir, familyID, covlower, covupper)
        homodiffcnv = '%s5_cnv/homodiff.ffest.tsv' % (hapdir)
        cfcnv = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, cfID)
        task8 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s.%s_%sx.hetmat %s png %s\n' % (script_plotgtr, patcbs, matcbs, homodiffcnv, cfcnv, cfID, cfcovlower, cfcovupper, hapdir, ffest)
        task9 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s.%s_%sx.hetmatboth %s png %s\n\n' % (script_plotgtr, patcbs, matbothcbs, homodiffcnv, cfcnv, cfID, cfcovlower, cfcovupper, hapdir, ffest)
        fo.write(task8)
        fo.write(task9)

        patcnv = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, patID)
        matcnv = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, matID)
        figdir = '%s6_fig/' % (hapdir)
        figscriptdir = '%sfig' % (scriptdir)
        task10 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s %s %s %s %s %s %s\n\n' % (script_plotchr, patcbs, matbothcbs, cfcnv, homodiffcnv, patcnv, matcnv, familyID, figdir, ffest, figscriptdir, ffsex)
        fo.write(task10)

        hetpatfh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteropat.tsv' % (hapdir, familyID)
        hetmatfh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteromat.tsv' % (hapdir, familyID)
        hetbothfh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteroboth.tsv' % (hapdir, familyID)
        homodifffh = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.homodiff.tsv' % (hapdir, familyID)
        hetpatsibbp = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteropat.sib.bp.tsv' % (hapdir, familyID)
        hetmatsibbp = '%s1_parentshap/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.10_500x.heteromat.sib.bp.tsv' % (hapdir, familyID)
        patcov = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, patID)
        matcov = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, matID)
        sibcov = '%s5_cnv/%s.totalcov.gccor.noref.logR.tsv' % (hapdir, sibID)
        task11 = '/cm/shared/apps/R/3.2.4/bin/Rscript %s %s %s %s %s %s %s %s %s %s %s %s\n\n' % (script_plotparsib, hetpatfh, hetmatfh, hetbothfh, homodifffh, hetpatsibbp, hetmatsibbp, patcov, matcov, sibcov, familyID, figdir)
        fo.write(task11)
        

if __name__ == '__main__':
    mainhap(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
