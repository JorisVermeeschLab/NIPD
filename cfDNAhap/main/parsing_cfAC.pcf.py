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

def Sum_FAR_Stats(fhtype, cffhapdir):
    pref = basename(fhtype)
    sumstats = '%s%s.summary' % (cffhapdir, pref)

    (tcnt, subA0, subA1, subB0, subB1, prepos, dist) = (0, 0, 0, 0, 0, 0, 0)
    (chrdistl, alldistl) = ([], [])
    distarr = np.empty(0)
    prechr = '0'

    with open(sumstats, 'w+') as fo:
        fo.write('##Summary for SNPs\n')
        with open(fhtype, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    (CHROM, POS, RSID, REF, ALT, PH1, PH2, MH1, MH2, FAR, TYPE, TCNT) = line.strip().split('\t')
                    tcnt += 1
                    if (TYPE == 'P1' and PH1 == '0') or (TYPE == 'M1' and MH1 == '0'):
                        subA0 += 1
                    elif (TYPE == 'P1' and PH1 == '1') or (TYPE == 'M1' and MH1 == '1'):
                        subA1 += 1
                    elif (TYPE == 'P2' and PH1 == '0') or (TYPE == 'M2' and MH1 == '0'):
                        subB0 += 1
                    elif (TYPE == 'P2' and PH1 == '1') or (TYPE == 'M2' and MH1 == '1'):
                        subB1 += 1
                    
                    if CHROM != prechr:
                        if prechr != '0':
                            distarr = np.append(distarr, prechr)
                            chrstat = cal_basic_stats(chrdistl)
                            distarr = np.append(distarr, chrstat)
                            chrdistl = []
                    elif CHROM == prechr:
                        dist = int(POS) - prepos
                        chrdistl.append(dist)
                        alldistl.append(dist)
                    prepos = int(POS)
                    prechr = CHROM
            #when loop ends, calculate the stats for the last chromosme
            distarr = np.append(distarr, prechr)
            chrstat = cal_basic_stats(chrdistl)
            dim1 = len(chrstat) + 1
            distarr = np.append(distarr, chrstat)

            dim0 = int(len(distarr)/dim1)
            distarr = np.reshape(distarr, (dim0, dim1))

            alldiststat = cal_basic_stats(alldistl)

            fo.write('#TotalNoSNPs\tP1_0/M1_0\tP1_1/M1_1\tP2_0/M2_0\tP2_1/M2_1\n')
            fo.write('%s\t%s\t%s\t%s\t%s\n' % (tcnt, subA0, subA1, subB0, subB1))
            fo.write('#Distance between SNPs\n#Mean\tS.D.\tMedian\tMax\tMin\noverall\t')
            print >>fo, '\t'.join(map(str, alldiststat))
            for i in range(dim0):
                print >>fo, '\t'.join(distarr[i])


def cal_heteropat_FAR(aerc, HeteroPat, mindepth, maxdepth, mediancov, cffhapdir):
    pref = basename(HeteroPat)
    farfo = '%s%s.aerc.tsv' % (cffhapdir, pref)
    
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    mediancov = int(mediancov)
    
    cffeval = dict()
    (PH1, PH2, MH1, MH2, TYPE) = ('9', '9', '9', '9', 'unknown')
    FAR = 9
    (prevline, prechrom, prevpos, dist, prevtype, prevcov) = ('na', '0', 0, 9999, 'na', 0)
    cnter = 0
    distlimit = 200
    
    with open(HeteroPat, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                val = '%s:%s:%s:%s' % (PH1, PH2, MH1, TYPE)
                cffeval[refID] = val
                
    with open(farfo, 'w+') as fo0:
        fo0.write('#CHROM\tPOS\tRSID\tREF\tALT\tPH1\tPH2\tMH1\tMH2\tFAR\tTYPE\tTCNT\n')
        with open(aerc, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    checkID = '%s:%s' % (contig, position)
                    if checkID in cffeval:
                        if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                            (PH1, PH2, MH1, TYPE) = cffeval[checkID].split(':')
                            MH2 = MH1
                            if contig != 'Y':
                                if MH1 == '0':
                                    FAR = 2*int(altCount)/int(totalCount)
                                elif MH1 == '1':
                                    FAR = 2*int(refCount)/int(totalCount)
                                    
                                #evaluate the distance of two same type SNPs; if the distance is smaller than distlimit, keep the one that has higher coverage
                                tmpline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, TYPE, totalCount)
                            
                                dist = int(position) - int(prevpos)
                                if contig == prechrom and dist < distlimit and TYPE == prevtype:
                                    if int(prevcov) > int(totalCount):
                                        pfline = prevline
                                        fo0.write(pfline)
                                else:
                                    if cnter == 0:
                                        tmpline = tmpline

                                    else:
                                        pfline = prevline
                                        fo0.write(pfline)
                            
                        elif contig == 'Y' and int(totalCount) <= maxdepth:
                            FAR = 2*int(altCount)/mediancov
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, TYPE, totalCount)
                            fo0.write(pfline)
                            
                        prevline = tmpline
                        prechrom = contig
                        prevpos = position
                        prevtype = TYPE
                        prevcov = totalCount
                        cnter += 1
                        
    Sum_FAR_Stats(farfo, cffhapdir)



def cal_heteromat_FAR(aerc, HeteroMat, mindepth, maxdepth, correctionft, cffhapdir):
    pref = basename(HeteroMat)
    farfo = '%s%s.aerc.tsv' % (cffhapdir, pref)

    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    correctionft = float(correctionft)
 
    adjref = 0.5/correctionft
    adjalt = 0.5/(1-correctionft)

    cffeval = dict()
    (PH1, PH2, MH1, MH2, TYPE) = ('9', '9', '9', '9', 'unknown')
    FAR = 9
    
    (prevline, prechrom, prevpos, dist, prevtype, prevcov) = ('na', '0', 0, 9999, 'na', 0)
    cnter = 0
    distlimit = 200
    
    with open(HeteroMat, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                val = '%s:%s:%s:%s' % (PH1, MH1, MH2, TYPE)
                cffeval[refID] = val

    with open(farfo, 'w+') as fo0:
        fo0.write('#CHROM\tPOS\tRSID\tREF\tALT\tPH1\tPH2\tMH1\tMH2\tFAR\tTYPE\tTCNT\n')
        with open(aerc, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('contig'):
                    (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                    checkID = '%s:%s' % (contig, position)
                    if checkID in cffeval:
                        if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                            (PH1, MH1, MH2, TYPE) = cffeval[checkID].split(':')
                            PH2 = PH1
                            
                            #need improvement, female and male fetus correction difference
                            if MH1 == '0':
                                FAR = (int(refCount)*adjref-int(altCount)*adjalt)/(int(refCount)*adjref+int(altCount)*adjalt) 
                            elif MH1 == '1':
                                FAR = (int(altCount)*adjalt-int(refCount)*adjref)/(int(refCount)*adjref+int(altCount)*adjalt)
                            
                            tmpline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, TYPE, totalCount)
                            
                            dist = int(position) - int(prevpos)
                            if contig == prechrom and dist < distlimit and TYPE == prevtype:
                                if int(prevcov) > int(totalCount):
                                    pfline = prevline
                                    fo0.write(pfline)
                            else:
                                if cnter == 0:
                                    tmpline = tmpline
                                else:
                                    pfline = prevline
                                    fo0.write(pfline)
                                    
                            prevline = tmpline
                            prechrom = contig
                            prevpos = position
                            prevtype = TYPE
                            prevcov = totalCount
                            cnter += 1
                            
    Sum_FAR_Stats(farfo, cffhapdir)



def Infer_Paternal_Inheritance(pat_fetalar, cffhapdir):
    pre_chr = '0'
    (pre_pos, spos, epos, valshift) = (-1, 0, 0, 0)
    (pre_homolog, homolog) = ('H0', 'H0')
    preType = 'P0'
    cutoff = 0.02
    pref=basename(pat_fetalar)
    pat_inherit = '%s%s.paternal_inheritance.tsv' % (cffhapdir, pref)
    multiarr = np.empty(0)
    
    with open(pat_inherit, 'w+') as fo:
        fo.write('chrom\tstart\tend\thomolog\tAR\tAlarm\n')
        with open(pat_fetalar, 'rU') as fh0:
            for line in fh0:
                if not line.strip().startswith('chr'):
                    (chrom, pos, ar, yhat, Type) = line.strip().split('\t')
                    #starting point of each chromosome
                    if chrom != pre_chr:
                        #when starting from other chromsome rather than chr1, specify end pos of previous chromsome
                        if pre_pos == -1:
                            pass
                        else:
                            if (preType == 'P2' and valshift > cutoff) or (preType == 'P1' and valshift < -cutoff):
                                epos = pre_pos
                                #pfline = '%s\t%s\t%s\t%s\t%s\n' % (pre_chr, spos, epos, pre_homolog, valshift)
                                #fo.write(pfline)
                                multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                        #initialize starting pos of each chromosome; 
                        #when P1 type value is ~FAR, homolog2 is inherited and when P2 type value is ~-FAR, homolog1 is inherited
                        if Type == 'P2' and float(yhat) > cutoff:
                            homolog = 'H1'
                            spos = int(pos)
                        elif Type == 'P1' and float(yhat) < -cutoff:
                            homolog = 'H2'
                            spos = int(pos)
                        else:
                            pass

                    #change point in a chromosome
                    elif chrom == pre_chr and pre_pos != -1:
                        if float(yhat) != valshift:
                            #within P1 or P2 type
                            if preType == Type:
                                #when the value changes from fetalAR to near 0, end of homolog
                                if (valshift > cutoff and float(yhat) < cutoff) or (valshift < -cutoff and float(yhat) > -cutoff):
                                    epos = pre_pos
                                    #pfline = '%s\t%s\t%s\t%s\t%s\n' % (chrom, spos, epos, pre_homolog, valshift)
                                    #fo.write(pfline)
                                    multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                                #when the value changes from near 0 to fetalAR, start of homolog
                                elif valshift < cutoff and float(yhat) > cutoff:
                                    homolog = 'H1'
                                    spos = int(pos)
                                elif valshift > -cutoff and float(yhat) < -cutoff:
                                    homolog = 'H2'
                                    spos = int(pos)
                                #when the value changes from one fetalAR to another fetalAR; remain the same homolog, both start and end happens
                                elif (valshift > cutoff and float(yhat) > cutoff) or (valshift < -cutoff and float(yhat) < -cutoff):
                                    epos = pre_pos
                                    multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                                    spos = int(pos)
 
                            #when P1 change to P2 or P2 to P1        
                            else:
                                #both start and end happens
                                if (valshift > cutoff and float(yhat) < -cutoff) or (valshift < -cutoff and float(yhat) > cutoff):
                                    epos = pre_pos
                                    #pfline = '%s\t%s\t%s\t%s\t%s\n' % (chrom, spos, epos, pre_homolog, valshift)
                                    #fo.write(pfline)
                                    multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                                    if float(yhat) > cutoff:
                                        homolog = 'H1'
                                        spos = int(pos)
                                    elif float(yhat) < -cutoff:
                                        homolog = 'H2'
                                        spos = int(pos)
                                #end of homolog
                                elif (valshift > cutoff and float(yhat) > -cutoff) or (valshift < -cutoff and float(yhat) < cutoff):
                                    epos = pre_pos
                                    #pfline = '%s\t%s\t%s\t%s\t%s\n' % (chrom, spos, epos, pre_homolog, valshift)
                                    #fo.write(pfline)
                                    multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                                #start of homolog
                                elif valshift > -cutoff and float(yhat) > cutoff:
                                    homolog = 'H1'
                                    spos = int(pos)
                                elif valshift < cutoff and float(yhat) < -cutoff:
                                    homolog = 'H2'
                                    spos = int(pos)
                                
                    valshift = float(yhat)
                    pre_pos = int(pos)
                    pre_chr = chrom
                    preType = Type
                    pre_homolog = homolog

            dim1 = int(len(multiarr)/5)
            multiarr = np.reshape(multiarr, (dim1, 5))
            #multiarr[:,1] = multiarr[:,1].astype(int)
            #sort array by chromosome and start position
            for i in range(0,22):
                temparr = multiarr[multiarr[:,0]==str(i+1)]
                if i == 0:
                    sortarrcur = temparr[temparr[:,1].astype(int).argsort()]
                    sortarr = sortarrcur
                                      
                else:
                    sortarrcur = temparr[temparr[:,1].astype(int).argsort()]
                    sortarr = np.append(sortarr, sortarrcur)
                                       
            sortarr = np.reshape(sortarr, (dim1, 5))
            #mark conflicts (overlap) in breakpoint change
            #need more codes here to resolve conflicts
            prechr = '0'
            for j in range(0,dim1):
                #start of chromosome
                if sortarr[j, 0] != prechr:
                    print >>fo, '\t'.join(np.append(sortarr[j], 'pass'))
                
                #same chromosome with no conflict
                elif sortarr[j, 0] == prechr and sortarr[j,1] >= sortarr[j-1,2]:
                    print >>fo, '\t'.join(np.append(sortarr[j], 'pass'))
                
                #same chromosome with conflict
                elif sortarr[j, 0] == prechr and sortarr[j,1] < sortarr[j-1,2]:
                    print >>fo, '\t'.join(np.append(sortarr[j], 'mark'))
                prechr = sortarr[j, 0]



def Update_Hetboth_Type(paternal_inheritance, phased_hetboth, cffhapdir):
    multiarr = np.empty(0)
    pref = basename(phased_hetboth)
    hetbothorder = '%s%s.paternal_infer.hetboth.order.tsv' % (cffhapdir, pref)
    Type = 'M0'

    with open(paternal_inheritance, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('chrom'):
                (chrom, spos, epos, homolog, ar, alarm) = line.strip().split('\t')
                multiarr = np.append(multiarr, [[chrom, spos, epos, homolog, ar, alarm]])
    dim1 = int(len(multiarr)/6)
    multiarr = np.reshape(multiarr, (dim1, 6))
    
    with open(hetbothorder, 'w+') as fo:
        fo.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\n')
        with open(phased_hetboth, 'rU') as fh1:
            i = 0
            refchr = multiarr[i,0]
            refstart = multiarr[i,1].astype(int)
            refend = multiarr[i,2].astype(int)
            refhomolog = multiarr[i,3]
            for line in fh1:
                if not line.strip().startswith('#'):
                    (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE) = line.strip().split('\t')
                    
                    #SNP resides in identified haplotype block                
                    if CHROM == refchr and int(POS) >= refstart and int(POS) <= refend:
                        if refhomolog == 'H1':
                            Type = 'M1'
                        elif refhomolog == 'H2':
                            Type = 'M2'
                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, Type)
                        fo.write(pfline)
                    elif CHROM == refchr and int(POS) > refend:
                        if i < dim1-1:
                            i += 1
                            refchr = multiarr[i,0]
                            refstart = multiarr[i,1].astype(int)
                            refend = multiarr[i,2].astype(int)
                            refhomolog = multiarr[i,3]
                            if CHROM == refchr and int(POS) >= refstart and int(POS) <= refend:
                                if refhomolog == 'H1':
                                    Type = 'M1'
                                elif refhomolog == 'H2':
                                    Type = 'M2'
                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, Type)
                                fo.write(pfline)
                    elif CHROM != refchr:
                        if int(CHROM) < int(refchr):
                            pass
                        else:
                            if i < dim1-1:
                                i += 1
                                refchr = multiarr[i,0]
                                refstart = multiarr[i,1].astype(int)
                                refend = multiarr[i,2].astype(int)
                                refhomolog = multiarr[i,3]
                                if CHROM == refchr and int(POS) >= refstart and int(POS) <= refend:
                                    if refhomolog == 'H1':
                                        Type = 'M1'
                                    elif refhomolog == 'H2':
                                        Type = 'M2'
                                    pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, Type)
                                    fo.write(pfline)


def combine_hetmat_hetboth(hetmataerc, hetbothaerc, pref, cffhapdir):
    multiarr = np.empty(0)
    hetmatboth = '%s%s.heteromatboth.aerc.tsv' % (cffhapdir, pref)

    with open(hetmataerc, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, PH1, PH2, MH1, MH2, FAR, TYPE, TCNT) = line.strip().split('\t')
                if CHROM == 'X':
                    CHROM = '23'
                multiarr = np.append(multiarr, [[CHROM, POS, FAR, TYPE]])
                       
    with open(hetbothaerc, 'rU') as fh1:
        for line in fh1:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, PH1, PH2, MH1, MH2, FAR, TYPE, TCNT) = line.strip().split('\t')
                if CHROM == 'X':
                    CHROM = '23'
                multiarr = np.append(multiarr, [[CHROM, POS, FAR, TYPE]])
    
    dim1 = int(len(multiarr)/4)
    multiarr = np.reshape(multiarr, (dim1, 4))
    #sort chromosome order first
    sortchrarr = multiarr[multiarr[:,0].astype(int).argsort()]
    #sort position in each chromosome
    for i in range(0,24):
        sortposarr = sortchrarr[sortchrarr[:,0]==str(i)]
        if i == 0:
            sortposarr = sortposarr[sortposarr[:,1].astype(int).argsort()]
            sortarr = sortposarr
        else:
            sortposarr = sortposarr[sortposarr[:,1].astype(int).argsort()]
            sortarr = np.append(sortarr, sortposarr)

    sortarr = np.reshape(sortarr, (dim1, 4))

    with open(hetmatboth, 'w+') as fo0:
        fo0.write('#CHROM\tPOS\tFAR\tTYPE\n')
        for j in range(0, dim1):
            if sortarr[j, 0] == '23':
                sortarr[j, 0] = 'X'
            print >>fo0, '\t'.join(sortarr[j])




def get_tasktype(tasktype):
    if tasktype.lower() == 'heteropat':
        tasktype = 'heteropat'
    elif tasktype.lower() == 'patinherit':
        tasktype = 'patinherit'
    elif tasktype.lower() == 'heteromat':
        tasktype = 'heteromat'
    elif tasktype.lower() == 'heteromatboth':
        tasktype = 'heteromatboth'
    elif tasktype.lower() == 'heteropatboth':
        tasktype = 'heteropatboth'
    else:
        raise ValueError('File not exists.\n')
    return(tasktype)
        



def calAR_main(aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir):
    try:
        tasktype = get_tasktype(tasktype)
    except ValueError:
        print 'SNP file not exists. please check filename.'

    if tasktype == 'heteropat':
        infh = parentshapdir + pref + '.heteropat.tsv'
        cal_heteropat_FAR(aerc, infh, mindepth, maxdepth, mediancov, cffhapdir)
    elif tasktype == 'patinherit':    
        bpinfh = cffhapdir + pref + '.heteropat.aerc.bp.tsv'
        Infer_Paternal_Inheritance(bpinfh, cffhapdir)
        phased_hetboth = parentshapdir + pref + '.heteroboth.tsv'
        pat_inherit = '%s%s.heteropat.aerc.bp.paternal_inheritance.tsv' % (cffhapdir, pref)
        Update_Hetboth_Type(pat_inherit, phased_hetboth, cffhapdir)
    elif tasktype == 'heteromat':
        infh = parentshapdir + pref + '.heteromat.tsv'
        cal_heteromat_FAR(aerc, infh, mindepth, maxdepth, correctionft, cffhapdir)
    elif tasktype == 'heteromatboth':
        heterobothorder = '%s%s.heteroboth.paternal_infer.hetboth.order.tsv' % (cffhapdir, pref)
        cal_heteromat_FAR(aerc, heterobothorder, mindepth, maxdepth, correctionft, cffhapdir)
        hetmataerc = '%s%s.heteromat.aerc.tsv' % (cffhapdir, pref)
        hetbothorderaerc = '%s%s.heteroboth.paternal_infer.hetboth.order.aerc.tsv' % (cffhapdir, pref)
        combine_hetmat_hetboth(hetmataerc, hetbothorderaerc, pref, cffhapdir)
    elif tasktype == 'heteropatboth':
        infh = cffhapdir + pref + '.heteropatboth.tsv'



if __name__ == '__main__':
    calAR_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
