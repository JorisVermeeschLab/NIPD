#!/usr/bin/env python

'''
Use: 
based on calulated paternal haplotype (using paternal inherited alleles), deducing the inherited paternal homolog in fetal;
infer allele inherited in both parental heterozygous sites
'''

from __future__ import division
import sys
import os
import numpy as np

def Infer_Paternal_Inheritance(pat_fetalar, workdir):
    pre_chr = '0'
    (pre_pos, spos, epos, valshift) = (-1, 0, 0, 0)
    (pre_homolog, homolog) = ('H0', 'H0')
    preType = 'P0'
    cutoff = 0.02
    pat_inherit = '%spaternal_inheritance.tsv' % (workdir)
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
                            if (preType == 'P1' and valshift > cutoff) or (preType == 'P2' and valshift < -cutoff):
                                epos = pre_pos
                                #pfline = '%s\t%s\t%s\t%s\t%s\n' % (pre_chr, spos, epos, pre_homolog, valshift)
                                #fo.write(pfline)
                                multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                        #initialize starting pos of each chromosome; 
                        #when P1 type value is ~FAR, homolog2 is inherited and when P2 type value is ~-FAR, homolog1 is inherited
                        if Type == 'P1' and float(yhat) > cutoff:
                            homolog = 'H2'
                            spos = int(pos)
                        elif Type == 'P2' and float(yhat) < -cutoff:
                            homolog = 'H1'
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
                                    homolog = 'H2'
                                    spos = int(pos)
                                elif valshift > -cutoff and float(yhat) < -cutoff:
                                    homolog = 'H1'
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
                                        homolog = 'H2'
                                        spos = int(pos)
                                    elif float(yhat) < -cutoff:
                                        homolog = 'H1'
                                        spos = int(pos)
                                #end of homolog
                                elif (valshift > cutoff and float(yhat) > -cutoff) or (valshift < -cutoff and float(yhat) < cutoff):
                                    epos = pre_pos
                                    #pfline = '%s\t%s\t%s\t%s\t%s\n' % (chrom, spos, epos, pre_homolog, valshift)
                                    #fo.write(pfline)
                                    multiarr = np.append(multiarr, [[pre_chr, spos, epos, pre_homolog, valshift]])
                                #start of homolog
                                elif valshift > -cutoff and float(yhat) > cutoff:
                                    homolog = 'H2'
                                    spos = int(pos)
                                elif valshift < cutoff and float(yhat) < -cutoff:
                                    homolog = 'H1'
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



def Update_Hetboth_Type(paternal_inheritance, phased_hetboth, workdir):
    multiarr = np.empty(0)
    hetbothorder = '%spaternal_infer.hetboth.order.tsv' % (workdir)
    Type = 'M0'

    with open(paternal_inheritance, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('chrom'):
                (chrom, spos, epos, homolog, ar, alarm) = line.strip().split('\t')
                multiarr = np.append(multiarr, [[chrom, spos, epos, homolog, ar, alarm]])
    dim1 = int(len(multiarr)/6)
    multiarr = np.reshape(multiarr, (dim1, 6))
    
    with open(hetbothorder, 'w+') as fo:
        fo.write('#CHROM\tPOS1\tPOS2\tID\tREF\tALT\tP1\tP2\tM1\tM2\tPROF\tPROM\tTYPE\n')
        with open(phased_hetboth, 'rU') as fh1:
            i = 0
            refchr = multiarr[i,0]
            refstart = multiarr[i,1].astype(int)
            refend = multiarr[i,2].astype(int)
            refhomolog = multiarr[i,3]
            for line in fh1:
                if not line.strip().startswith('#'):
                    (CHROM, POS1, POS2, ID, REF, ALT, P1, P2, M1, M2, PROF, PROM, TYPE) = line.strip().split('\t')
                    
                    #SNP resides in identified haplotype block                
                    if CHROM == refchr and int(POS2) >= refstart and int(POS2) <= refend:
                        if refhomolog == 'H1':
                            Type = 'M1'
                        elif refhomolog == 'H2':
                            Type = 'M2'
                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS1, POS2, ID, REF, ALT, P1, P2, M1, M2, PROF, PROM, Type)
                        fo.write(pfline)
                    elif CHROM == refchr and int(POS2) > refend:
                        if i < dim1-1:
                            i += 1
                            refchr = multiarr[i,0]
                            refstart = multiarr[i,1].astype(int)
                            refend = multiarr[i,2].astype(int)
                            refhomolog = multiarr[i,3]
                            if CHROM == refchr and int(POS2) >= refstart and int(POS2) <= refend:
                                if refhomolog == 'H1':
                                    Type = 'M1'
                                elif refhomolog == 'H2':
                                    Type = 'M2'
                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS1, POS2, ID, REF, ALT, P1, P2, M1, M2, PROF, PROM, Type)
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
                                if CHROM == refchr and int(POS2) >= refstart and int(POS2) <= refend:
                                    if refhomolog == 'H1':
                                        Type = 'M1'
                                    elif refhomolog == 'H2':
                                        Type = 'M2'
                                    pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS1, POS2, ID, REF, ALT, P1, P2, M1, M2, PROF, PROM, Type)
                                    fo.write(pfline)


if __name__ == '__main__':
#    Infer_Paternal_Inheritance(sys.argv[1], sys.argv[2])
    Update_Hetboth_Type(sys.argv[1], sys.argv[2], sys.argv[3])
