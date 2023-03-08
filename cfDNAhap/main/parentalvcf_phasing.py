#!/usr/bin/env python

'''
scripts merged from genotype (parsing raw unphased jointcall gvcf) and phasedvcf (parsing phased jointcall gvcf based on genotype)
file parsing starts from phased jointcall gvcf
'''

from __future__ import division
import sys
from Helper import mkdir, check_dirslash
from Utility_stats import cal_basic_stats
import re
import numpy as np


def Get_SiblingGender(pedfile):
    sib_gender = 'unknown'
    with open(pedfile, 'rU') as fh:
        for line in fh:
            if line.strip().split()[2] != '0':
                sib_gender = line.strip().split()[4]
                if sib_gender == '1':
                    sib_gender = 'male'
                elif sib_gender == '2':
                    sib_gender = 'female'
                else:
                    sib_gender = 'unknown'
    if sib_gender == 'unknown':
        raise ValueError('The sibling\'s gender is unknown.')
    return(sib_gender)
                

def Sum_Parsedvcf_Stats(fhtype):
    (fo_prefix, snptype) = (fhtype.rsplit('.', 2)[0], fhtype.rsplit('.', 2)[1])
    statsfh = '%s.summary' % (fo_prefix)

    (tcnt, subA0, subA1, subB0, subB1, prepos, dist) = (0, 0, 0, 0, 0, 0, 0)
    (prafl, mrafl, chrdistl, alldistl) = ([], [], [], [])
    distarr = np.empty(0)
    prechr = '0'

    with open(statsfh, 'a') as fo:
        fo.write('##Summary for %s SNPs\n' % (snptype))
        with open(fhtype, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                    tcnt += 1
                    if PATRAF != '.':
                        prafl.append(float(PATRAF))
                    if MATRAF != '.':
                        mrafl.append(float(MATRAF))
                    if TYPE != '.' or TYPE != 'unknown':
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
            prafstat = cal_basic_stats(prafl)
            mrafstat = cal_basic_stats(mrafl)

            fo.write('#TotalNoSNPs\tP1_0/M1_0\tP1_1/M1_1\tP2_0/M2_0\tP2_1/P2_1\n')
            fo.write('%s\t%s\t%s\t%s\t%s\n' % (tcnt, subA0, subA1, subB0, subB1))
            fo.write('#PaternalRAF\n#MaternalRAF\n#Mean\tS.D.\tMedian\tMax\tMin\n')
            print >>fo, '\t'.join(map(str, prafstat))
            print >>fo, '\t'.join(map(str, mrafstat))
            fo.write('#Distance between SNPs\n#Mean\tS.D.\tMedian\tMax\tMin\noverall\t')
            print >>fo, '\t'.join(map(str, alldiststat))
            for i in range(dim0):
                print >>fo, '\t'.join(distarr[i])

        



def Parsing_Phasedrawvcf(phasedrawvcf, pedfile, mindepth, maxdepth, workdir, ffsex):
    check_dirslash(workdir)
    parentaldir = workdir + '1_parentshap'
    mkdir(parentaldir)
    
    try:
        sib_gender = Get_SiblingGender(pedfile)
    except ValueError:
        print 'Program Interrupted!\n\nThe sibling\'s gender is unknown. Please check it.'

    triopatID = 'NA'
    triomatID = 'NA'
    triochildID = 'NA'
    pedfh = open(pedfile, 'rU')
    pedfhlines = pedfh.readlines()
    pedinfoline = pedfhlines[2]
    pedinfol = pedinfoline.strip().split()
    (triopatID, triomatID, triochildID) = (pedinfol[2], pedinfol[3], pedinfol[1])

    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    #may need to apply strict threshold when maternal has CNV (long stretch of imbalance); will affect plasma allele solving
    raflower = 0.01
    rafupper = 0.99
    
    infh = phasedrawvcf.rsplit('/', 1)[1]
    fo_prefix = infh.rsplit('.', 1)[0]
    #Y chromosome snps are included in HeteroPat file; while X chromosome snps are included in HeteroMat file
    HomoDiff = '%s/%s.%s_%sx.homodiff.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroPat = '%s/%s.%s_%sx.heteropat.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroMat = '%s/%s.%s_%sx.heteromat.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroBoth = '%s/%s.%s_%sx.heteroboth.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HomoSame = '%s/%s.%s_%sx.homosame.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    #SNPqc = '%s/%s.%s_%sx.qc.tsv' % (parentaldir, mindepth, maxdepth, fo_prefix)
    #StatsSum = '%s/%s.%s_%sx.stats.out' % (parentaldir, mindepth, maxdepth, fo_prefix)
    
    with open(phasedrawvcf, 'rU') as fh:
        with open(HomoDiff, 'w+') as fo0, open(HeteroPat, 'w+') as fo1, open(HeteroMat, 'w+') as fo2, open(HeteroBoth, 'w+') as fo3, open(HomoSame, 'w+') as fo4:
            headerline = '#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n'
            fo0.write('##SNPs category (Type 1) - parents of different homozygous alleles\n')
            fo0.write(headerline)
            fo1.write('##SNPs category (Type 2) - paternal heterozygous alleles and maternal homozygous alleles\n')
            fo1.write(headerline)
            fo2.write('##SNPs category (Type 3) - paternal homozygous alleles and maternal heterozygous alleles\n')
            fo2.write(headerline)
            fo3.write('##SNPs category (Type 4) - both parents of heterozygous alleles\n')
            fo3.write(headerline)
            fo4.write('##SNPs category (Type 5) - parents of the same mutant alleles\n')
            fo4.write(headerline)
            
            
            for line in fh:
                if not line.strip().startswith('##'):
                    #handling three sample (parents-child trio) gvcf
                    if line.strip().startswith('#CHROM'):
                        infoline = line.strip().split('\t')
                        triochild_pos = int(infoline.index(triochildID))
                        triopat_pos = int(infoline.index(triopatID))
                        triomat_pos = int(infoline.index(triomatID))
                    else:
                        spline = line.strip().split('\t')
                        (CHROM, POS, ID, REF, ALT, QUAL, FORMAT, paternal, maternal, sibling) = (spline[0], spline[1], spline[2], spline[3], spline[4], spline[5], spline[8], spline[triopat_pos], spline[triomat_pos], spline[triochild_pos])
                    
                        #Filter 
                        #- depth threshold
                        #- heterozygous site with reference allele frequency theshold
                        #- only SNPs are considered, indels are filtered as AERC can only compile single position (may further look into this problem)
                        #- may add more filtering
                        if len(REF) == 1 and len(ALT) == 1:
                            dp_search = re.search(r'DP', FORMAT)
                            ad_search = re.search(r'AD', FORMAT)

                            if dp_search and ad_search:
                                #get position of DP in FORMAT, in order to search for the DP values of samples
                                #format might be GT:AD:DP:GQ:PL or GT:AB:AD:DP:GQ:PL

                                FORMATl = FORMAT.split(':')
                                dp_pos = int(FORMATl.index('DP'))
                                ad_pos = int(FORMATl.index('AD'))

                                #if Y chromosome, handle seperately for male fetus
                                if (CHROM == 'Y') and (paternal != './.') and (ffsex == 'male'):
                                    pat_dp = paternal.split(':')[dp_pos]
                                    if pat_dp != '.' and int(pat_dp) >= mindepth and int(pat_dp) <= maxdepth:
                                        pat_ad = paternal.split(':')[ad_pos]
                                        patphaseval = list(paternal)[1]
                                        #check whether the site has been phased
                                        if (patphaseval == '|') and (pat_ad != '.'):
                                            (pat_ad_ref, pat_ad_alt) = pat_ad.split(',')
                                            pattotad = int(pat_ad_ref) + int(pat_ad_alt)
                                            if pattotad != 0:
                                                PATRAF = '{0:.4f}'.format(int(pat_ad_ref)/pattotad)
                                                #need to handle male sibling here (not handled currently)
                                                MATRAF = SIBRAF = '.'
                                                pat_gt = paternal.split(':')[0]
                                                (PH1, PH2) = pat_gt.split('|')
                                                MH1 = MH2 = SH1 = SH2 = '.'
                                                TYPE = 'P2'
                                                MATDP = SIBDP = '.'
                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, MATDP, SIBDP)
                                                fo1.write(pfline)

                                elif (paternal != './.') and (maternal != './.') and (sibling != './.'):
                                    pat_dp = paternal.split(':')[dp_pos]
                                    mat_dp = maternal.split(':')[dp_pos]
                                    sib_dp = sibling.split(':')[dp_pos]
                                    pat_gt = paternal.split(':')[0]
                                    mat_gt = maternal.split(':')[0]
                                    sib_gt = sibling.split(':')[0]
                                    patphaseval = list(paternal)[1]
                                    matphaseval = list(maternal)[1]
                                    sibphaseval = list(sibling)[1]
                                    pat_ad = paternal.split(':')[ad_pos]
                                    mat_ad = maternal.split(':')[ad_pos]
                                    sib_ad = sibling.split(':')[ad_pos]

                                    #check whether the site has been phased and has allele depth info
                                    if (patphaseval == '|') and (matphaseval == '|') and (sibphaseval == '|') and (pat_ad != '.') and (mat_ad != '.') and (sib_ad != '.'):
                                        (PH1, PH2) = pat_gt.split('|')
                                        (MH1, MH2) = mat_gt.split('|')
                                        (SH1, SH2) = sib_gt.split('|')
                                        TYPE = '.'
                                        (pat_ad_ref, pat_ad_alt) = pat_ad.split(',')
                                        (mat_ad_ref, mat_ad_alt) = mat_ad.split(',')
                                        (sib_ad_ref, sib_ad_alt) = sib_ad.split(',')
                                        pattotad = (int(pat_ad_ref)+int(pat_ad_alt))
                                        mattotad = (int(mat_ad_ref)+int(mat_ad_alt))
                                        sibtotad = (int(sib_ad_ref)+int(sib_ad_alt))

                                        if pattotad != 0 and mattotad != 0 and sibtotad != 0:
                                            PATRAF = '{0:.4f}'.format(int(pat_ad_ref)/pattotad)
                                            MATRAF = '{0:.4f}'.format(int(mat_ad_ref)/mattotad)
                                            SIBRAF = '{0:.4f}'.format(int(sib_ad_ref)/sibtotad)


                                            if pat_dp != '.' and int(pat_dp) >= mindepth and int(pat_dp) <= maxdepth and mat_dp != '.' and int(mat_dp) >= mindepth and int(mat_dp) <= maxdepth and sib_dp != '.' and int(sib_dp) >= mindepth and int(sib_dp) <= maxdepth:

                                                #!indels are also included here, need evaluating whether to include indels sites
                                                #!only biallelic sites are considered here

                                                #Case 1 - homozygous but of different alleles
                                                if (pat_gt == '0|0' and mat_gt == '1|1') or (pat_gt == '1|1' and mat_gt == '0|0'):

                                                    if CHROM != 'Y':
                                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                        fo0.write(pfline)

                                                #Case 2 - type2
                                                elif (pat_gt == '0|1' and mat_gt == '0|0') or (pat_gt == '1|0' and mat_gt == '1|1') or (pat_gt == '0|1' and mat_gt == '1|1') or (pat_gt == '1|0' and mat_gt == '0|0'):

                                                    if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper:
                                                        if sib_gt == '0|0' or sib_gt == '1|1':
                                                            TYPE = 'P1'
                                                        elif sib_gt == '0|1' or sib_gt == '1|0':
                                                            TYPE = 'P2'
                                                        if CHROM != 'X' and CHROM != 'Y':
                                                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                            fo1.write(pfline)

                                                #Case3 - type3
                                                elif (pat_gt == '0|0' and mat_gt == '0|1') or (pat_gt == '1|1' and mat_gt == '1|0') or (pat_gt == '1|1' and mat_gt == '0|1') or (pat_gt == '0|0' and mat_gt == '1|0'):

                                                    if MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                        if sib_gt == '0|0' or sib_gt == '1|1':
                                                            TYPE = 'M1'
                                                        elif sib_gt == '1|0' or sib_gt == '0|1':
                                                            TYPE = 'M2'
                                                        if CHROM != 'Y':
                                                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                            fo2.write(pfline)

                                                #Case4 - type4
                                                elif (pat_gt == '0|1' and mat_gt == '0|1') or (pat_gt == '1|0' and mat_gt == '1|0'):

                                                   if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper and MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                       #TYPE to be assigned after phasing paternal category snps
                                                       TYPE = 'unknown'
                                                       if CHROM != 'X' and CHROM != 'Y':
                                                           pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                           fo3.write(pfline)

                                                #Case5 - type5
                                                elif pat_gt == '1|1' and mat_gt == '1|1' and sib_gt == '1|1':
                                                    if CHROM != 'Y':
                                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                        fo4.write(pfline)

                            else:
                                #if no depth information available, discard the site
                                pass
    Sum_Parsedvcf_Stats(HomoDiff)
    Sum_Parsedvcf_Stats(HeteroPat)
    Sum_Parsedvcf_Stats(HeteroMat)
    Sum_Parsedvcf_Stats(HeteroBoth)
    Sum_Parsedvcf_Stats(HomoSame)



if __name__ == '__main__':
    Parsing_Phasedrawvcf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]) 
    
    
    

'''
def Parsing_Phasedrawvcf(phasedrawvcf, pedfile, mindepth, maxdepth, workdir, ffsex):
    check_dirslash(workdir)
    parentaldir = workdir + '1_parentshap'
    mkdir(parentaldir)
    
    try:
        sib_gender = Get_SiblingGender(pedfile)
    except ValueError:
        print 'Program Interrupted!\n\nThe sibling\'s gender is unknown. Please check it.'
        
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    raflower = 0.3
    rafupper = 0.7
    
    infh = phasedrawvcf.rsplit('/', 1)[1]
    fo_prefix = infh.rsplit('.', 1)[0]
    #Y chromosome snps are included in HeteroPat file; while X chromosome snps are included in HeteroMat file
    HomoDiff = '%s/%s.%s_%sx.homodiff.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroPat = '%s/%s.%s_%sx.heteropat.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroMat = '%s/%s.%s_%sx.heteromat.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HeteroBoth = '%s/%s.%s_%sx.heteroboth.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    HomoSame = '%s/%s.%s_%sx.homosame.tsv' % (parentaldir, fo_prefix, mindepth, maxdepth)
    #SNPqc = '%s/%s.%s_%sx.qc.tsv' % (parentaldir, mindepth, maxdepth, fo_prefix)
    #StatsSum = '%s/%s.%s_%sx.stats.out' % (parentaldir, mindepth, maxdepth, fo_prefix)
    
    with open(phasedrawvcf, 'rU') as fh:
        with open(HomoDiff, 'w+') as fo0, open(HeteroPat, 'w+') as fo1, open(HeteroMat, 'w+') as fo2, open(HeteroBoth, 'w+') as fo3, open(HomoSame, 'w+') as fo4:
            headerline = '#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\n'
            fo0.write('##SNPs category (Type 1) - parents of different homozygous alleles\n')
            fo0.write(headerline)
            fo1.write('##SNPs category (Type 2) - paternal heterozygous alleles and maternal homozygous alleles\n')
            fo1.write(headerline)
            fo2.write('##SNPs category (Type 3) - paternal homozygous alleles and maternal heterozygous alleles\n')
            fo2.write(headerline)
            fo3.write('##SNPs category (Type 4) - both parents of heterozygous alleles\n')
            fo3.write(headerline)
            fo4.write('##SNPs category (Type 5) - parents of the same mutant alleles\n')
            fo4.write(headerline)
            
            
            for line in fh:
                if not line.strip().startswith('#'):
                    #handling three sample (parents-child trio) gvcf
                    (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, paternal, maternal, sibling) = line.strip().split('\t')
                    
                    #Filter 
                    #- depth threshold
                    #- heterozygous site with reference allele frequency theshold
                    #- only SNPs are considered, indels are filtered as AERC can only compile single position (may further look into this problem)
                    #- may add more filtering
                    if len(REF) == 1 and len(ALT) == 1:
                        dp_search = re.search(r'DP', FORMAT)

                        if dp_search:
                            #get position of DP in FORMAT, in order to search for the DP values of samples
                            #format might be GT:AD:DP:GQ:PL or GT:AB:AD:DP:GQ:PL

                            FORMATl = FORMAT.split(':')
                            dp_pos = int(FORMATl.index('DP'))
                            #if Y chromosome, handle seperately for male fetus
                            if (CHROM == 'Y') and (paternal != './.') and (ffsex == 'male'):
                                pat_dp = paternal.split(':')[dp_pos]
                                if pat_dp != '.' and int(pat_dp) >= mindepth and int(pat_dp) <= maxdepth:
                                    patphaseval = list(paternal)[1]
                                    #check whether the site has been phased
                                    if patphaseval == '|':
                                        PATRAF = MATRAF = SIBRAF = '.'
                                        pat_gt = paternal.split(':')[0]
                                        (PH1, PH2) = pat_gt.split('|')
                                        MH1 = MH2 = SH1 = SH2 = '.'
                                        TYPE = 'P2'
                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                        fo1.write(pfline)

                            elif (paternal != './.') and (maternal != './.') and (sibling != './.'):
                                pat_dp = paternal.split(':')[dp_pos]
                                mat_dp = maternal.split(':')[dp_pos]
                                sib_dp = sibling.split(':')[dp_pos]
                                pat_gt = paternal.split(':')[0]
                                mat_gt = maternal.split(':')[0]
                                sib_gt = sibling.split(':')[0]
                                patphaseval = list(paternal)[1]
                                matphaseval = list(maternal)[1]
                                sibphaseval = list(sibling)[1]
                                #check whether the site has been phased
                                if patphaseval == '|' and matphaseval == '|' and sibphaseval == '|':
                                    (PH1, PH2) = pat_gt.split('|')
                                    (MH1, MH2) = mat_gt.split('|')
                                    (SH1, SH2) = sib_gt.split('|')
                                    PATRAF = MATRAF = SIBRAF = TYPE = '.'


                                    if pat_dp != '.' and int(pat_dp) >= mindepth and int(pat_dp) <= maxdepth and mat_dp != '.' and int(mat_dp) >= mindepth and int(mat_dp) <= maxdepth and sib_dp != '.' and int(sib_dp) >= mindepth and int(sib_dp) <= maxdepth:
                                        ab_search = re.search(r'AB', FORMAT)

                                        #!indels are also included here, need evaluating whether to include indels sites
                                        #!only biallelic sites are considered here

                                        #Case 1 - homozygous but of different alleles
                                        if (pat_gt == '0|0' and mat_gt == '1|1') or (pat_gt == '1|1' and mat_gt == '0|0'):
                                            if ab_search:
                                                ab_pos = int(FORMATl.index('AB'))
                                                #RAF threshold not applied here, only applied on parents
                                                SIBRAF = sibling.split(':')[ab_pos]
                                            if CHROM != 'Y':
                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                                fo0.write(pfline)

                                        #Case2 - type2
                                        elif (pat_gt == '0|1' and mat_gt == '0|0') or (pat_gt == '1|0' and mat_gt == '1|1') or (pat_gt == '0|1' and mat_gt == '1|1') or (pat_gt == '1|0' and mat_gt == '0|0'):
                                            if ab_search:
                                                ab_pos = int(FORMATl.index('AB'))
                                                PATRAF = paternal.split(':')[ab_pos]
                                                SIBRAF = sibling.split(':')[ab_pos]
                                                if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper:
                                                    if sib_gt == '0|0' or sib_gt == '1|1':
                                                        TYPE = 'P1'
                                                    elif sib_gt == '0|1' or sib_gt == '1|0':
                                                        TYPE = 'P2'
                                                    if CHROM != 'X' and CHROM != 'Y':
                                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                                        fo1.write(pfline)

                                        #Case3 - type3
                                        elif (pat_gt == '0|0' and mat_gt == '0|1') or (pat_gt == '1|1' and mat_gt == '1|0') or (pat_gt == '1|1' and mat_gt == '0|1') or (pat_gt == '0|0' and mat_gt == '1|0'):
                                            if ab_search:
                                                ab_pos = int(FORMATl.index('AB'))
                                                MATRAF = maternal.split(':')[ab_pos]
                                                SIBRAF = sibling.split(':')[ab_pos]
                                                if MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                    if sib_gt == '0|0' or sib_gt == '1|1':
                                                        TYPE = 'M1'
                                                    elif sib_gt == '1|0' or sib_gt == '0|1':
                                                        TYPE = 'M2'
                                                    if CHROM != 'Y':
                                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                                        fo2.write(pfline)

                                        #Case4 - type4
                                        elif (pat_gt == '0|1' and mat_gt == '0|1') or (pat_gt == '1|0' and mat_gt == '1|0'):
                                            if ab_search:
                                               ab_pos = int(FORMATl.index('AB'))
                                               PATRAF = paternal.split(':')[ab_pos]
                                               MATRAF = maternal.split(':')[ab_pos]
                                               SIBRAF = sibling.split(':')[ab_pos]
                                               if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper and MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                   #TYPE to be assigned after phasing paternal category snps
                                                   TYPE = 'unknown'
                                                   if CHROM != 'X' and CHROM != 'Y':
                                                       pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                                       fo3.write(pfline)

                                        #Case5 - type5
                                        elif pat_gt == '1|1' and mat_gt == '1|1' and sib_gt == '1|1':
                                            if CHROM != 'Y':
                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE)
                                                fo4.write(pfline)

                        else:
                            #if no depth information available, discard the site
                            pass
    Sum_Parsedvcf_Stats(HomoDiff)
    Sum_Parsedvcf_Stats(HeteroPat)
    Sum_Parsedvcf_Stats(HeteroMat)
    Sum_Parsedvcf_Stats(HeteroBoth)
    Sum_Parsedvcf_Stats(HomoSame)
'''
