#!/usr/bin/env python

'''
scripts for parsing grandparent-father(/mother) trio vcf
'''

from __future__ import division
import sys
from Helper import mkdir, check_dirslash, get_familyinfo, get_familyID
from Utility_stats import cal_basic_stats
import re
import numpy as np

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


def Evaluate_GrandparTrio(phasedtriovcf, familyID, triopatID, triomatID, triochildID, mindepth, maxdepth, parentaldir):
    grandparTriosum = '%s/%s.%s_%sx.grandpartrio.sum' % (parentaldir, familyID, mindepth, maxdepth)
    t1cnt = t2cnt = t3cnt = t4cnt = t5cnt = 0
    
    with open(grandparTriosum, 'w+') as fo:
        with open(phasedtriovcf, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('##'):
                    if line.strip().startswith('#CHROM'):
                        infoline = line.strip().split('\t')
                        triochild_pos = int(infoline.index(triochildID))
                        triopat_pos = int(infoline.index(triopatID))
                        triomat_pos = int(infoline.index(triomatID))
                    else:
                        spline = line.strip().split('\t')
                        (CHROM, POS, REF, ALT, FORMAT, triopat, triomat, triochild) = (spline[0], spline[1], spline[3], spline[4], spline[8], spline[triopat_pos], spline[triomat_pos], spline[triochild_pos])
                        if len(REF) == 1 and len(ALT) == 1:
                            dp_search = re.search(r'DP', FORMAT)
                            
                            if dp_search:
                                FORMATl = FORMAT.split(':')
                                dp_pos = int(FORMATl.index('DP'))
                                if triopat != './.' and triomat != './.' and triochild != './.':
                                    triopat_dp = triopat.split(':')[dp_pos]
                                    triomat_dp = triomat.split(':')[dp_pos]
                                    triochild_dp = triochild.split(':')[dp_pos]
                                    if triopat_dp != '.' and int(triopat_dp) >= mindepth and int(triopat_dp) <= maxdepth and triomat_dp != '.' and int(triomat_dp) >= mindepth and int(triomat_dp) <= maxdepth and triochild_dp != '.' and int(triochild_dp) >= mindepth and int(triochild_dp) <= maxdepth:
                                        triopat_gt = triopat.split(':')[0]
                                        triomat_gt = triomat.split(':')[0]
                                        triochild_gt = triochild.split(':')[0]
                                        if (triopat_gt == '0|0' and triomat_gt == '1|1' and triochild_gt == '0|1') or (triopat_gt == '1|1' and triomat_gt == '0|0' and triochild_gt == '1|0'):
                                            t1cnt += 1
                                        elif (triopat_gt == '0|1' and triomat_gt == '0|0' and triochild_gt == '0|0') or (triopat_gt == '0|1' and triomat_gt == '1|1' and triochild_gt == '0|1') or (triopat_gt == '1|0' and triomat_gt == '0|0' and triochild_gt == '1|0') or (triopat_gt == '1|0' and triomat_gt == '1|1' and triochild_gt == '1|1'):
                                            t2cnt += 1
                                        elif (triomat_gt == '0|1' and triopat_gt == '0|0' and triochild_gt == '0|0') or (triomat_gt == '0|1' and triopat_gt == '1|1' and triochild_gt == '1|0') or (triomat_gt == '1|0' and triopat_gt == '0|0' and triochild_gt == '0|1') or (triomat_gt == '1|0' and triopat_gt == '1|1' and triochild_gt == '1|1'):
                                            t3cnt += 1
                                        elif (triopat_gt == '0|1' and triomat_gt == '0|1' and triochild_gt == '0|0') or (triopat_gt == '1|0' and triomat_gt == '1|0' and triochild_gt == '1|1'):
                                            t4cnt += 1
                                        elif (triopat_gt == '1|1' and triomat_gt == '1|1' and triochild_gt == '1|1'):
                                            t5cnt += 1
        fo.write('No. Type1 SNPs - grandparents of different homozygous alleles\n')
        fo.write('%s\n' % t1cnt)
        fo.write('No. Type2 SNPs - grandfather heterozygous and grandmother homozygous sites\n')
        fo.write('%s\n' % t2cnt)
        fo.write('No. Type3 SNPs - grandfather homozygous and grandmother heterozygous sites\n')
        fo.write('%s\n' % t3cnt)
        fo.write('No. Type4 SNPs - phased heterozygous sites for both grandparents\n')
        fo.write('%s\n' % t4cnt)
        fo.write('No. Type5 SNPs - variant sites for both grandparents\n')
        fo.write('%s\n' % t5cnt)
                


def Parsing_Grandpat_Phasedrawvcf(phasedtriovcf, allrawvcf, infolist, patgrand, mindepth, maxdepth, workdir):
    workdir = check_dirslash(workdir)
    parentaldir = workdir + '1_parentshap'
    mkdir(parentaldir)
   
    familyinfo = get_familyinfo(infolist)
    triopatID = familyinfo['trio1']
    triomatID = familyinfo['trio2']
    triochildID = familyinfo['trio3']
    unphasedID = familyinfo['unphasedpar']
    familyID = get_familyID(infolist)
           
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)

    #evaluate grandparents trio phasing
    Evaluate_GrandparTrio(phasedtriovcf, familyID, triopatID, triomatID, triochildID, mindepth, maxdepth, parentaldir)

    
    #may need to apply strict threshold when maternal has CNV (long stretch of imbalance); will affect plasma allele solving
    raflower = 0.01
    rafupper = 0.99
    
    infhtrio = phasedtriovcf.rsplit('/', 1)[1]
    infhall = allrawvcf.rsplit('/', 1)[1]
    #Y chromosome snps are included in HeteroPat file; while X chromosome snps are included in HeteroMat file
    HomoDiff = '%s/%s.%s_%sx.homodiff.tsv' % (parentaldir, familyID, mindepth, maxdepth)
    HeteroPat = '%s/%s.%s_%sx.heteropat.tsv' % (parentaldir, familyID, mindepth, maxdepth)
    HeteroMat = '%s/%s.%s_%sx.heteromat.tsv' % (parentaldir, familyID, mindepth, maxdepth)
    HeteroBoth = '%s/%s.%s_%sx.heteroboth.tsv' % (parentaldir, familyID, mindepth, maxdepth)
    HomoSame = '%s/%s.%s_%sx.homosame.tsv' % (parentaldir, familyID, mindepth, maxdepth)
    
    #First read in the parent data that hasn't beed phased to extract homozygous sites
    unphased_par = dict()
    (triopat_pos, triomat_pos, triochild_pos, unphased_pos) = (0, 0, 0, 0)

    with open(allrawvcf, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('##'):
                if line.strip().startswith('#CHROM'):
                    infoline = line.strip().split('\t')
                    unphased_pos = int(infoline.index(unphasedID))
                else:
                    spline = line.strip().split('\t')
                    (CHROM, POS, FORMAT, sampleinfo) = (spline[0], spline[1], spline[8], spline[unphased_pos])
                    
                    #check depth
                    dp_search = re.search(r'DP', FORMAT)
                    if dp_search and sampleinfo != './.':
                        FORMATl = FORMAT.split(':')
                        dp_pos = int(FORMATl.index('DP'))
                        sp_dp = sampleinfo.split(':')[dp_pos]
                        if sp_dp != '.' and int(sp_dp) >= mindepth and int(sp_dp) <= maxdepth:
                            samplegt = sampleinfo.split(':')[0]
                            if samplegt == '0/0' or samplegt == '1/1':
                                refID = '%s:%s' % (CHROM, POS)
                                refval = '%s:%s' % (sp_dp, samplegt)
                                unphased_par[refID] = refval


    with open(phasedtriovcf, 'rU') as fh:
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
                    #(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, paternal, maternal, sibling) = line.strip().split('\t')
                    if line.strip().startswith('#CHROM'):
                        infoline = line.strip().split('\t')
                        triochild_pos = int(infoline.index(triochildID))
                        triopat_pos = int(infoline.index(triopatID))
                        triomat_pos = int(infoline.index(triomatID))
                    else:
                        spline = line.strip().split('\t')
                        (CHROM, POS, ID, REF, ALT, QUAL, FORMAT, triopat, triomat, triochild) = (spline[0], spline[1], spline[2], spline[3], spline[4], spline[5], spline[8], spline[triopat_pos], spline[triomat_pos], spline[triochild_pos])

                        (PATRAF, MATRAF, SIBRAF) = ('.', '.', '.')
                        (PH1, PH2, MH1, MH2, SH1, SH2) = ('.', '.', '.', '.', '.', '.')
                        (TYPE, pat_dp, mat_dp, sib_dp) = ('.', '.', '.', '.')

                        checkID = '%s:%s' % (CHROM, POS)
                        #get unphased parent site
                        if checkID in unphased_par:
                            getval = unphased_par[checkID].split(':')
                            if (patgrand == 'yes'):
                                (mat_dp, mat_gt) = (getval[0], getval[1])
                            else:
                                (pat_dp, pat_gt) = (getval[0], getval[1])

                            #Filter 
                            #- depth threshold
                            #- heterozygous site with reference allele frequency theshold
                            #- only SNPs are considered, indels are filtered as AERC can only compile single position (may further look into this problem)
                            #- may add more filtering
                            if len(REF) == 1 and len(ALT) == 1:
                                dp_search = re.search(r'DP', FORMAT)
                                ad_search = re.search(r'AD', FORMAT)
                                #Y chromosome removed temporarily
                                if dp_search and ad_search and CHROM != 'Y':
                                    #get position of DP in FORMAT, in order to search for the DP values of samples
                                    #format might be GT:AD:DP:GQ:PL or GT:AB:AD:DP:GQ:PL

                                    FORMATl = FORMAT.split(':')
                                    dp_pos = int(FORMATl.index('DP'))
                                    ad_pos = int(FORMATl.index('AD'))

                                    if (triopat != './.') and (triomat != './.') and (triochild != './.'):
                                        triopat_dp = triopat.split(':')[dp_pos]
                                        triomat_dp = triomat.split(':')[dp_pos]
                                        child_dp = triochild.split(':')[dp_pos]
                                        child_gt = triochild.split(':')[0]
                                        childphaseval = list(triochild)[1]
                                        child_ad = triochild.split(':')[ad_pos]

                                        #check whether the site has been phased and has allele depth info
                                        if (childphaseval == '|') and (child_ad != '.'):
                                            (CH1, CH2) = child_gt.split('|')
                                            TYPE = '.'
                                            (child_ad_ref, child_ad_alt) = child_ad.split(',')
                                            childtotad = (int(child_ad_ref)+int(child_ad_alt))

                                            if childtotad != 0:
                                                if (patgrand == 'yes'):
                                                    PATRAF = '{0:.4f}'.format(int(child_ad_ref)/childtotad)
                                                    pat_dp = child_dp
                                                else:
                                                    MATRAF = '{0:.4f}'.format(int(child_ad_ref)/childtotad)
                                                    mat_dp = child_dp


                                                if triopat_dp != '.' and int(triopat_dp) >= mindepth and int(triopat_dp) <= maxdepth and triomat_dp != '.' and int(triomat_dp) >= mindepth and int(triomat_dp) <= maxdepth and child_dp != '.' and int(child_dp) >= mindepth and int(child_dp) <= maxdepth:

                                                    #!indels are also included here, need evaluating whether to include indels sites
                                                    #!only biallelic sites are considered here


                                                    #Case 1 - homozygous but of different alleles &
                                                    #Case 5 - type5 error estimation; in grandparents phasing case, father and mother might be 0/0 and 0/0
                                                    if (child_gt == '0|0' or child_gt == '1|1'):
                                                        if (patgrand == 'yes'):
                                                            if (child_gt == '0|0' and mat_gt == '1/1') or (child_gt == '1|1' and mat_gt == '0/0'):
                                                                PH1 = PH2 = child_gt.split('|')[0]
                                                                MH1 = MH2 = mat_gt.split('/')[0]
                                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                                fo0.write(pfline)
                                                            elif (child_gt == '1|1' and mat_gt == '1/1'):
                                                                PH1 = PH2 = MH1 = MH2 = child_gt.split('|')[0]
                                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                                fo4.write(pfline)

                                                        else:
                                                            if (child_gt == '0|0' and pat_gt == '1/1') or (child_gt == '1|1' and pat_gt == '0/0'):
                                                                MH1 = MH2 = child_gt.split('|')[0]
                                                                PH1 = PH2 = pat_gt.split('/')[0]                                                            
                                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                                fo0.write(pfline)
                                                            elif (child_gt == '1|1' and pat_gt == '1/1'):
                                                                PH1 = PH2 = MH1 = MH2 =child_gt.split('|')[0]
                                                                pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                                fo4.write(pfline)

                                                    #Case 2 - type2: exclude X chromosome when phasing father, father's genotype is unlikely to be heterozygous
                                                    elif ((child_gt == '0|1' or child_gt == '1|0') and patgrand == 'yes' and CHROM != 'X'):
                                                        if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper:
                                                            if (child_gt == '0|1' and mat_gt == '0/0') or (child_gt == '1|0' and mat_gt == '1/1'):
                                                                TYPE = 'P1'
                                                            elif (child_gt == '0|1' and mat_gt == '1/1') or (child_gt == '1|0' and mat_gt == '0/0'):
                                                                TYPE = 'P2'
                                                          
                                                            (PH1, PH2) = (child_gt.split('|')[0], child_gt.split('|')[1])
                                                            MH1 = MH2 = mat_gt.split('/')[0]
                                                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                            fo1.write(pfline)

                                                    #Case3 - type3
                                                            #set mark
                                                    elif ((child_gt == '0|1' or child_gt == '1|0') and patgrand == 'no'):
                                                        if MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                            if (child_gt == '0|1' and pat_gt == '0/0') or (child_gt == '1|0' and pat_gt == '1/1'):
                                                                TYPE = 'M1'
                                                            elif (child_gt == '0|1' and pat_gt == '1/1') or (child_gt == '1|0' and pat_gt == '0/0'):
                                                                TYPE = 'M2'

                                                            (MH1, MH2) = (child_gt.split('|')[0], child_gt.split('|')[1])
                                                            PH1 = PH2 = pat_gt.split('/')[0]
                                                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                            fo2.write(pfline)

                                                    #Case4 - type4 (removed from single grandparents case)
                                                    #elif (child_gt == '0|1' or child_gt == '1|0'):
                                                    #   if PATRAF != '.' and float(PATRAF) >= raflower and float(PATRAF) <= rafupper and MATRAF != '.' and float(MATRAF) >= raflower and float(MATRAF) <= rafupper:
                                                           #TYPE to be assigned after phasing paternal category snps
                                                           #TYPE = 'unknown'
                                                           #if CHROM != 'X' and CHROM != 'Y':
                                                               #pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, pat_dp, mat_dp, sib_dp)
                                                               #fo3.write(pfline)

                                                    #Case5 - type5

                            else:
                                #if no depth information available, discard the site
                                pass
    Sum_Parsedvcf_Stats(HomoDiff)
    Sum_Parsedvcf_Stats(HeteroPat)
    Sum_Parsedvcf_Stats(HeteroMat)
    Sum_Parsedvcf_Stats(HeteroBoth)
    Sum_Parsedvcf_Stats(HomoSame)



if __name__ == '__main__':
    Parsing_Grandpat_Phasedrawvcf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
