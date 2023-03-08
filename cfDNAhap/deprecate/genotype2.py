#!/usr/bin/env python

'''
Author: Huiwen Che
Created: December 29 2017
Use: parsing gatk haplotypecaller variant calling joint vcf file, extracting raw genotype information of 2 samples (parents) 
'''

import sys
import os
from Helper import mkdir, check_dirslash
import re
import numpy as np


def cal_basic_stats(inlist):
    inlist_mean = np.mean(inlist)
    inlist_sd = np.std(inlist)
    inlist_median = np.median(inlist)
    inlist_max = np.nanmax(inlist)
    inlist_min = np.nanmin(inlist)
    inlist_stats = '#Mean\tSD\tMedian\tMax\tMin\n%s\t%s\t%s\t%s\t%s\n' % (inlist_mean, inlist_sd, inlist_median, inlist_max, inlist_min)
    return(inlist_stats)


def get_allele_cate(infile):
    infh = infile.rsplit('/', 1)[1]
    fo_prefix = infh.split('.', 1)[0]
    infile_dir = os.path.dirname(os.path.realpath(infile))
    '''mkdir for genotype of parents data'''
    infile_dir = check_dirslash(infile_dir)    
    mindepth = 15
    maxdepth = 500
    gt_dir = '%sgenotype_parents_%sx_%sx' % (infile_dir, mindepth, maxdepth)
    mkdir(gt_dir)    

    HomoDiff = '%s/%s.homodiff.bed' % (gt_dir, fo_prefix)
    HeteroPat = '%s/%s.heteropat.bed' % (gt_dir, fo_prefix)
    HeteroMat = '%s/%s.heteromat.bed' % (gt_dir, fo_prefix)
    HomoSame = '%s/%s.homosame.bed' % (gt_dir, fo_prefix)
    HeteroBoth = '%s/%s.heteroboth.bed' % (gt_dir, fo_prefix)
    StatsSum = '%s/%s.stats.out' % (gt_dir, fo_prefix)

    with open(infile, 'rU') as fh:

        '''output bed format, pos1 = pos-1 where pos is the POS in vcf file
        '''
        with open(HomoDiff, 'w+') as fo0:
            fo0.write('##Different Homozygous Alleles\n')
            fo0.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n')
            with open(HeteroPat, 'w+') as fo1:
                fo1.write('##Heterozygous Paternal; Homozygous Maternal\n')
                fo1.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n')
                with open(HeteroMat, 'w+') as fo2:
                    fo2.write('##Homozygous Paternal; Heterozygous Maternal\n')
                    fo2.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n')
                    with open(HomoSame, 'w+') as fo3:
                        fo3.write('##Variant Homozygous Paternal; Variant Homozygous Maternal\n')
                        fo3.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n')
                        with open(HeteroBoth, 'w+') as fo4:
                            fo4.write('##Heterozygous Paternal; Heterozygous Maternal\n')
                            fo4.write('#CHROM\tPOS\tRSID\tREF\tALT\tQUAL\tPATRAF\tMATRAF\tSIBRAF\tPH1\tPH2\tMH1\tMH2\tSH1\tSH2\tTYPE\tPATDP\tMATDP\tSIBDP\n')
                            with open(StatsSum, 'w+') as fo5:
                                fo5.write('##Summary stats of parents genotype data\n')
                                (homodiff_cnt, heteropat_cnt, heteromat_cnt, homosame_cnt, heteroboth_cnt) = (0, 0, 0, 0, 0)
                                (heteropat_abl, heteromat_abl, heteroboth_patabl, heteroboth_matabl) = ([], [], [], [])

                                for line in fh:
                                    if not line.strip().startswith('#'):
                                        #handling two sample (parents) vcf
                                        (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, paternal, maternal) = line.strip().split('\t')
                                
                                        '''Filter - filter out sites where either of the sample has depth < 15 and depth > 500
                                        '''
                                        #possibly more filtering or better filtering criterion?
                                        dp_search = re.search(r'DP', FORMAT)
                                        
                                        if dp_search:
                                            '''get position of DP in FORMAT, in order to search for the DP values of samples
                                               format might be GT:AD:DP:GQ:PL or GT:AB:AD:DP:GQ:PL
                                            '''                                  
                                            FORMATl = FORMAT.split(':')
                                            dp_pos = int(FORMATl.index('DP'))
                                            if (paternal != './.') and (maternal != './.'):
                                                pat_dp = paternal.split(':')[dp_pos]
                                                mat_dp = maternal.split(':')[dp_pos]
                                                
                                                if pat_dp != '.' and int(pat_dp) >= mindepth and int(pat_dp) <= maxdepth and mat_dp != '.' and int(mat_dp) >= mindepth and int(mat_dp) <= maxdepth:
                                                    pat_gt = paternal.split(':')[0]
                                                    pat_gtl = pat_gt.split('/')
                                                    mat_gt = maternal.split(':')[0]
                                                    mat_gtl = mat_gt.split('/')
                                                    pos1 = int(POS) - 1
                                                    printline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (CHROM, POS, ID, REF, ALT, QUAL, '.', '.', '.')
                                                
                                                    #indels are also included here, need evaluating whether to include indels sites
                                                    #only biallelic sites considered at the moment
                                                    '''Case 1 - homozygous but of different alleles
                                                    '''
                                                    if (pat_gt == '0/0' and mat_gt == '1/1') or (pat_gt == '1/1' and mat_gt == '0/0'):
                                                        PH1 = PH2 = pat_gt.split('/')[0]
                                                        MH1 = MH2 = mat_gt.split('/')[0]
                                                        fpline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (printline, PH1, PH2, MH1, MH2, '.', '.', '.', '.', '.', '.')
                                                        fo0.write(fpline)
                                                        homodiff_cnt += 1

                                                    #Case 2 - heterozygous paternal and homozygous maternal
                                                    elif (pat_gt == '0/1' and mat_gt == '0/0') or (pat_gt == '0/1' and mat_gt == '1/1'):
                                            
                                                        '''extract allele balance values from heterozygous calls
                                                        '''
                                                        ab_search = re.search(r'AB', FORMAT)
                                                        if ab_search:
                                                            ab_pos = int(FORMATl.index('AB'))
                                                            heteropat_ab = paternal.split(':')[ab_pos]
                                                            heteropat_abl.append(float(heteropat_ab))
                                                            fpline = '%s\t%s\n' % (printline, heteropat_ab)
                                                        else:
                                                            fpline = '%s\t%s\n' % (printline, 'N/A')
                                                        fo1.write(fpline)
                                                    
                                                        heteropat_cnt += 1

                                                    #Case 3 - homozygous paternal; heterozygous maternal
                                                    elif (pat_gt == '0/0' and mat_gt == '0/1') or (pat_gt == '1/1' and mat_gt == '0/1'):
                                                        ab_search = re.search(r'AB', FORMAT)
                                                        if ab_search:
                                                            ab_pos = int(FORMATl.index('AB'))
                                                            heteromat_ab = maternal.split(':')[ab_pos]
                                                            heteromat_abl.append(float(heteromat_ab))
                                                            fpline = '%s\t%s\n' % (printline, heteromat_ab)
                                                        else:
                                                            fpline = '%s\t%s\n' % (printline, 'N/A')
                                                        fo2.write(fpline)
                                                        heteromat_cnt += 1

                                                    #Case 4 - variant homozygous paternal; variant homozygous maternal
                                                    elif (pat_gt == '1/1' and mat_gt == '1/1'):
                                                        fpline = '%s\t%s\n' % (printline, 'N/A')
                                                        fo3.write(fpline)
                                                        homosame_cnt += 1

                                                    #Case 5 - heterozygous paternal; heterozygous maternal; cff_infer_gt needs modification pat0/1 mat0/1 result in 3 types
                                                    else:
                                                        ab_search = re.search(r'AB', FORMAT)
                                                        if ab_search:
                                                            ab_pos = int(FORMATl.index('AB'))
                                                            heteroboth_patab = paternal.split(':')[ab_pos]
                                                            heteroboth_matab = maternal.split(':')[ab_pos]
                                                            heteroboth_patabl.append(float(heteroboth_patab))
                                                            heteroboth_matabl.append(float(heteroboth_matab))
                                                            heteroboth_ab = '%s:%s' % (heteroboth_patab, heteroboth_matab)
                                                            fpline = '%s\t%s\n' % (printline, heteroboth_ab)
                                                        else:
                                                            fpline = '%s\t%s\n' % (printline, 'N/A')
                                                        fo4.write(fpline)
                                                        heteroboth_cnt += 1
                                            else:
                                                pass
                                        else:
                                            print line
                                                                                                                                 
 
                                fo5.write('#Total number of different homozygous sites:\n')
                                fo5.write('%s\n' % (homodiff_cnt))
                                fo5.write('#Total number of hetero paternal and homo maternal sites:\n')
                                fo5.write('%s\n' % (heteropat_cnt))
                                fo5.write('#Total number of homo paternal and hetero maternal sites:\n')
                                fo5.write('%s\n' % (heteromat_cnt))
                                fo5.write('#Total number of variant homozygous sites:\n')
                                fo5.write('%s\n' % (homosame_cnt))
                                fo5.write('#Total number of hetero paternal and hetero maternal sites:\n')
                                fo5.write('%s\n' % (heteroboth_cnt))
                                
                                fo5.write('###Allele Balance Stats\n')
                                fo5.write('##Heterozygous paternal:\n')
                                fo5.write(cal_basic_stats(heteropat_abl))

                                fo5.write('##Heterozygous maternal:\n')
                                fo5.write(cal_basic_stats(heteromat_abl))

                                fo5.write('##Hetero both paternal:\n')
                                fo5.write(cal_basic_stats(heteroboth_patabl))

                                fo5.write('##Hetero both maternal:\n')
                                fo5.write(cal_basic_stats(heteroboth_matabl))


if __name__ == '__main__':
    get_allele_cate(sys.argv[1])
