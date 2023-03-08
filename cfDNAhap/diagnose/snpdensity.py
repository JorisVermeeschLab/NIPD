#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np


def cal_basic_stats(inlist):
    if not inlist:
        inlist_stats = [0, 0, 0, 0]
    else:
        inlist_mean = np.mean(inlist)
        inlist_median = np.median(inlist)
        inlist_max = np.nanmax(inlist)
        inlist_min = np.nanmin(inlist)
        inlist_stats = [inlist_mean, inlist_median, inlist_max, inlist_min]
    return(inlist_stats)


def DS_SNPdensity(roilist, infh, workdir):
    outfh = workdir + '/SNPdensity.ds.out' 
    roi = list()

    (PM1distl, PM2distl, PM1DPl, PM2DPl, PM1FARl, PM2FARl) = (list(), list(), list(), list(), list(), list())
    
    with open(roilist, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('chr'):
                roi.append(line.strip().split('\t'))
                

    with open(outfh, 'w+') as fo:
        for i in range(0, len(roi)):
            (chrom, start, end) = roi[i]
            start = int(start)
            end = int(end)
            (preposT1, preposT2) = (0, 0)
            (cntT1, cntT2) = (1, 1)
            with open(infh, 'rU') as fh1:
                for line in fh1:
                    if not line.strip().startswith('#'):
                        (CHROM,POS,RSID,REF,ALT,PH1,PH2,MH1,MH2,FAR,SITEINFER,TYPE,TCNT) = line.strip().split('\t')
                        if int(CHROM) <= int(chrom):
                            if CHROM == chrom:
                                if int(POS) >= int(start) and int(POS) <= int(end):
                                    if TYPE == 'P1' or TYPE == 'M1':
                                        if preposT1 != 0:
                                            snpdist = int(POS) - preposT1
                                            PM1distl.append(snpdist)
                                            PM1DPl.append(int(TCNT))
                                            PM1FARl.append(float(FAR))
                                            type1 = TYPE
                                            cntT1 += 1
                                        preposT1 = int(POS)
                                    elif TYPE == 'P2' or TYPE == 'M2':
                                        if preposT2 != 0:
                                            snpdist = int(POS) - preposT2
                                            PM2distl.append(snpdist)
                                            PM2DPl.append(int(TCNT))
                                            PM2FARl.append(float(FAR))
                                            type2 = TYPE
                                            cntT2 += 1
                                        preposT2 = int(POS)


                        else:
                            break

                PM1diststats = cal_basic_stats(PM1distl)
                PM1DPstats = cal_basic_stats(PM1DPl)
                PM1FARstats = cal_basic_stats(PM1FARl)
                PM2diststats = cal_basic_stats(PM2distl)
                PM2DPstats = cal_basic_stats(PM2DPl)
                PM2FARstats = cal_basic_stats(PM2FARl)
                fo.write('################\n')
                fo.write('chr\tstart\tend\n')
                fo.write('%s\t%s\t%s\n' % (chrom, start, end))
                fo.write('Measure\tMean\tMedian\tMax\tMin\n')
                fo.write('%s\n' % type1)
                fo.write('No.SNPs -- %s\n' % cntT1)
                fo.write('Distance\t')
                print >>fo, '\t'.join(map(str, PM1diststats))
                fo.write('Depth\t')
                print >>fo, '\t'.join(map(str, PM1DPstats))
                fo.write('FAR\t')
                print >>fo, '\t'.join(map(str, PM1FARstats))
                fo.write('%s\n' % type2)
                fo.write('No.SNPs -- %s\n' % cntT2)
                fo.write('Distance\t')
                print >>fo, '\t'.join(map(str, PM2diststats))
                fo.write('Depth\t')
                print >>fo, '\t'.join(map(str, PM2DPstats))
                fo.write('FAR\t')
                print >>fo, '\t'.join(map(str, PM2FARstats))
                
                
DS_SNPdensity(sys.argv[1], sys.argv[2], sys.argv[3])
