#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np

def cal_basic_stats(inlist):
    if not inlist:
        inlist_stats = [0, 0, 0, 0, 0]
    else:
        inlist_mean = np.mean(inlist)
        inlist_median = np.median(inlist)
        inlist_max = np.nanmax(inlist)
        inlist_min = np.nanmin(inlist)

        #cal N50 of hap block size                                                                                                                                                   
        inlist.sort(reverse=True)
        halflen = sum(inlist)/2
        N50cut = 0
        inlist_N50 = -1
        for ele in inlist:
            if N50cut > halflen:
                inlist_N50 = ele

                break
            else:
                N50cut += ele

        inlist_stats = [inlist_mean, inlist_median, inlist_max, inlist_min, inlist_N50]
    return(inlist_stats)


def cal_switch_res(resolvedsumfh):
    pref = resolvedsumfh.rsplit('.', 1)[0]
    outfh = pref + '.switchover.resolution'
    (preS, preE, preHL, prechr) = (0, 0, '0', '0')
    (resdist, blocklen) = (0, 0)
    reslist = list()
    blocklist = list()

    with open(outfh, 'w+') as fo:
        fo.write('#----------Predicted Haplotype Block-----------\n')
        fo.write('#chr\tstart\tend\tblocksize\tHL\n')
        with open(resolvedsumfh, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('chr'):
                    (chrom, start, end, resHL, supp) = line.strip().split('\t')
                    if resHL != 'inconclusive':
                        blocklen = int(end) - int(start)
                        pfline = '%s\t%s\t%s\t%s\t%s\n' % (chrom, start, end, blocklen, resHL)
                        fo.write(pfline)
                        blocklist.append(blocklen)

                        if chrom == prechr and resHL != preHL:
                            resdist = int(start) - int(preE)
                            reslist.append(resdist)

                        prechr = chrom
                        preHL = resHL
                        preS = start
                        preE = end

        blockstat = cal_basic_stats(blocklist)
        fo.write('#Haplotype Block Size Summary:\n')
        fo.write('#Mean\tMedian\tMax\tMin\tN50\n')
        print >>fo, '\t'.join(map(str, blockstat))

        resstat = cal_basic_stats(reslist)
        fo.write('#Switchovers Resolution Summary:\n')
        fo.write('#Mean\tMedian\tMax\tMin\tN50\n')
        print >>fo, '\t'.join(map(str, resstat))



if __name__ == '__main__':
    cal_switch_res(sys.argv[1])
