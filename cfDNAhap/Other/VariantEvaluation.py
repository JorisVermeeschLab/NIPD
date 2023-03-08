#!/usr/bin/env python

'''
Author: Huiwen Che
Created: Jan 24 2018
Use: evaluate variant calling result
- ti/tv ratio
- calculate SNP marker density across each chromosome
'''

import sys
import os
from Helper import mkdir, check_dirslash
import re
import numpy as np

def cal_basic_stats(inlist, chrom):
    inlist_mean = np.mean(inlist)
    inlist_sd = np.std(inlist)
    inlist_median = np.median(inlist)
    inlist_max = np.nanmax(inlist)
    inlist_min = np.nanmin(inlist)
    inlist_stats = '%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, inlist_mean, inlist_sd, inlist_median, inlist_max, inlist_min)
    return(inlist_stats)


def get_distance_perchr(infercff):
    '''input file: sample aerc infercff.tsv 
       currently different type of SNP are in seperate files
    '''
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    infh = infercff.rsplit('/', 1)[1]
    fo_prefix = infh.split('.', 1)[0]
    infile_dir = os.path.dirname(os.path.realpath(infercff))
    infile_dir = check_dirslash(infile_dir)
    sumstats = '%s/%s.snpdist.stats.out' % (infile_dir, fo_prefix)

    with open(sumstats, 'w+') as fo0:
        fo0.write('#Mean\tSD\tMedian\tMax\tMin\n')
        with open(infercff, 'rU') as fh0:
            prev_pos = 0
            snp_distl = []
            all_distl = []
            i = 0
            for line in fh0:
                if not line.strip().startswith('#'):
                    (chro, position) = (line.strip().split()[0], line.strip().split()[1])
                    if chro == chroms[i]:
                        current_pos = int(position)
                        snp_dist = current_pos - prev_pos
                        snp_distl.append(snp_dist)
                        all_distl.append(snp_dist)
                        prev_pos = current_pos
                    else:
                        fo0.write(cal_basic_stats(snp_distl, chroms[i]))
                        i += 1
                        prev_pos = 0
                        snp_distl = []
            fo0.write(cal_basic_stats(all_distl, 'overall'))



if __name__ == '__main__':
    get_distance_perchr(sys.argv[1])
