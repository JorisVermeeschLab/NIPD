#!/usr/bin/env python

'''
for stats use
'''

import numpy as np


def cal_basic_stats(inlist):
    if not inlist:
        inlist_stats = [0, 0, 0, 0, 0]
    else:
        inlist_mean = np.mean(inlist)
        inlist_sd = np.std(inlist)
        inlist_median = np.median(inlist)
        inlist_max = np.nanmax(inlist)
        inlist_min = np.nanmin(inlist)
        inlist_stats = [inlist_mean, inlist_sd, inlist_median, inlist_max, inlist_min]
    return(inlist_stats)
