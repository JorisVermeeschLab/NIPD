#!/usr/bin/env python

from __future__ import division
import sys
import os
from subprocess import call

def cal_acc2snparr(snparrfh, cffh, workdir, samplename):
    hbolp = workdir + '/' + samplename + '.hb2snparr.olp'
    accfo = workdir + '/' + samplename + '.acc2snparr.out'
    chrlist = map(str, range(1, 23))
    if snparrfh.endswith('mat.resolve.2.block'):
        chrlist.append('X')
    totrefhb = 0
    chrrefhb = dict()
    
    hp1arr = workdir + '/snparr.hp1.tmp'
    hp2arr = workdir + '/snparr.hp2.tmp'
    with open(snparrfh, 'rU') as fh1, open(hp1arr, 'w+') as fo1, open(hp2arr, 'w+') as fo2:
        for line in fh1:
            if not line.startswith('chr'):
                (chrom, start, end, hp) = line.strip().split('\t')
                if chrom in chrlist:
                    if int(hp) == 1:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, hp)
                        fo1.write(pfline)
                    elif int(hp) == 2:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, hp)
                        fo2.write(pfline)
                    blocksize = int(end) - int(start)
                    totrefhb += blocksize
                    if chrom in chrrefhb:
                        chrrefhb[chrom] += blocksize
                    else:
                        chrrefhb[chrom] = blocksize
    
    hp1cf = workdir + '/cf.hp1.tmp'
    hp2cf = workdir + '/cf.hp2.tmp'
    with open(cffh, 'rU') as fh2, open(hp1cf, 'w+') as fo3, open(hp2cf, 'w+') as fo4:
        for line in fh2:
            if not line.startswith('chr'):
                (chrom, start, end, hp) = line.strip().split('\t')
                if chrom in chrlist:
                    #if hp == 'MH1' or hp == 'PH1':
                    if int(hp) == 1:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, hp)
                        fo3.write(pfline)
                    #elif hp == 'MH2' or hp == 'PH2':
                    elif int(hp) == 2:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, hp)
                        fo4.write(pfline)
           
    
    hp1arr_sort = workdir + '/snparr.hp1.sorted.tmp'
    hp2arr_sort = workdir + '/snparr.hp2.sorted.tmp'
    hp1cf_sort = workdir + '/cf.hp1.sorted.tmp'
    hp2cf_sort = workdir + '/cf.hp2.sorted.tmp'

    shcmd1 = 'bedtools sort -i %s > %s' % (hp1arr, hp1arr_sort)
    shcmd2 = 'bedtools sort -i %s > %s' % (hp1cf, hp1cf_sort)
    shcmd12 = 'bedtools intersect -a %s -b %s > %s' % (hp1arr_sort, hp1cf_sort, hbolp)
    call(shcmd1, shell=True)
    call(shcmd2, shell=True)
    call(shcmd12, shell=True)
    
    shcmd3 = 'bedtools sort -i %s > %s' % (hp2arr, hp2arr_sort)
    shcmd4 = 'bedtools sort -i %s > %s' % (hp2cf, hp2cf_sort)
    shcmd34 = 'bedtools intersect -a %s -b %s >> %s' % (hp2arr_sort, hp2cf_sort, hbolp)
    call(shcmd3, shell=True)
    call(shcmd4, shell=True)
    call(shcmd34, shell=True)
    
    shcmd5 = 'rm *.tmp'
    call(shcmd5, shell=True)
    
    olptotsize = 0
    olpchrsize = dict()
    
    with open(hbolp, 'rU') as fh3:
        for line in fh3:
            (chrom, start, end, val) = line.strip().split('\t')
            olpsize = int(end) - int(start)
            olptotsize += olpsize
            chrom = str(chrom)
            if chrom in olpchrsize:
                olpchrsize[chrom] += olpsize
            else:
                olpchrsize[chrom] = olpsize
    
    with open(accfo, 'w+') as fo5:
        fo5.write("#Sample: %s\n" % (samplename))
        pfline = 'Overall accuracy: %s\n' % (olptotsize/totrefhb)
        fo5.write(pfline)
        for i in olpchrsize:
            pfline = 'Chr%s accuracy: %s\n' % (i, olpchrsize[i]/chrrefhb[i])
            fo5.write(pfline)
            

if __name__ == '__main__':
    cal_acc2snparr(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
