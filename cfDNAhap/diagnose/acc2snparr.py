#!/usr/bin/env python

from __future__ import division
import sys
import os
from subprocess import call

def cal_acc2snparr(snparrfh, cffh, workdir):
    snptemp1 = workdir + '/snp1.bed.temp'
    snptemp2 = workdir + '/snp2.bed.temp'
    cftemp1 = workdir + '/cf1.bed.temp'
    cftemp2 = workdir + '/cf2.bed.temp'
    olp1 = workdir + '/olp1.bed'
    olp2 = workdir + '/olp2.bed'
    chrlist = map(str, range(1, 23))
    totrefhl = 0
    chrrefhl = dict()
    
    with open(snptemp1, 'w+') as snp1, open(snptemp2, 'w+') as snp2:
        with open(snparrfh, 'rU') as fh1:
            for line in fh1:
                (chrom, start, end, length, cumlength, HL, infsnps, infsnpsdisc, frac, fitted, resid, prob) = line.strip().split('\t')
                if chrom in chrlist:
                    if int(HL) == 1:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, HL)
                        snp1.write(pfline)
                    elif int(HL) == 2:
                        pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, HL)
                        snp2.write(pfline)
                    blocksize = int(end)-int(start)
                    totrefhl += blocksize
                    if chrom in chrrefhl:
                        chrrefhl[chrom] += blocksize
                    else:
                        chrrefhl[chrom] = blocksize

                    
    with open(cftemp1, 'w+') as cf1, open(cftemp2, 'w+') as cf2:
        with open(cffh, 'rU') as fh2:
            for line in fh2:
                if not line.strip().startswith('#'):
                    (chrom,start, end, blocksize, HL) = line.strip().split('\t')
                    if chrom in chrlist:
                        if HL == "MH1" or HL == "PH1":
                            pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, HL)
                            cf1.write(pfline)
                        elif HL == "MH2" or HL == "PH2":
                            pfline = '%s\t%s\t%s\t%s\n' % (chrom, start, end, HL)
                            cf2.write(pfline)
                            
    snptemp1st = workdir + '/snp1.sorted.bed.temp'
    cftemp1st = workdir + '/cf1.sorted.bed.temp'
    shcmd1 = '/cm/shared/apps/bedtools/current/bin/bedtools sort -i %s > %s' % (snptemp1, snptemp1st)
    shcmd2 = '/cm/shared/apps/bedtools/current/bin/bedtools sort -i %s > %s' % (cftemp1, cftemp1st)
    shcmd12 = '/cm/shared/apps/bedtools/current/bin/bedtools intersect -a %s -b %s > %s' % (snptemp1st, cftemp1st, olp1)
    call(shcmd1, shell=True)
    call(shcmd2, shell=True)
    call(shcmd12, shell=True)
    
    snptemp2st = workdir + '/snptemp2.sorted.bed.temp'
    cftemp2st = workdir + '/cftemp2.sorted.bed.temp'
    shcmd3 = '/cm/shared/apps/bedtools/current/bin/bedtools sort -i %s > %s' % (snptemp2, snptemp2st)
    shcmd4 = '/cm/shared/apps/bedtools/current/bin/bedtools sort -i %s > %s' % (cftemp2, cftemp2st)
    shcmd34 = '/cm/shared/apps/bedtools/current/bin/bedtools intersect -a %s -b %s > %s' % (snptemp2st, cftemp2st, olp2)
    call(shcmd3, shell=True)
    call(shcmd4, shell=True)
    call(shcmd34, shell=True)
    
    shcmd9 = 'rm *.temp'
    call(shcmd9, shell=True)

    cfname = cffh.rsplit('/',1)[1]
    outpref = cfname.rsplit('.', 1)[0]
    sumstat = workdir + '/'+ outpref +'.acc2snp.summary'
    olptotsize = 0
    olpchrsize = dict()

    with open(olp1, 'rU') as fh5:
        for line in fh5:
            (chrom, start, end, val) = line.strip().split('\t')
            olpsegsize = int(end) - int(start)
            olptotsize += olpsegsize
            chrom = str(chrom)
            if chrom in olpchrsize:
                olpchrsize[chrom] += olpsegsize
            else:
                olpchrsize[chrom] = olpsegsize

    with open(olp2, 'rU') as fh6:
        for line in fh6:
            (chrom, start, end, val) = line.strip().split('\t')
            olpsegsize = int(end) - int(start)
            olptotsize += olpsegsize
            chrom = str(chrom)
            if chrom in olpchrsize:
                olpchrsize[chrom] += olpsegsize
            else:
                olpchrsize[chrom] = olpsegsize

    with open(sumstat, 'w+') as fo7:
        pfline = 'total.acc.: %s\n' % (olptotsize/totrefhl)
        fo7.write(pfline)
        for i in olpchrsize:
            pfline = 'chr%s.acc.: %s\n' % (i, olpchrsize[i]/chrrefhl[i])
            fo7.write(pfline)

            
    

if __name__ == '__main__':
    cal_acc2snparr(sys.argv[1], sys.argv[2], sys.argv[3])
