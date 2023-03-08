#!/usr/bin/env python

import sys
import os

def AB2ACTGconversion(cyto12AB, commonvcfbim, outfh):
    actgdb = dict()

    with open(commonvcfbim, 'r') as fh0:
        for line0 in fh0:
            line0info = line0.strip().split('\t')
            (chrom, rsid, pos, refallele, altallele) = (line0info[0], line0info[1], line0info[3], line0info[4], line0info[5])
            actgval = '%s:%s:%s:%s' % (chrom, pos, refallele, altallele)
            actgdb[rsid] = actgval

    with open(outfh, 'w') as fo:
        fo.write('CHR\tPOS\tRef\tAlt\tRSID\tSNP\tIllumStrd\tSourceStrd\tRefStrd\tOA\tOB\tFA\tFB\n')
        with open(cyto12AB, 'r') as fh1:
            for line1 in fh1:
                if not line1.startswith('Name'):
                    line1info = line1.strip().split('\t')
                    (mrsid, SNP, illumstrand, sourcestrand, refstrand) = (line1info[0], line1info[1], line1info[2], line1info[3], line1info[4])
                    if mrsid in actgdb:
                        (chrom, pos, refallele, altallele) = actgdb[mrsid].split(':')
                        Aallele = 'N'
                        Ballele = 'N'
                        FA = 'N'
                        FB = 'N'
                        if SNP == '[A/G]':
                            if illumstrand == 'TOP':
                                Aallele = 'A'
                                Ballele = 'G'
                            else:
                                Aallele = 'G'
                                Ballele = 'A'
                        elif SNP == '[T/C]':
                            if illumstrand == 'TOP':
                                Aallele = 'C'
                                Ballele = 'T'
                            else:
                                Aallele = 'T'
                                Ballele = 'C'
                        elif SNP == '[A/C]':
                            if illumstrand == 'TOP':
                                Aallele = 'A'
                                Ballele = 'C'
                            else:
                                Aallele = 'C'
                                Ballele = 'A'
                        elif SNP == '[T/G]':
                            if illumstrand == 'TOP':
                                Aallele = 'G'
                                Ballele = 'T'
                            else:
                                Aallele = 'T'
                                Ballele = 'G'

                        if (Aallele == refallele and (Ballele in altallele)):
                            FA = Aallele
                            FB = Ballele
                            altallele = Ballele
                        elif (Ballele == refallele and (Aallele in altallele)):
                            FA = Aallele
                            FB = Ballele
                            altallele = Aallele
                        else:
                            if Aallele == 'A':
                                Aallele = 'T'
                            elif Aallele == 'G':
                                Aallele = 'C'
                            elif Aallele == 'C':
                                Aallele = 'G'
                            elif Aallele == 'T':
                                Aallele = 'A'

                            if Ballele == 'A':
                                Ballele = 'T'
                            elif Ballele == 'G':
                                Ballele = 'C'
                            elif Ballele == 'C':
                                Ballele = 'G'
                            elif Ballele == 'T':
                                Ballele = 'A'

                            if (Aallele == refallele and (Ballele in altallele)):
                                FA = Aallele
                                FB = Ballele
                                altallele = Ballele
                            elif (Ballele == refallele and (Aallele in altallele)):
                                FA = Aallele
                                FB = Ballele
                                altallele = Aallele
                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, pos, refallele, altallele, mrsid, SNP, illumstrand, sourcestrand, refstrand, Aallele, Ballele, FA, FB)
                        fo.write(pfline)

if __name__ == '__main__':
    AB2ACTGconversion(sys.argv[1], sys.argv[2], sys.argv[3])
