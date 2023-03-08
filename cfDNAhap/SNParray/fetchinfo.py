#!/usr/bin/env python

import sys
import os

def fetchAB(cytocombine, abreport, outfh):
    actgdict = dict()
    with open(cytocombine, 'r') as fh0:
        for line0 in fh0:
            if not line0.startswith('CHR'):
                line0info = line0.strip().split('\t')
                (rsid, ref, alt, FA, FB) = (line0info[4], line0info[2], line0info[3], line0info[11], line0info[12])
                actgval = '%s:%s:%s:%s' % (ref, alt, FA, FB)
                actgdict[rsid] = actgval

    with open(outfh, 'w') as fo:
        fo.write('CHROM\tPOSITION\tRSID\tREF\tALT\tSIBGT\tPATGT\tMATGT\tCSIBGT\tCPATGT\tCMATGT\n')
        with open(abreport, 'r') as fh1:
            for line1 in fh1:
                if not line1.startswith('Name'):
                    line1info = line1.strip().split('\t')
                    (mRSID,Chrom,Position,sibgt,fathergt, mothergt) = (line1info[0], line1info[1], line1info[2], line1info[3], line1info[6], line1info[9])
                    csibgt = cfathergt = cmothergt = 'NN'
                    if mRSID in actgdict:
                        (ref, alt, FA, FB) = actgdict[mRSID].split(':')
                        if ref == FA:
                            if (sibgt == 'AA'):
                                csibgt = '%s/%s' % (0,0)
                            elif (sibgt == 'BB'):
                                csibgt = '%s/%s' % (1,1)
                            elif (sibgt == 'AB'):
                                csibgt = '%s/%s' % (0,1)

                            if (fathergt == 'AA'):
                                cfathergt = '%s/%s' % (0,0)
                            elif (fathergt == 'BB'):
                                cfathergt = '%s/%s' % (1,1)
                            elif (fathergt == 'AB'):
                                cfathergt = '%s/%s' % (0,1)

                            if (mothergt == 'AA'):
                                cmothergt = '%s/%s' % (0,0)
                            elif (mothergt == 'BB'):
                                cmothergt = '%s/%s' % (1,1)
                            elif (mothergt == 'AB'):
                                cmothergt = '%s/%s' % (0,1)

                        elif alt == FA:
                            if (sibgt == 'AA'):
                                csibgt = '%s/%s' % (1,1)
                            elif (sibgt == 'BB'):
                                csibgt = '%s/%s' % (0,0)
                            elif (sibgt == 'AB'):
                                csibgt = '%s/%s' % (0,1)

                            if (fathergt == 'AA'):
                                cfathergt = '%s/%s' % (1,1)
                            elif (fathergt == 'BB'):
                                cfathergt = '%s/%s' % (0,0)
                            elif (fathergt == 'AB'):
                                cfathergt = '%s/%s' % (0,1)

                            if (mothergt == 'AA'):
                                cmothergt = '%s/%s' % (1,1)
                            elif (mothergt == 'BB'):
                                cmothergt = '%s/%s' % (0,0)
                            elif (mothergt == 'AB'):
                                cmothergt = '%s/%s' % (0,1)

                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (Chrom,Position,mRSID,ref,alt,sibgt,fathergt, mothergt,csibgt,cfathergt,cmothergt)
                        fo.write(pfline)

if __name__ == '__main__':
    fetchAB(sys.argv[1], sys.argv[2], sys.argv[3])
