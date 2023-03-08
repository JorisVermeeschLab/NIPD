#!/usr/bin/env python

import sys
import os

def array2vcf(convertfh, outfh):
    with open(outfh, 'w') as fo:
        fo.write('##fileformat=VCFv4.2\n')
        fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIBGT\tPATGT\tMATGT\n')
        with open(convertfh, 'r') as fh:
            for line in fh:
                if not line.startswith('CHROM'):
                    PPATGT = PMATGT = PSIBGT = './.'
                    (CHROM,POSITION,RSID,REF,ALT,SIBGT,PATGT,MATGT,CSIBGT,CPATGT,CMATGT) = line.strip().split('\t')
                    if CHROM != 'XY' and CHROM !='Y':
                        if (CPATGT == '0/0' and CMATGT == '1/1'):
                            PPATGT = '0|0'
                            PMATGT = '1|1'
                            FPATGT = '%s:%s,%s:%s' % (PPATGT,100,0,100)
                            FMATGT = '%s:%s,%s:%s' % (PMATGT,0,100,100)

                            if (CSIBGT == '0/1'):
                                PSIBGT = '0|1'
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == 'NC'):
                                PSIBGT = './.'
                                FSIBGT = '%s' % (PSIBGT)
                            elif (CSIBGT == '1/1'):
                                PSIBGT = CSIBGT
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == '0/0'):
                                PSIBGT = CSIBGT
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)
                        elif (CPATGT == '1/1' and CMATGT == '0/0'):
                            PPATGT = '1|1'
                            PMATGT = '0|0'
                            FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                            FMATGT = '%s:%s,%s:%s' % (PMATGT,100,0,100)
                            if (CSIBGT == '0/1'):
                                PSIBGT = '1|0'
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == 'NC'):
                                PSIBGT = './.'
                                FSIBGT = '%s' % (PSIBGT)
                            elif (CSIBGT == '1/1'):
                                PSIBGT = CSIBGT
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == '0/0'):
                                PSIBGT = CSIBGT
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)

                        elif (CPATGT == '0/1' and CMATGT == '0/0'):
                            if (CSIBGT == '0/0'):
                                PPATGT = '0|1'
                                PMATGT = '0|0'
                                PSIBGT = '0|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,100,0,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            elif (CSIBGT == '0/1'):
                                PPATGT = '1|0'
                                PMATGT = '0|0'
                                PSIBGT = '1|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,100,0,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == 'NN'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,100,0,100)
                                FSIBGT = './.'
                            else:
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,0,100,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)
                        elif (CPATGT == '0/1' and CMATGT == '1/1'):
                            if (CSIBGT == '0/1'):
                                PPATGT = '0|1'
                                PMATGT = '1|1'
                                PSIBGT = '0|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == '1/1'):
                                PPATGT = '1|0'
                                PMATGT = '1|1'
                                PSIBGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == 'NN'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = './.'
                            else:
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,100,0,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)

                        elif (CPATGT == '0/0' and CMATGT == '0/1'):
                            if (CSIBGT == '0/0'):
                                PPATGT = '0|0'
                                PMATGT = '0|1'
                                PSIBGT = '0|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,100,0,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            elif (CSIBGT == '0/1'):
                                PPATGT = '0|0'
                                PMATGT = '1|0'
                                PSIBGT = '0|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,100,0,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == 'NN'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = './.'
                            else:
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,100,0,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,0,100,100)

                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)
                        elif (CPATGT == '1/1' and CMATGT == '0/1'):
                            if (CSIBGT == '0/1'):
                                PPATGT = '1|1'
                                PMATGT = '0|1'
                                PSIBGT = '1|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,50,50,100)
                            elif (CSIBGT == '1/1'):
                                PPATGT = '1|1'
                                PMATGT = '1|0'
                                PSIBGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == 'NN'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = './.'
                            else:
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,100,0,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)

                        elif (CPATGT == '0/1' and CMATGT == '0/1'):
                            if (CSIBGT == '0/0'):
                                PPATGT = '0|1'
                                PMATGT = '0|1'
                                PSIBGT = '0|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            elif (CSIBGT == '1/1'):
                                PPATGT = '1|0'
                                PMATGT = '1|0'
                                PSIBGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == '0/1'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,50,50,100)
                            elif (CSIBGT == 'NN'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,50,50,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,50,50,100)
                                FSIBGT = './.'
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)

                        elif (CPATGT == '1/1' and CMATGT == '1/1'):
                            if (CSIBGT == '1/1'):
                                PPATGT = '1|1'
                                PMATGT = '1|1'
                                PSIBGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            elif (CSIBGT == 'NN'):
                                PPATGT = '1|1'
                                PMATGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (PMATGT,0,100,100)
                                FSIBGT = './.'
                            elif (CSIBGT == '0/0'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,100,0,100)
                            elif (CSIBGT == '0/1'):
                                FPATGT = '%s:%s,%s:%s' % (CPATGT,0,100,100)
                                FMATGT = '%s:%s,%s:%s' % (CMATGT,0,100,100)
                                FSIBGT = '%s:%s,%s:%s' % (CSIBGT,50,50,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)

                    elif (CHROM == 'Y'):
                        if (CPATGT == '1/1' and CMATGT == 'NN'):
                            if (CSIBGT == 'NN'):
                                PPATGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = './.'
                                FSIBGT = './.'
                            elif (CIBGT == '1/1'):
                                PPATGT = '1|1'
                                PSIBGT = '1|1'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,0,100,100)
                                FMATGT = './.'
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,0,100,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)
                        elif (CPATGT == '0/0' and CMATGT == 'NN'):
                            if (CSIBGT == 'NN'):
                                PPATGT = '0|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,100,0,100)
                                FMATGT = './.'
                                FSIBGT = './.'
                            elif (CSIBGT == '0/0'):
                                PPATGT = '0|0'
                                PSIBGT = '0|0'
                                FPATGT = '%s:%s,%s:%s' % (PPATGT,100,0,100)
                                FMATGT = './.'
                                FSIBGT = '%s:%s,%s:%s' % (PSIBGT,100,0,100)
                            pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM,POSITION,RSID,REF,ALT,'.','.','.','GT:AD:DP',FSIBGT,FPATGT,FMATGT)
                            fo.write(pfline)


if __name__ == '__main__':
    array2vcf(sys.argv[1], sys.argv[2])
