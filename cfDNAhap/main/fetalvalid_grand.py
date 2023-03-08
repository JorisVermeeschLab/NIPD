#!/usr/bin/env python

from __future__ import division
import sys
from Helper import mkdir, check_dirslash
from Utility_stats import cal_basic_stats
import re
import numpy as np

def Parsing_grand_s1(phasedgrand, phasedtrio, grandfatherID, grandmotherID, fatherID, motherID, fetalID, mindepth, maxdepth, workdir):
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    grandtrio = workdir + '/' + fatherID + '.grand_trio.phasing.tsv'
    fetushap = workdir + '/' + fetalID + '.fetal.genome.tsv'
    pathap = dict()
    
    with open(grandtrio, 'w+') as fo1:
        with open(phasedgrand, 'rU') as fh1:
            for line in fh1:
                if not line.strip().startswith('##'):
                    if line.strip().startswith('#CHROM'):
                        infoline = line.strip().split('\t')
                        fatherpos = int(infoline.index(fatherID))
                        grandfatherpos = int(infoline.index(grandfatherID))
                        grandmotherpos = int(infoline.index(grandmotherID))
                    else:
                        spline = line.strip().split('\t')
                        (CHROM, POS, ID, REF, ALT, QUAL, FORMAT, grandfather, grandmother, father) = (spline[0], spline[1], spline[2], spline[3], spline[4], spline[5], spline[8], spline[grandfatherpos], spline[grandmotherpos], spline[fatherpos])
                    
                        (GPH1, GPH2, GMH1, GMH2, PH1, PH2) = ('.', '.', '.', '.', '.', '.')
                        (grandpat_dp, grandmat_dp, pat_dp) = ('.', '.', '.')
                    
                        if len(REF) == 1 and len(ALT) == 1:
                            dp_search = re.search(r'DP', FORMAT)
                        
                            if dp_search and CHROM != 'Y':
                                FORMATl = FORMAT.split(':')
                                dp_pos = int(FORMATl.index('DP'))
                            
                                if (grandfather != './.') and (grandmother != './.') and (father != './.'):
                                    grandfather_dp = grandfather.split(':')[dp_pos]
                                    grandmother_dp = grandmother.split(':')[dp_pos]
                                    father_dp = father.split(':')[dp_pos]
                                    grandfather_gt = grandfather.split(':')[0]
                                    grandmother_gt = grandmother.split(':')[0]
                                    father_gt = father.split(':')[0]
                                    grandfphaseval = list(grandfather)[1]
                                    grandmphaseval = list(grandmother)[1]
                                    fatherphaseval = list(father)[1]

                                
                                    if grandfather_dp != '.' and grandmother_dp != '.' and father_dp != '.' and int(grandfather_dp) >= mindepth and int(grandfather_dp) <= maxdepth and int(grandmother_dp) >= mindepth and int(grandmother_dp) <= maxdepth and int(father_dp) >= mindepth and int(father_dp) <= maxdepth and grandfphaseval == '|' and grandmphaseval == '|' and fatherphaseval == '|':
                                        (GPH1, GPH2) = grandfather_gt.split('|')
                                        (GMH1, GMH2) = grandmother_gt.split('|')
                                        (PH1, PH2) = father_gt.split('|')
                                        posindex = '%s:%s' % (CHROM, POS)
                                        pathap[posindex] = father_gt
                                        #mathap[posindex] = mother_gt
                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, GPH1, GPH2, GMH1, GMH2, PH1, PH2)
                                        fo1.write(pfline)
                                        
                                        
    
                                    
    with open(fetushap, 'w+') as fo2:
        with open(phasedtrio, 'rU') as fh2:
            for line in fh2:
                if not line.strip().startswith('##'):
                    if line.strip().startswith('#CHROM'):
                        infoline = line.strip().split('\t')
                        fetalpos = int(infoline.index(fetalID))
                        fatherpos = int(infoline.index(fatherID))
                        motherpos = int(infoline.index(motherID))
                    else:
                        spline = line.strip().split('\t')
                        (CHROM, POS, ID, REF, ALT, QUAL, FORMAT, father, mother, fetal) = (spline[0], spline[1], spline[2], spline[3], spline[4], spline[5], spline[8], spline[fatherpos], spline[motherpos], spline[fetalpos])
                    
                        (PH1, PH2, MH1, MH2, FT1, FT2) = ('.', '.', '.', '.', '.', '.')
                        (PATTYPE, MATTYPE, pat_dp, mat_dp, fetal_dp) = ('.', '.', '.', '.', '.')
                    
                        if len(REF) == 1 and len(ALT) == 1:
                            dp_search = re.search(r'DP', FORMAT)
                        
                            if dp_search and CHROM != 'Y':
                                FORMATl = FORMAT.split(':')
                                dp_pos = int(FORMATl.index('DP'))
                            
                                if (father != './.') and (mother != './.') and (fetal != './.'):
                                    father_dp = father.split(':')[dp_pos]
                                    mother_dp = mother.split(':')[dp_pos]
                                    fetal_dp = fetal.split(':')[dp_pos]
                                    father_gt = father.split(':')[0]
                                    mother_gt = mother.split(':')[0]
                                    fetal_gt = fetal.split(':')[0]
                                    fphaseval = list(father)[1]
                                    mphaseval = list(mother)[1]
                                    ftphaseval = list(fetal)[1]
                                
                                    if father_dp != '.' and mother_dp != '.' and fetal_dp != '.' and int(father_dp) >= mindepth and int(father_dp) <= maxdepth and int(mother_dp) >= mindepth and int(mother_dp) <= maxdepth and int(fetal_dp) >= mindepth and int(fetal_dp) <= maxdepth and fphaseval == '|' and mphaseval == '|' and ftphaseval == '|':
                                        (CPH1, CPH2) = father_gt.split('|')
                                        (CMH1, CMH2) = mother_gt.split('|')
                                        (FT1, FT2) = fetal_gt.split('|')

                                        checkid = '%s:%s' % (CHROM, POS)
                                        if checkid in pathap:
                                            (PH1, PH2) = pathap[checkid].split('|')
                                            if (PH1 != PH2):
                                                if (FT1 == PH1):
                                                    PATTYPE = 'PH1'
                                                elif (FT1 == PH2):
                                                    PATTYPE = 'PH2'
                                            else:
                                                PATTYPE = 'Undetermined'
                                            
                                            #(MH1, MH2) = mathap[checkid].split('|')
                                            #if (MH1 != MH2):
                                            #    if (FT2 == MH1):
                                            #        MATTYPE = 'MH1'
                                            #    elif (FT2 == MH2):
                                            #        MATTYPE = 'MH2'
                                            #else:
                                            #    MATTYPE = 'Undetermined'
                                        pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (CHROM, POS, ID, REF, ALT, CPH1, CPH2, CMH1, CMH2, FT1, FT2, PATTYPE, '.')
                                        fo2.write(pfline)
            
    
if __name__ == '__main__':
    Parsing_grand_s1(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
