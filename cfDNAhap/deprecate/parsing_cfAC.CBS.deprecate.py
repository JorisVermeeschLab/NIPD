#!/usr/bin/env python

'''
Author: Huiwen Che
Created: Mar 05 2018
Merging scripts from parsingAERC_heteropat.py, parsingAERC_heteromat.py, parsingAERC_hetboth.py; modifying input and output handling according to modification in previous steps
Use: 
parsing cfDNA pileup (gatk aerc count file) according to different SNPs types
'''

from __future__ import division
import sys
import os
from Helper import mkdir, check_dirslash, basename
from Utility_stats import cal_basic_stats
import numpy as np

def Sum_FAR_Stats(fhtype, cffhapdir):
    pref = basename(fhtype)
    sumstats = '%s%s.summary' % (cffhapdir, pref)

    (tcnt, subA0, subA1, subB0, subB1) = (0, 0, 0, 0, 0)
    (T1chrdistl, T1alldistl, T1prechr, T1prepos, T1dist) = ([], [], '0', 0, 0)
    (T2chrdistl, T2alldistl, T2prechr, T2prepos, T2dist) = ([], [], '0', 0, 0)
    T1distarr = np.empty(0)
    T2distarr = np.empty(0)

    with open(sumstats, 'w+') as fo:
        fo.write('##Summary for SNPs\n')
        with open(fhtype, 'rU') as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    (CHROM, POS, RSID, REF, ALT, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, TCNT) = line.strip().split('\t')
                    tcnt += 1
                    if (TYPE == 'P1' and PH1 == '0') or (TYPE == 'M1' and MH1 == '0'):
                        subA0 += 1
                    elif (TYPE == 'P1' and PH1 == '1') or (TYPE == 'M1' and MH1 == '1'):
                        subA1 += 1
                    elif (TYPE == 'P2' and PH1 == '0') or (TYPE == 'M2' and MH1 == '0'):
                        subB0 += 1
                    elif (TYPE == 'P2' and PH1 == '1') or (TYPE == 'M2' and MH1 == '1'):
                        subB1 += 1
                    
                    
                    if TYPE == 'P1' or TYPE == 'M1':
                        if CHROM != T1prechr:
                            if T1prechr != '0':
                                T1distarr = np.append(T1distarr, T1prechr)
                                T1chrstat = cal_basic_stats(T1chrdistl)
                                T1distarr = np.append(T1distarr, T1chrstat)
                                T1chrdistl = []
                        elif CHROM == T1prechr:
                            T1dist = int(POS) - T1prepos
                            T1chrdistl.append(T1dist)
                            T1alldistl.append(T1dist)
                        T1prepos = int(POS)
                        T1prechr = CHROM
                        
                    elif TYPE == 'P2' or TYPE == 'M2':
                        if CHROM != T2prechr:
                            if T2prechr != '0':
                                T2distarr = np.append(T2distarr, T2prechr)
                                T2chrstat = cal_basic_stats(T2chrdistl)
                                T2distarr = np.append(T2distarr, T2chrstat)
                                T2chrdistl = []
                        elif CHROM == T2prechr:
                            T2dist = int(POS) - T2prepos
                            T2chrdistl.append(T2dist)
                            T2alldistl.append(T2dist)
                        T2prepos = int(POS)
                        T2prechr = CHROM
                        
            #when loop ends, calculate the stats for the last chromosme
            T1distarr = np.append(T1distarr, T1prechr)
            T1chrstat = cal_basic_stats(T1chrdistl)
            T1dim1 = len(T1chrstat) + 1
            T1distarr = np.append(T1distarr, T1chrstat)

            T1dim0 = int(len(T1distarr)/T1dim1)
            T1distarr = np.reshape(T1distarr, (T1dim0, T1dim1))

            T1alldiststat = cal_basic_stats(T1alldistl)
            
            T2distarr = np.append(T2distarr, T2prechr)
            T2chrstat = cal_basic_stats(T2chrdistl)
            T2dim1 = len(T2chrstat) + 1
            T2distarr = np.append(T2distarr, T2chrstat)

            T2dim0 = int(len(T2distarr)/T2dim1)
            T2distarr = np.reshape(T2distarr, (T2dim0, T2dim1))

            T2alldiststat = cal_basic_stats(T2alldistl)

            fo.write('#TotalNoSNPs\tP1_0/M1_0\tP1_1/M1_1\tP2_0/M2_0\tP2_1/M2_1\n')
            
            fo.write('%s\t%s\t%s\t%s\t%s\n' % (tcnt, subA0, subA1, subB0, subB1))
            fo.write('#Distance between P1/M1 SNPs\n#Mean\tS.D.\tMedian\tMax\tMin\noverall\t')
            print >>fo, '\t'.join(map(str, T1alldiststat))
            for i in range(T1dim0):
                print >>fo, '\t'.join(T1distarr[i])
                
            fo.write('#Distance between P2/M2 SNPs\n#Mean\tS.D.\tMedian\tMax\tMin\noverall\t')
            print >>fo, '\t'.join(map(str, T2alldiststat))
            for i in range(T2dim0):
                print >>fo, '\t'.join(T2distarr[i])


def cal_heteropat_FAR(aerc, HeteroPat, mindepth, maxdepth, mediancov, cffhapdir, ffest):
    pref = basename(HeteroPat)
    farfo = '%s%s.aerc.tsv' % (cffhapdir, pref)
    
    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    mediancov = int(mediancov)
    
    cffeval = dict()
    (PH1, PH2, MH1, MH2, TYPE) = ('9', '9', '9', '9', 'unknown')
    FAR = 9
    distlimit = 150
    
    fflimit = float(ffest)*0.3
    SITEinf = 'na'
    ffthresh = float(ffest)*5
    
    multilistT1 = list()
    multilistT2 = list()
    deletelistT1 = list()
    deletelistT2 = list()
    
    with open(HeteroPat, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                val = '%s:%s:%s:%s' % (PH1, PH2, MH1, TYPE)
                cffeval[refID] = val
                
    with open(aerc, 'rU') as fh1:
        for line in fh1:
            if not line.strip().startswith('contig'):
                (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                checkID = '%s:%s' % (contig, position)
                if checkID in cffeval:
                    if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                        (PH1, PH2, MH1, TYPE) = cffeval[checkID].split(':')
                        MH2 = MH1
                        if contig != 'Y':
                            if TYPE == 'P1':
                                if MH1 == '0':
                                    FAR = 2*int(altCount)/int(totalCount)
                                elif MH1 == '1':
                                    FAR = 2*int(refCount)/int(totalCount)
                                    
                                if FAR > fflimit:
                                    SITEinf = 'PH2'
                                else:
                                    SITEinf = 'PH1'
                                    
                                #filter out extreme FAR        
                                if FAR <= ffthresh:
                                    #add tuple to create structure
                                    multilistT1.append((contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, totalCount))
                                    
                            elif TYPE == 'P2':
                                if MH1 == '0':
                                    FAR = 2*int(altCount)/int(totalCount)
                                elif MH1 == '1':
                                    FAR = 2*int(refCount)/int(totalCount)
                                    
                                if FAR > fflimit:
                                    SITEinf = 'PH1'
                                else:
                                    SITEinf = 'PH2'
                                    
                                if FAR <= ffthresh:
                                    multilistT2.append((contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, totalCount))
                                                     
                    elif contig == 'Y' and int(totalCount) <= maxdepth:
                        FAR = 2*int(altCount)/mediancov
                        if FAR <= ffthresh:
                            SITEinf = 'PH1'
                            multilistT2.append((contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, totalCount))
                            
    dtp = np.dtype([('CHROM', 'S10'), ('POS', 'int32'), ('RSID', 'S10'), ('REF', 'S10'), ('ALT', 'S10'), ('PH1', 'S10'), ('PH2', 'S10'), ('MH1', 'S10'), ('MH2', 'S10'), ('FAR', 'float32'), ('SITEINFER', 'S10'), ('TYPE', 'S10'), ('TCNT', 'int32')])
    multiarrT1 = np.array(multilistT1, dtype=dtp)
    multiarrT2 = np.array(multilistT2, dtype=dtp)
    
    #for (chrom in CHROM):
    #    multiarrT1chr = multiarrT1[multiarrT1['CHROM'] == chrom]
    #    multiarrT1chrsub = multiarrT1chr[['POS', 'TCNT']]
        #get index 
        #searchindex = np.where(np.diff(multiarrT1.chr.sub['POS'])<=150)[0]
        #nextindex = searchindex + 1
    for i in range(1, multiarrT1.shape[0]):
        j = i - 1
        while (j in deletelistT1):
            j = j - 1    
        
        dist = multiarrT1[i]['POS']-multiarrT1[j]['POS']
        if (multiarrT1[i]['CHROM'] == multiarrT1[j]['CHROM'] and dist < distlimit):
            if (multiarrT1[i]['TCNT'] > multiarrT1[j]['TCNT']):
                deletelistT1.append(j)
            else:
                deletelistT1.append(i)
                
    multiarrT1 = np.delete(multiarrT1, deletelistT1)
    
    for i in range(1, multiarrT2.shape[0]):
        j = i - 1
        while (j in deletelistT2):
            j = j-1    
        
        dist = multiarrT2[i]['POS']-multiarrT2[j]['POS']
        if (multiarrT2[i]['CHROM'] == multiarrT2[j]['CHROM'] and dist < distlimit):
            if (multiarrT2[i]['TCNT'] > multiarrT2[j]['TCNT']):
                deletelistT2.append(j)
            else:
                deletelistT2.append(i)

    multiarrT2 = np.delete(multiarrT2, deletelistT2)
    
    multiarr = np.append(multiarrT1, multiarrT2)
    multiarr = np.sort(multiarr, order=['CHROM', 'POS'])

    with open(farfo, 'w+') as fo0:
        fo0.write('#CHROM\tPOS\tRSID\tREF\tALT\tPH1\tPH2\tMH1\tMH2\tFAR\tSITEINFER\tTYPE\tTCNT\n')
        for k in range(0, multiarr.shape[0]):
            print >>fo0, '\t'.join(map(str, multiarr[k]))
                                                                       
    Sum_FAR_Stats(farfo, cffhapdir)



def cal_heteromat_FAR(aerc, HeteroMat, mindepth, maxdepth, correctionft, cffhapdir, ffest):
    pref = basename(HeteroMat)
    farfo = '%s%s.aerc.tsv' % (cffhapdir, pref)

    mindepth = int(mindepth)
    maxdepth = int(maxdepth)
    correctionft = float(correctionft)
 
    adjref = 0.5/correctionft
    adjalt = 0.5/(1-correctionft)

    cffeval = dict()
    (PH1, PH2, MH1, MH2, TYPE) = ('9', '9', '9', '9', 'unknown')
    FAR = 9
    
    distlimit = 150
    
    fflimit = float(ffest)*0.3
    SITEinf = 'na'
    ffthresh = float(ffest)*5
    
    multilistT1 = list()
    multilistT2 = list()
    deletelistT1 = list()
    deletelistT2 = list()
    
    with open(HeteroMat, 'rU') as fh0:
        for line in fh0:
            if not line.strip().startswith('#'):
                (CHROM, POS, RSID, REF, ALT, QUAL, PATRAF, MATRAF, SIBRAF, PH1, PH2, MH1, MH2, SH1, SH2, TYPE, PATDP, MATDP, SIBDP) = line.strip().split('\t')
                refID = '%s:%s' % (CHROM, POS)
                val = '%s:%s:%s:%s' % (PH1, MH1, MH2, TYPE)
                cffeval[refID] = val

    with open(aerc, 'rU') as fh1:
        for line in fh1:
            if not line.strip().startswith('contig'):
                (contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases, improperPairs) = line.strip().split('\t')
                checkID = '%s:%s' % (contig, position)
                if checkID in cffeval:
                    if int(totalCount) >= mindepth and int(totalCount) <= maxdepth:
                        (PH1, MH1, MH2, TYPE) = cffeval[checkID].split(':')
                        PH2 = PH1
                        
                        #need improvement, female and male fetus correction difference on X
                        if TYPE == 'M1':
                            if MH1 == '0':
                                FAR = (int(refCount)*adjref-int(altCount)*adjalt)/(int(refCount)*adjref+int(altCount)*adjalt)
                            elif MH1 == '1':
                                FAR = (int(altCount)*adjalt-int(refCount)*adjref)/(int(refCount)*adjref+int(altCount)*adjalt)
                                
                            if FAR > fflimit:
                                SITEinf = 'MH1'
                            else:
                                SITEinf = 'MH2'
                                
                            if abs(FAR) <= ffthresh:
                                
                                multilistT1.append((contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, totalCount))
                            
                        elif TYPE == 'M2':
                            if MH1 == '0':
                                FAR = (int(refCount)*adjref-int(altCount)*adjalt)/(int(refCount)*adjref+int(altCount)*adjalt)
                            elif MH1 == '1':
                                FAR = (int(altCount)*adjalt-int(refCount)*adjref)/(int(refCount)*adjref+int(altCount)*adjalt)
                                
                            if FAR < (0-fflimit):
                                SITEinf = 'MH2'
                            else:
                                SITEinf = 'MH1'
                                
                            if abs(FAR) <= ffthresh:
                                
                                multilistT2.append((contig, position, variantID, refAllele, altAllele, PH1, PH2, MH1, MH2, FAR, SITEinf, TYPE, totalCount))
                                
                            
    dtp = np.dtype([('CHROM', 'S10'), ('POS', 'int32'), ('RSID', 'S10'), ('REF', 'S10'), ('ALT', 'S10'), ('PH1', 'S10'), ('PH2', 'S10'), ('MH1', 'S10'), ('MH2', 'S10'), ('FAR', 'float32'), ('SITEINFER', 'S10'), ('TYPE', 'S10'), ('TCNT', 'int32')])
    multiarrT1 = np.array(multilistT1, dtype=dtp)
    multiarrT2 = np.array(multilistT2, dtype=dtp)
    
    for i in range(1, multiarrT1.shape[0]):
        j = i - 1
        while (j in deletelistT1):
            j = j - 1    
        
        dist = multiarrT1[i]['POS']-multiarrT1[j]['POS']
        if (multiarrT1[i]['CHROM'] == multiarrT1[j]['CHROM'] and dist < distlimit):
            if (multiarrT1[i]['TCNT'] > multiarrT1[j]['TCNT']):
                deletelistT1.append(j)
            else:
                deletelistT1.append(i)
                
    multiarrT1 = np.delete(multiarrT1, deletelistT1)
    
    for i in range(1, multiarrT2.shape[0]):
        j = i - 1
        while (j in deletelistT2):
            j = j - 1    
        
        dist = multiarrT2[i]['POS']-multiarrT2[j]['POS']
        if (multiarrT2[i]['CHROM'] == multiarrT2[j]['CHROM'] and dist < distlimit):
            if (multiarrT2[i]['TCNT'] > multiarrT2[j]['TCNT']):
                deletelistT2.append(j)
            else:
                deletelistT2.append(i)

    multiarrT2 = np.delete(multiarrT2, deletelistT2)
    
    multiarr = np.append(multiarrT1, multiarrT2)
    multiarr = np.sort(multiarr, order=['CHROM', 'POS'])

    with open(farfo, 'w+') as fo0:
        fo0.write('#CHROM\tPOS\tRSID\tREF\tALT\tPH1\tPH2\tMH1\tMH2\tFAR\tSITEINFER\tTYPE\tTCNT\n')
        for k in range(0, multiarr.shape[0]):
            print >>fo0, '\t'.join(map(str, multiarr[k]))
    
    Sum_FAR_Stats(farfo, cffhapdir)



def get_tasktype(tasktype):
    if tasktype.lower() == 'heteropat':
        tasktype = 'heteropat'
    elif tasktype.lower() == 'heteromat':
        tasktype = 'heteromat'
    elif tasktype.lower() == 'heteromatboth':
        tasktype = 'heteromatboth'
    elif tasktype.lower() == 'heteropatboth':
        tasktype = 'heteropatboth'
    else:
        raise ValueError('File not exists.\n')
    return(tasktype)
        



def calAR_main(aerc, tasktype, pref, mindepth, maxdepth, mediancov, correctionft, parentshapdir, cffhapdir, ffest):
    try:
        tasktype = get_tasktype(tasktype)
    except ValueError:
        print 'SNP file not exists. please check filename.'

    if tasktype == 'heteropat':
        infh = parentshapdir + pref + '.heteropat.tsv'
        cal_heteropat_FAR(aerc, infh, mindepth, maxdepth, mediancov, cffhapdir, ffest)
    elif tasktype == 'heteromat':
        infh = parentshapdir + pref + '.heteromat.tsv'
        cal_heteromat_FAR(aerc, infh, mindepth, maxdepth, correctionft, cffhapdir, ffest)
    elif tasktype == 'heteromatboth':
        infh = parentshapdir + pref + '.heteromatboth.tsv'
        cal_heteromat_FAR(aerc, infh, mindepth, maxdepth, correctionft, cffhapdir, ffest)
    elif tasktype == 'heteropatboth':
        infh = parentshapdir + pref + '.heteropatboth.tsv'
        cal_heteropat_FAR(aerc, infh, mindepth, maxdepth, mediancov, cffhapdir, ffest)



if __name__ == '__main__':
    calAR_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
