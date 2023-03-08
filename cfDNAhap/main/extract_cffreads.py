#!/usr/bin/env python

'''
Author: Huiwen Che
Created: December 30 2017
Use: parsing gatk pileup txt (pile up cfDNA sites where parents are of different homozygous alleles), extracting fetal specific reads names and shared reads names
'''
import sys
import re
from Helper import basename

def get_gpile_homodiff_readnames(homodiffpile, piledir):
    pref = basename(homodiffpile)
    fetalreads = '%s%s.fetalreads.txt' % (piledir, pref)
    sharedreads = '%s%s.sharedreads.txt' % (piledir, pref)
    errreads = '%s%s.err.reads.txt' % (piledir, pref)
    
    with open(homodiffpile, 'rU') as fh:
        with open(fetalreads, 'w+') as fo0:
            with open(sharedreads, 'w+') as fo1:
                with open(errreads, 'w+') as fo2:
                    
                    for line in fh:
                        if not line.strip().startswith('#') and not line.strip().startswith('['):
                            #get parternal genotype
                            pat = line.strip().split('\t')[1]
                            pat_gt = pat.split('/')[0]
                            findallele = re.search(r'alleles=\[(\w+)\*,\s(\w+)\]', line)
                            try:
                                if findallele:
                                    refallele = findallele.group(1)
                                    altallele = findallele.group(2)
                                    pileupstr = line.strip().split(' ',5)[3]
                                    readsstr = line.strip().rsplit(' ',1)[1]
                                    pileupl = list(pileupstr)
                                    readsl = readsstr.split(',')

                                    if len(pileupl) == len(readsl):
                                        for i in range(0,len(pileupl)):
                                            readsname = readsl[i].split('@', 1)[0]
                                            if (pileupl[i] == refallele and pat_gt == '0') or (pileupl[i] == altallele and pat_gt == '1'):
                                                print >> fo0, readsname
                                            elif (pileupl[i] == refallele and pat_gt == '1') or (pileupl[i] == altallele and pat_gt == '0'):
                                                print >> fo1, readsname
                                            else:
                                                print >> fo2, readsname
                                    else:
                                        print >> fo2, line
                            except:
                                print >> fo2, 'Allele not found'
                            

def get_pileup_heteropat_readnames(heteropatpile):
    fo_prefix = heteropatpile.rsplit('.', 1)[0]
    fetalreads = '%s.fetalreads.txt' % (fo_prefix)
    sharedreads = '%s.sharedreads.txt' % (fo_prefix)
    errreads = '%s.err.reads.txt' % (fo_prefix)
    
    with open(heteropatpile, 'rU') as fh:
        with open(fetalreads, 'w+') as fo0:
            with open(sharedreads, 'w+') as fo1:
                with open(errreads, 'w+') as fo2:
                    
                    for line in fh:
                        if not line.strip().startswith('#') and not line.strip().startswith('\['):
                            #get maternal genotype
                            mat = line.strip().split('\t')[2]
                            mat_gt = mat.split('/')[0]
                            findallele = re.search(r'alleles=\[(\w+)\*,\s(\w+)\]', line)
                            try:
                                if findallele:
                                    refallele = findallele.group(1)
                                    altallele = findallele.group(2)
                                    pileupstr = line.strip().split(' ',5)[3]
                                    readsstr = line.strip().rsplit(' ',1)[1]
                                    pileupl = list(pileupstr)
                                    readsl = readsstr.split(',')

                                    if len(pileupl) == len(readsl):
                                        for i in range(0,len(pileupl)):
                                            readsname = readsl[i].split('@', 1)[0]
                                            if (pileupl[i] == refallele and mat_gt == '0') or (pileupl[i] == altallele and mat_gt == '1'):
                                                print >> fo1, readsname
                                            elif (pileupl[i] == refallele and mat_gt == '1') or (pileupl[i] == altallele and mat_gt == '0'):
                                                print >> fo0, readsname
                                            else:
                                                print >> fo2, readsname
                                    else:
                                        print >> fo2, line
                            except:
                                print >> fo2, 'Allele not found'
                            
    
    
    

if __name__ == '__main__':
    get_gpile_homodiff_readnames(sys.argv[1], sys.argv[2])
