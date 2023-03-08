#!/usr/bin/env python

import sys
import os
import ConfigParser

curdir = os.getcwd()

def mkdir(sdir):
    if not os.path.isdir(sdir):
        os.makedirs(sdir)


def check_fexist(f):
    if os.path.isfile(f) == True:
        return(1)
    else:
        return(0)


def check_dirslash(d):
    if d.strip()[-1] == '/':
        return(d)
    else:
        d = d + '/'
        return(d)


def get_config():
    parser = ConfigParser.RawConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')
    
    fastqdir = parser.get('mydir', 'fastqdir')
    workdir = parser.get('mydir', 'workdir')
    samtools = parser.get('mytool', 'samtools')
    bwa = parser.get('mytool', 'bwa')
    java = parser.get('mytool', 'java')
    picard = parser.get('mytool', 'picard')
    fastqc = parser.get('mytool', 'fastqc')
    gatk = parser.get('mytool', 'gatk')
    reference = parser.get('myfile', 'reference')
    target = parser.get('myfile', 'target')
    targetflank50 = parser.get('myfile', 'targetflank50')
    targetflank100 = parser.get('myfile', 'targetflank100')
    targetSNPref = parser.get('myfile', 'targetSNPref')
    dbsnpref = parser.get('myfile', 'dbsnpref')
    pedfile = parser.get('myfile', 'pedfile')

    fastqdir = check_dirslash(fastqdir)
    workdir = check_dirslash(workdir)    

    return(fastqdir, workdir, samtools, bwa, java, picard, fastqc, gatk, reference, target, targetflank50, targetflank100, targetSNPref, dbsnpref, pedfile)

    
def get_sample():
    parser = ConfigParser.RawConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')
    
    samplelist = parser.get('myfile', 'samplelist')
    sampleID = []

    with open(samplelist, 'rU') as fh:
        for line in fh:
            sampleID.append(line.rstrip('\n'))
    fh.close()

    return(sampleID)
