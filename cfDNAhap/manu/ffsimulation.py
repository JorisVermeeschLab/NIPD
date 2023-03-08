#!/usr/bin/env python

import sys
import os
from subprocess import call

def mkdir(sdir):
    if not os.path.isdir(sdir):
        os.makedirs(sdir)

def write_mix_script(scriptdir, workdir, simucode):
    simudir = workdir + '/' + simucode
    mkdir(simudir)
    bamdir = simudir + '/bam'
    aercdir = simudir + '/aerc'
    mkdir(bamdir)
    mkdir(aercdir)

    motherbam = '/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/4_Jun2018/3_align/GC059006_GTCGTAGA_S2/GC059006_GTCGTAGA_S2.sorted.merged.md.filtered.bam'
    childbam = '/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/simulation/ffdilute/mixture/bam/GC063925_TTCACGCA_S3.sorted.merged.md.filtered.bam'

    motherper = [85, 875, 90, 92, 93, 95, 956, 963, 97, 978, 99, 100]
    childper = ['137', '12', '10', '09', '08', '07', '065', '06', '055', '05', '045', '04']
    ffest = [19.3, 17, 14.2, 12.8, 11.4, 10, 9.28, 8.5, 7.86, 7.1, 6.4, 5.5]
    pername = ['85perADD13d7per', '875perADD12per','90perADD10per','92perADD9per','93perADD8per','95perADD7per','956perADD6d5per','96perADD6per','97perADD5d5per','97perADD5per','99perADD4d5per','100perADD4per']

    for prop in pername:
        subdir = simudir + '/' + prop
        mkdir(subdir)
        shcmd = 'cp -r /uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/simulation/ffdilute/mixture/%s/5_variant/ %s' % (prop, subdir)
        call(shcmd, shell=True)

    simucode = int(simucode)
    #mixscript = '%s/mix.%s.sh' % (scriptdir, simucode)

    for i in range(0, len(motherper)):
        mixscript = '%s/mix.%s.%s.sh' % (scriptdir, simucode, pername[i])
        with open(mixscript, 'w+') as fo:
            motherinper = '%s.%s' % (simucode, motherper[i])
            childinper = '%s.%s' % (simucode, childper[i])

            motheroutbam = bamdir + '/mother.' + pername[i] + '.bam'
            childoutbam = bamdir + '/child.' + pername[i] + '.bam'

            if (motherper[i] == 100):
                motheroutbam = motherbam
            else:
                pfline1 = '/cm/shared/apps/samtools/1.3.1/samtools view -bS -s %s %s > %s\n' % (motherinper, motherbam, motheroutbam)
                fo.write(pfline1)

            pfline2 = '/cm/shared/apps/samtools/1.3.1/samtools view -bS -s %s %s > %s\n' % (childinper, childbam, childoutbam)
            fo.write(pfline2)

            mixbam = bamdir + '/mix.' + pername[i] + '.bam' 
            pfline3 = '/cm/shared/apps/jdk/1.8.0/bin/java -Xmx16G -jar /cm/shared/apps/picard/2.9.0/picard.jar MergeSamFiles TMP_DIR=/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/9_mixOct2018/tmp I=%s I=%s O=%s VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ASSUME_SORTED=true\n' % (motheroutbam, childoutbam, mixbam)
            fo.write(pfline3)

            pfline4 = '/cm/shared/apps/samtools/1.3.1/samtools index %s\n' % (mixbam)
            fo.write(pfline4)
            
            aercout = aercdir + '/mix.' + pername[i] + '.aerc.csv'
            pfline5 = '/cm/shared/apps/jdk/1.8.0/bin/java -Xmx4G -jar /cm/shared/apps/gatk/3.7/GenomeAnalysisTK.jar -T ASEReadCounter -R /home/hche0/00_database/GRCh37/hs37d5/hs37d5.fa -I %s -sites /uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/uz/nipt_research/Huiwen/09_hap/4_Jun2018/5_variant/0_rawvcf/GC059006_parental_sib.joint.raw.snps.indels.vcf -o %s --countOverlapReadsType COUNT_FRAGMENTS_REQUIRE_SAME_BASE -U ALLOW_N_CIGAR_READS -minDepth 2 --minMappingQuality 20 --minBaseQuality 2\n' % (mixbam, aercout)
            fo.write(pfline5)

if __name__ == '__main__':
    write_mix_script(sys.argv[1], sys.argv[2], sys.argv[3])
