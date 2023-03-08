#!/usr/bin/env python

import sys
import os
import time
from helper import mkdir, get_config, get_sample, check_fexist

(fastqdir, workdir, samtools, bwa, java, picard, fastqc, gatk, reference, target, targetflank50, targetflank100, targetSNPref, dbsnpref, pedfile) = get_config()
SampleIDs = get_sample()
mytime = time.strftime('%y_%m_%d_%H_%M')
scriptdir = workdir + 'scripts_' + mytime
mkdir(scriptdir)
logdir = workdir + 'log'
mkdir(logdir)


def write_sub_script_qsub(taskscriptdir, taskname, taskworkdir, stderrdir):
    '''
    mpinum added in case of parallel jobs needed
    taskname added in case of multiple qsub in one scriptdir
    stderrdir added in case output to terminal or stdout
    '''
    qsub = 'qsub'

    task_qsub = taskscriptdir + '/' + taskname + '.qsub.sh'
    with open(task_qsub, 'w+') as fo:
        fo.write('#!/bin/bash -l\n\n')

        for script in os.listdir(taskscriptdir):
            tasksearch = taskname + '.task.pbs'
            if script.endswith(tasksearch):
                #myqsub = '%s -d %s -e %s -o %s %s\n' % (qsub, taskworkdir, stderrdir, logdir, script)
                myqsub = '%s %s\n' % (qsub, script)
                fo.write(myqsub)
    return(task_qsub)


def write_sub_script_fastqc():
    fastqc_scriptdir = scriptdir + '/fastqc'
    fastqc_workdir = workdir + '1_qc'
    mkdir(fastqc_scriptdir)
    mkdir(fastqc_workdir)
    
    for sampleid in SampleIDs:

        samplescript = fastqc_scriptdir + '/' + sampleid + '.fastqc.task.pbs'
        with open(samplescript, 'w+') as fo:
            fo.write('#!/bin/bash -l\n\n')
            fo.write('#PBS -l walltime=01:00:00\n')
            fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo.write('#PBS -l mem=64gb\n')
            fo.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo.write('#PBS -m abe\n')
            fo.write('#PBS -N fastqc_%s\n' % (sampleid))
            fo.write('#PBS -A lp_biogenomics\n\n')

            for fastqfile in os.listdir(fastqdir):
                if fastqfile.startswith(sampleid) and fastqfile.endswith('fastq.gz'):
                    qc_task = '%s %s%s -o %s\n' % (fastqc, fastqdir, fastqfile, fastqc_workdir)
                    fo.write(qc_task)
            
    write_sub_script_qsub(fastqc_scriptdir, 'fastqc', fastqc_workdir, logdir)


'''
change UMI to seqeunce name comment position, in order to pass the UMI to bam tag
'''
def write_sub_script_rename():
    UMI_scriptdir = scriptdir + '/UMIrename'
    UMI_workdir = workdir + '2_UMImod'
    mkdir(UMI_scriptdir)
    mkdir(UMI_workdir)
    
    for sampleid in SampleIDs:
        
        samplescript = UMI_scriptdir + '/' + sampleid + '.umimodname.task.pbs'
        UMI_sub_workdir = UMI_workdir + '/' + sampleid
        mkdir(UMI_sub_workdir)
        with open(samplescript, 'w+') as fo:
            fo.write('#!/bin/bash -l\n\n')
            fo.write('#PBS -l walltime=05:00:00\n')
            fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo.write('#PBS -l mem=64gb\n')
            fo.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo.write('#PBS -m abe\n')
            fo.write('#PBS -N UMIrename_%s\n' % (sampleid))
            fo.write('#PBS -A lp_biogenomics\n\n')

            for fastqfile in os.listdir(fastqdir):
                if fastqfile.startswith(sampleid) and fastqfile.endswith('fastq.gz'):
                   fqprefix = fastqfile.split('.', 1)[0]
                   fastqumifile = '%s.rename.fastq.gz' % (fqprefix) 
                   umimod_task = 'zcat %s%s | sed \'s/\\(^@.*:\\([A-Z]\\{10\\}\\)\\)\\s\\(.*\\)/\\1 RX:Z:\\2/g\' | gzip -c > %s/%s\n' % (fastqdir, fastqfile, UMI_sub_workdir, fastqumifile)
                   fo.write(umimod_task)

    write_sub_script_qsub(UMI_scriptdir, 'umimodname', UMI_workdir, logdir)


def write_sub_script_align():
    align_scriptdir = scriptdir + '/align'
    align_workdir = workdir + '3_align'
    mkdir(align_scriptdir)
    mkdir(align_workdir)
    tmpdir = workdir + 'tmp'
    mkdir(tmpdir)

    for sampleid in SampleIDs:

        samplescript = align_scriptdir + '/' + sampleid + '.align.task.pbs'
        align_sub_workdir = align_workdir + '/' + sampleid
        mkdir(align_sub_workdir)
        lanebams = []
        rmlanebams = []
        with open(samplescript, 'w+') as fo:
            fo.write('#!/bin/bash -l\n\n')
            fo.write('#PBS -l walltime=05:00:00\n')
            fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo.write('#PBS -l mem=64gb\n')
            fo.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo.write('#PBS -m abe\n')
            fo.write('#PBS -N align_%s\n' % (sampleid))
            fo.write('#PBS -A lp_biogenomics\n\n')
            fo.write('source switch_to_2015a\n')
            fo.write('module load SAMtools/1.7-foss-2015a\n')
            fo.write('module load GATK/3.7\n\n')
            
            subfastqumidir = workdir + '2_UMImod/' + sampleid
            #lanenum = len(os.listdir(subfastqumidir))/2
            lanenum = 4 #temporary; may need to put into ConfigFile
            #align reads by lane
            for lane in range(1, lanenum+1):
                readgroup = '\'@RG\\tID:%s_L00%s\\tSM:%s\\tPL:Illumina_nextseq\\\'' % (sampleid, lane, sampleid)
                laneread1 = '%s/%s_L00%s_R1_001.rename.fastq.gz' % (subfastqumidir, sampleid, lane)
                laneread2 = '%s/%s_L00%s_R2_001.rename.fastq.gz' % (subfastqumidir, sampleid, lane)
                bamfile = '%s_L00%s.sorted.bam' % (sampleid, lane)
                bwamem_task = '%s mem -R %s -t 12 -M -C %s %s %s | %s sort -@ 12 -O bam -o %s/%s\n' % (bwa, readgroup, reference, laneread1, laneread2, samtools, align_sub_workdir, bamfile)
                fo.write(bwamem_task)
                
                inlanebam = 'I=' + align_sub_workdir + '/' + bamfile
                lanebams.append(inlanebam)
                rmlanebam = align_sub_workdir + '/' + bamfile
                rmlanebams.append(rmlanebam)
            
            #merge all lane bams
            concatbams = ' '.join(lanebams)
            mergebam = '%s.sorted.merged.bam' % (sampleid)
            merge_task = '%s -Xmx16G -jar %s MergeSamFiles TMP_DIR=%s %s O=%s/%s VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ASSUME_SORTED=true\n\n' % (java, picard, tmpdir, concatbams, align_sub_workdir, mergebam)
            fo.write(merge_task)

            '''add function here to check if merge bam succeded, if did, then remove lanebams'''

            #validate merged bam
            validatebam = '%s.sorted.merged.bam.validation' % (sampleid)
            valid_task = '%s -Xmx4G -jar %s ValidateSamFile TMP_DIR=%s I=%s/%s O=%s/%s MODE=SUMMARY\n\n' % (java, picard, tmpdir, align_sub_workdir, mergebam, align_sub_workdir, validatebam)
            fo.write(valid_task)

            #remove lane bams
            rmbams = ' '.join(rmlanebams)
            rm_task = 'rm %s\n\n' % (rmbams) 
            fo.write(rm_task)

            #markduplicate with barcode tag            
            mdbam = '%s.sorted.merged.md.bam' % (sampleid)
            metricfile = '%s.sorted.merged.md.metrics' % (sampleid)
            markdup_task = '%s -Xmx16G -jar %s MarkDuplicates TMP_DIR=%s I=%s/%s METRICS_FILE=%s/%s O=%s/%s BARCODE_TAG=RX ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT\n\n' % (java, picard, tmpdir, align_sub_workdir, mergebam, align_sub_workdir, metricfile, align_sub_workdir, mdbam)
            fo.write(markdup_task)

            #generate md bam stats
            statsfile = '%s.stats' % (mdbam)
            stats_task = '%s stats %s/%s > %s/%s\n\n' % (samtools, align_sub_workdir, mdbam, align_sub_workdir, statsfile)
            fo.write(stats_task)

            #filter md bam
            filteredbam = '%s.sorted.merged.md.filtered.bam' % (sampleid)
            filterstats = '%s.stats' % (sampleid)
            filter_task = '%s view -bS -@ 12 -F 4 -F 256 -F 1024 -q 20 %s/%s > %s/%s\n\n' % (samtools, align_sub_workdir, mdbam, align_sub_workdir, filteredbam)
            filterstats_task = '%s stats %s/%s > %s/%s\n\n' % (samtools, align_sub_workdir, filteredbam, align_sub_workdir, filterstats)
            createbai_task = '%s index %s/%s\n\n' % (samtools, align_sub_workdir, filteredbam)
            fo.write(filter_task)
            fo.write(filterstats_task)
            fo.write(createbai_task)

            #remove raw bam
            rmmergebam_task = 'rm %s/%s\n' % (align_sub_workdir, mergebam)
            fo.write(rmmergebam_task)
   
    write_sub_script_qsub(align_scriptdir, 'align', align_workdir, logdir)


def write_sub_script_coverage():
    cov_scriptdir = scriptdir + '/cov'
    cov_workdir = workdir + '4_cov'
    mkdir(cov_scriptdir)
    mkdir(cov_workdir)
    tmpdir = workdir + 'tmp'
    mkdir(tmpdir)

    for sampleid in SampleIDs:

        samplescript = cov_scriptdir + '/' + sampleid + '.cov.task.pbs'
        samplescript2 = cov_scriptdir + '/' + sampleid + '.readcount.task.pbs'
        cov_sub_workdir = cov_workdir + '/' + sampleid
        mkdir(cov_sub_workdir)

        with open(samplescript, 'w+') as fo:
            fo.write('#!/bin/bash -l\n\n')
            fo.write('#PBS -l walltime=05:00:00\n')
            fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo.write('#PBS -l mem=64gb\n')
            fo.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo.write('#PBS -m abe\n')
            fo.write('#PBS -N fragcov_%s\n' % (sampleid))
            fo.write('#PBS -A lp_biogenomics\n\n')
            fo.write('source switch_to_2015a\n')
            fo.write('module load SAMtools/1.7-foss-2015a\n')
            fo.write('module load GATK/3.7\n\n')
            
            #each step based on previous step; maybe wrap previous step into this step
            align_workdir = workdir + '3_align'
            infilteredbam = '%s/%s/%s.sorted.merged.md.filtered.bam' % (align_workdir, sampleid, sampleid)
            #count coverage in terms of fragments (overlapped pair-end reads)
            fragcov = '%s.sorted.merged.md.filtered.fragcov' % (sampleid)
            fragcov_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T DepthOfCoverage -R %s -o %s/%s -I %s -L %s --countType COUNT_FRAGMENTS\n\n' % (java, tmpdir, gatk, reference, cov_sub_workdir, fragcov, infilteredbam, target)                   
            fo.write(fragcov_task)
            #count coverage in terms of reads
            #readcov = '%s.sorted.merged.md.filtered.readcov' % (sampleid)
            #readcov_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T DepthOfCoverage -R %s -o %s/%s -I %s -L %s --countType COUNT_READS\n\n' % (java, tmpdir, gatk, reference, cov_sub_workdir, readcov, infilteredbam, target)
            #fo.write(readcov_task)


        ##output of CountReads is problematic, needs to be solved
        with open(samplescript2, 'w+') as fo2:
            fo2.write('#!/bin/bash -l\n\n')
            fo2.write('#PBS -l walltime=01:00:00\n')
            fo2.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo2.write('#PBS -l mem=64gb\n')
            fo2.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo2.write('#PBS -m abe\n')
            fo2.write('#PBS -N readcount_%s\n' % (sampleid))
            fo2.write('#PBS -A lp_biogenomics\n\n')
            fo2.write('source switch_to_2015a\n')
            fo2.write('module load SAMtools/1.7-foss-2015a\n')
            fo2.write('module load GATK/3.7\n\n')
                                    
            #count total filtered reads
            #CountReads write to stdout
            #readcountall = '%s.sorted.merged.md.filtered.readcountall.out' % (sampleid)
            readcountall_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s\n' % (java, tmpdir, gatk, reference, infilteredbam)
            fo2.write(readcountall_task)
            #count on target reads
            #readcountontar = '%s.sorted.merged.md.filtered.readcountOnTar.out' % (sampleid)
            readcountontar_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s\n' % (java, tmpdir, gatk, reference, infilteredbam, target)
            fo2.write(readcountontar_task)
            #count flanked on target reads
            #readcountontar50 = '%s.sorted.merged.md.filtered.readcountOnTar.50bpflank.out' % (sampleid)
            readcountontar50_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s\n' % (java, tmpdir, gatk, reference, infilteredbam, targetflank50)
            fo2.write(readcountontar50_task)

            #readcountontar100 = '%s.sorted.merged.md.filtered.readcountOnTar.100bpflank.out' % (sampleid)
            readcountontar100_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s\n' % (java, tmpdir, gatk, reference, infilteredbam, targetflank100)
            fo2.write(readcountontar100_task)

    write_sub_script_qsub(cov_scriptdir, 'cov', cov_workdir, logdir)
    write_sub_script_qsub(cov_scriptdir, 'readcount', cov_workdir, cov_workdir)
   

#haplotypecaller
def write_sub_script_varcalling():
    #work with parent-sibling family data
    var_scriptdir = scriptdir + '/varcalling'
    var_workdir = workdir + '5_variant/0_rawvcf'
    mkdir(var_scriptdir)
    mkdir(var_workdir)
    
    familyID = SampleIDs[0].split('_', 1)[0]
    samplescript = var_scriptdir + '/' + familyID + '.haplotypecaller.task.pbs'
    
    with open(pedfile, 'rU') as fh0:
        for line in fh0:
            (familyID, sampleID, paternalID, maternalID, gender, phenotype) = line.strip().split()
            if paternalID != '0':
                paternalID = paternalID
                maternalID = maternalID
                siblingID = sampleID 

    #joint call for parents and siblings
    #ped file needed
    with open(samplescript, 'w+') as fo:
        fo.write('#!/bin/bash -l\n\n')
        fo.write('#PBS -l walltime=24:00:00\n')
        fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
        fo.write('#PBS -l mem=64gb\n')
        fo.write('#PBS -M huiwen.che@kuleuven.be\n')
        fo.write('#PBS -m abe\n')
        fo.write('#PBS -N varcalling_%s\n' % (familyID))
        fo.write('#PBS -A lp_biogenomics\n\n')
        fo.write('source switch_to_2015a\n')
        fo.write('module load SAMtools/1.7-foss-2015a\n')
        fo.write('module load GATK/3.7\n\n')
        
        align_workdir = workdir + '3_align'
        paternalbam = '%s/%s/%s.sorted.merged.md.filtered.bam' % (align_workdir, paternalID, paternalID) 
        maternalbam = '%s/%s/%s.sorted.merged.md.filtered.bam' % (align_workdir, maternalID, maternalID)
        siblingbam = '%s/%s/%s.sorted.merged.md.filtered.bam' % (align_workdir, siblingID, siblingID)
        jointrawgvcf = '%s/%s_parental_sib.joint.raw.snps.indels.vcf' % (var_workdir, familyID)
        gvcf_task = '%s -Xmx4G -jar %s -T HaplotypeCaller -R %s -A AlleleBalanceBySample -I %s -I %s -I %s --dbsnp %s -L %s -o %s\n\n' % (java, gatk, reference, paternalbam, maternalbam, siblingbam, dbsnpref, targetSNPref, jointrawgvcf)
        fo.write(gvcf_task)
        
        phasedrawvcf = '%s/%s_parental_sib.joint.raw.snps.indels.phasedbytrans.vcf' % (var_workdir, familyID)
        phasing_task = '%s -Xmx4G -jar %s -T PhaseByTransmission -R %s -V %s -ped %s --FatherAlleleFirst -o %s\n' % (java, gatk, reference, jointrawgvcf, pedfile, phasedrawvcf)
        fo.write(phasing_task)
        
    write_sub_script_qsub(var_scriptdir, 'haplotypecaller', var_workdir, logdir)



def write_sub_script_AE():
    AERC_scriptdir = scriptdir + '/AE'
    AERC_workdir = workdir + '6_AERC'
    mkdir(AERC_scriptdir)
    mkdir(AERC_workdir)
    
    #taking vcf generated from previous step 
    familyID = SampleIDs[0].split('_', 1)[0]
    var_workdir = workdir + '5_variant/0_rawvcf'
    refvcf = '%s/%s_parental_sib.joint.raw.snps.indels.vcf' % (var_workdir, familyID)
        
    for sampleid in SampleIDs:

        samplescript = AERC_scriptdir + '/' + sampleid + '.aerc.task.pbs'
        AERC_sub_workdir = AERC_workdir + '/' + sampleid
        mkdir(AERC_sub_workdir)
        
        with open(samplescript, 'w+') as fo:
            fo.write('#!/bin/bash -l\n\n')
            fo.write('#PBS -l walltime=03:00:00\n')
            fo.write('#PBS -l nodes=1:ppn=10:ivybridge\n')
            fo.write('#PBS -l mem=64gb\n')
            fo.write('#PBS -M huiwen.che@kuleuven.be\n')
            fo.write('#PBS -m abe\n')
            fo.write('#PBS -N AERC_%s\n' % (sampleid))
            fo.write('#PBS -A lp_biogenomics\n\n')
            fo.write('source switch_to_2015a\n')
            fo.write('module load SAMtools/1.7-foss-2015a\n')
            fo.write('module load GATK/3.7\n\n')
      
            align_workdir = workdir + '3_align'
            infilteredbam = '%s/%s/%s.sorted.merged.md.filtered.bam' % (align_workdir, sampleid, sampleid)
            mindepth = 2 #specify depth here, when fetal is male, to detect Y chromosome mutant the number of reads is on average half of the FF
            AEfrag = '%s.sorted.merged.md.filtered.bam.frag.aerc.%sx.csv' % (sampleid, mindepth)
            AE_task = '%s -Xmx4G -jar %s -T ASEReadCounter -R %s -I %s -sites %s -o %s/%s --countOverlapReadsType COUNT_FRAGMENTS_REQUIRE_SAME_BASE -U ALLOW_N_CIGAR_READS -minDepth %s --minMappingQuality 20 --minBaseQuality 2\n' % (java, gatk, reference, infilteredbam, refvcf, AERC_sub_workdir, AEfrag, mindepth)
            fo.write(AE_task)

    write_sub_script_qsub(AERC_scriptdir, 'aerc', AERC_workdir, logdir)


def write_scripts():
    write_sub_script_fastqc()
    write_sub_script_rename()
    write_sub_script_align()
    write_sub_script_coverage()
    write_sub_script_varcalling()
    write_sub_script_AE()
