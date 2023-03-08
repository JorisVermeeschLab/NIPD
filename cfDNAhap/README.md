# cfDNAhaplotyping
Genome-wide haplotyping of cell-free DNA

The scripts are used for non-invasive prenatal detection of inherited monogenic diseases and aneuploidy detection from cfDNA. *This is not the release of a software package. Scripts are only provided to reproduce analyses and they will be contiuously developed. The scripts are not for commercial use. For commercial use a license is required, please contact the lab.*


## Software

Python 2.7

R

## Sequencing data processing

Sequences are aligned and duplicates are marked with UMI barcodes.HaplotypeCaller was used to call variants from family gDNA samples jointly and parental genotypes were phased by PhaseByTransmission. Maternal plasma samples were handled separately, and allele counts were collected using ASEReadCounter by counting pair-end fragments requiring the overlapping bases to be identical, minimal mapping quality greater than 20, and base quality greater than 2.


## Phase parental genotype

To phase parental genotype, a vcf that contains variants of the family and a pedgree is needed. Lower and higher threshold of coverage for the vairants to be specified. If phasing with a proband from the family, gender of the proband should be provided.

**usage**: python ./main/parentalvcf_phasing.py <joint.phased.vcf> <family.ped> <lower_coverage_threshold> <higher_coverage_threshold> <workdir> <proband_gender>

## Estimate fetal fraction from cfDNA

Input data for estimating fetal fraction from cfDNA are type1 and type2 SNPs identified from parental genotype phasing step. Allele counts of cfDNA generated from ASEReadCounter is required. Lower and higher coverage thresholds, and median coverage of the cfDNA sample need to be specified.

**usage**: python ./main/ffest.py <type1_SNP> <type2_SNP> <cfDNA_aerc> <workdir> <lower_coverage_threshold> <higher_coverage_threshold> <median_coverage_threshold>

## Infer cffDNA haplotypes

To perform cffDNA haplotying, allele counts of cfDNA generated from ASEReadCounter is required. Lower and higher coverage thresholds, and median coverage of the cfDNA sample need to be specified. A correction factor to adjust reference allele bias can be specified. Fetal fraction and fetal gender information infered from previous step are used in this step as input. 

**usage**: python ./main/cffhap.CBS.py <cfDNA_aerc> <prefix_name> <family_name> <lower_coverage_threshold> <higher_coverage_threshold> <median_coverage_threshold> <correction_factor> <workdir> <scriptdir> <fetal_fraction> <fetal_gender>

## Family member SNP array data conversion

If parental and family member cyto12 SNParray data is used, first convert array data coding (0/1) to Illumina TOP/BOT A/C/T/G coding using constructed conversion file. Second, convert the tsv output to vcf format.

**usage**: python ./SNParray/fetchinfo.py ./SNParray/cytocombine.v2.tsv <family_SNParray_data> <converted_SNPcoding_to_ACTG_coding_output.tsv>

**usage**: python ./SNParray/array2vcf.py <converted_SNPcoding_to_ACTG_coding_output.tsv> <family_data.vcf>