#!/bin/bash

module load bcftools
module load htslib
#GATK_ver=gatk4300
GATK_ver=gatk4200

timeID=$(timeID)

#$1: VCFList
#$2: TagList (Normal_Ctrls -> Samples)

VCFList=$1
TagList=$2

cat ${TagList} | xargs -i File_pickByTag_onlyExitStatus.sh ${VCFList} {} || exit 1

for i in $(cat ${TagList})
do
  TargetVCF=$(grep -e ${i} ${VCFList})
  cp ${TargetVCF} ./${i}.merged.unfiltered.vcf
done

#cp ../*/*.merged.unfiltered.vcf ./. # VCFList can be replaced for this part.

ls *.vcf | xargs -i BsubS -n 4 -M 8 "bgzip -@ 4 {} && bcftools index {}.gz" BGZIP_INDEX >> JobID_BGZIP_${timeID}

sleep 300

WaitSignal=$(WaitSignalMaker.sh JobID_BGZIP_${timeID})

# VCF_gzs=$(ls *.gz | xargs -i printf {}" ")
VCF_gzs=$(cat ${TagList} | xargs -i printf {}.merged.unfiltered.vcf.gz" ")
BsubL -n 12 -M 60 "bcftools merge --merge all --threads 24 --force-samples ${VCF_gzs} > SampleMerged_VCFs.vcf && touch ${timeID}_checkpoint.txt" BcfMerge ${WaitSignal}

##########################################################################
#The codes below can be used for the MergeVCF_ChrSep of 1st Mu2Call_VCFs.#
##########################################################################

while ! [[ -e ${timeID}_checkpoint.txt ]]
do
  sleep 60
done

rm ${timeID}_checkpoint.txt

#After SampleMergedVCFs, this was split by Chr.
#Code example
for i in $(seq 1 22) X Y
do
  CMD_ChrSep='cat SampleMerged_VCFs.vcf | grep -E ''"'"^${i}\s|^#"'"'' > '"Chr${i}_merged.vcf"
  CMD_FileIndex="${GATK_ver} IndexFeatureFile --input Chr${i}_merged.vcf"
  BsubS "${CMD_ChrSep} && ${CMD_FileIndex}"
done

### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^1\s|^#" > Chr1_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr1_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^2\s|^#" > Chr2_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr2_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^3\s|^#" > Chr3_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr3_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^4\s|^#" > Chr4_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr4_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^5\s|^#" > Chr5_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr5_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^6\s|^#" > Chr6_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr6_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^7\s|^#" > Chr7_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr7_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^8\s|^#" > Chr8_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr8_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^9\s|^#" > Chr9_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr9_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^10\s|^#" > Chr10_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr10_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^11\s|^#" > Chr11_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr11_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^11\s|^#" > Chr12_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr12_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^13\s|^#" > Chr13_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr13_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^14\s|^#" > Chr14_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr14_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^15\s|^#" > Chr15_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr15_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^16\s|^#" > Chr16_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr16_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^17\s|^#" > Chr17_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr17_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^18\s|^#" > Chr18_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr18_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^19\s|^#" > Chr19_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr19_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^20\s|^#" > Chr20_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr20_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^21\s|^#" > Chr21_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr21_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^22\s|^#" > Chr22_merged.vcf && ${GATK_ver} IndexFeatureFile --input Chr22_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^X\s|^#" > ChrX_merged.vcf && ${GATK_ver} IndexFeatureFile --input ChrX_merged.vcf'
### BsubS 'cat SampleMerged_VCFs.vcf | grep -E "^Y\s|^#" > ChrY_merged.vcf && ${GATK_ver} IndexFeatureFile --input ChrY_merged.vcf'

exit

#########################################################################
# Following are to call Mu2Call of forced call with the given VCF loci. #
#########################################################################

# Index is required to use as the interval of Mutect2 Call.
# for i in $(ls *.vcf); do   BsubS -n 4 -M 8 "${GATK_ver} IndexFeatureFile --input ${i}"; done

# based on this SampleMerged_ChrSep.vcf, try second Mu2Call -> (-L / --allele option for specific locus, ) 
for i in $(cat TagList)
do
  BsubS "ForceCallonMergedVCF_SepChr_MedQue.sh ${i} BamList $(cat TargetMergedVCFDir)"
done
