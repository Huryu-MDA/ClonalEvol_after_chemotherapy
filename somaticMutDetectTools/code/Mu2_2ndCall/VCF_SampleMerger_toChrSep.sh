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

