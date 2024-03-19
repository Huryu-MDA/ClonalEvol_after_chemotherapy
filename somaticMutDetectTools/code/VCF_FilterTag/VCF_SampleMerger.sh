#!/bin/bash

module load bcftools
module load htslib

TimeID=$(timeID)

#$1: VCF_List
#$2: TagList (Normal_Ctrls -> Samples)

VCF_List=$1
TagList=$2

#VCFsort() {
#  Tag=$1
#  VCF=$2
#  #echo ${Tag}" "${VCF}
#  cat ${VCF} | grep "^#" > ${Tag}.merged.unfiltered.vcf
#  cat ${VCF} | grep -v "^#" | sort -k1,1n -k2,2n | grep -v "^X" | grep -v "^Y" >> ${Tag}.merged.unfiltered.vcf
#  cat ${VCF} | grep -v "^#" | sort -k1,1n -k2,2n | grep "^X" >> ${Tag}.merged.unfiltered.vcf
#  cat ${VCF} | grep -v "^#" | sort -k1,1n -k2,2n | grep "^Y" >> ${Tag}.merged.unfiltered.vcf
#}

for i in $(cat ${TagList} | grep -v "^#")
do
  if [[ $(cat ${VCF_List} | grep -e ${i} | wc -l) -ne 1 ]]
  then
    echo ${i}
    echo "Some of the VCF file cannot identified uniquely by Tag."
    exit
  fi
done


# Parallel_prepare
timeID_temp=$(timeID)
for i in $(cat ${TagList} | grep -v "^#")
do
  VCF=$(grep ${i} ${VCF_List})
  #echo ${i}
  echo "VCFsort ${i}.merged.unfiltered ${VCF}" >> CMD_list_temp${timeID_temp}
  #echo "echo ${i}.merged.unfiltered ${VCF}" >> CMD_list_${timeID_temp}
done

# Parallel_done
#cat CMD_list_temp${timeID_temp} | xargs -P 12 -i bash -c {}

# For the following line, the JOB_NAME on bjobs will be "*${timeID_temp:5:9}"
cat CMD_list_temp${timeID_temp} | xargs -i BsubM_jobid_task -n 4 -M 40 ${timeID_temp} "{}" VCFsort

sleep 30
JobnameOnBJOBS=${timeID_temp:5:9}
BJOB_COUNT=$(bjobs | grep "${JobnameOnBJOBS}" | wc -l)
echo "VCF_sort started" > BJOB_COUNT_${timeID_temp}.txt
echo -e "${JobnameOnBJOBS}\t${BJOB_COUNT}" >> BJOB_COUNT_${timeID_temp}.txt

# Loop to wait for the completion of previous BJOBs
while [ ${BJOB_COUNT} -ne 0 ] 
do
  sleep 60
  BJOB_COUNT=$(bjobs | grep "${JobnameOnBJOBS}" | wc -l)
  echo -e "${JobnameOnBJOBS}\t${BJOB_COUNT}" >> BJOB_COUNT_${timeID_temp}.txt
done

echo "VCF_sort done!" >> BJOB_COUNT_${timeID_temp}.txt

# CommandList_delete
rm CMD_list_temp${timeID_temp}

#exit
#cp ../*/*.merged.unfiltered.vcf ./. # VCFList can be replaced for this part.

ls *.vcf | xargs -P 8 -i bgzip -@ 4 {}
ls *.vcf.gz | xargs -P 8 -i bcftools index {}


#for i in $(ls *.vcf)
#do
#  BsubS -n 4 -M 8 "bgzip -@ 4 ${i} && bcftools index ${i}.gz" BGZIP_INDEX >> JobID_BGZIP_${TimeID}
#done

#sleep 300

#WaitSignal=$(WaitSignalMaker.sh JobID_BGZIP_${TimeID})

# VCF_gzs=$(ls *.gz | xargs -i printf {}" ")
VCF_gzs=$(cat ${TagList} | xargs -i printf {}.merged.unfiltered.vcf.gz" ")
#BsubS -n24 -M 60 "bcftools merge --merge all --threads 24 --force-samples ${VCF_gzs} > SampleMerged_VCFs.vcf && touch ${TimeID}_checkpoint.txt" BcfMerge ${WaitSignal}
bcftools merge --merge all --threads 24 --force-samples ${VCF_gzs} > SampleMerged_VCFs.vcf

mkdir -p VCFgz
mv *vcf.gz* VCFgz/.

exit 0
