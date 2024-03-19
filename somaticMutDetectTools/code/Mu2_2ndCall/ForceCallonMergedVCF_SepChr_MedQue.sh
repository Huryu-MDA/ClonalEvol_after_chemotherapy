#!/bin/bash

module load openjdk/11.0.5-10
WorkDir=Mu2Call_ForceCallmVCF
# $1:Bam_Tag
# $2:BamList
# $3:TargetVCF_Interval(MergedSplitUnfilteredVCF)
# $4:Species hg19 / mm10

BamTag=$1
BamList=$2
TargetVCFDir=$3
Spe=$4
timeID=$5

RefFa_hg19="/rsrch3/scratch/reflib/iacs/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta"
RefFa_mm10="/rsrch3/home/leuk-rsrch/huryu/Database/mm10/RefFa/mm10/fastas/mm10.fa"

CH_List_hg19="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
CH_List_mm10="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y"

if [[ $# -lt 3 ]]
then
  echo "You need 3 args. #1 BamTag #2 BamList #3 folder pass that contains Chr_Split_Sample_Merged_unfiltered_VCF."
  exit
fi

File_pickByTag_onlyExitStatus.sh ${BamList} ${BamTag} || exit 1 
InputBam=$(File_pickByTag.sh ${BamList} ${BamTag})

#InputBam=$(cat ${BamList} | grep -e ${BamTag})
#InputBamNum=$(cat ${BamList} | grep -e ${BamTag} | wc -l)

#if [[ ${InputBamNum} -ne 1 ]]
#then
#  echo "Tag cannot identify unique bam file. Please check the tag / bam file list."
#  exit
#fi

# RefFa setting
if [[ $# -eq 3 ]]
then
  RefFa=${RefFa_hg19}
  CH_List=${CH_List_hg19}
fi

if [[ ${Spe} = "hg19" ]]
then
  RefFa=${RefFa_hg19}
  CH_List=${CH_List_hg19}
elif [[ ${Spe} = "mm10" ]]
then
  RefFa=${RefFa_mm10}
  CH_List=${CH_List_mm10}
fi

if [[ ${RefFa} = "" ]]
then
  echo "RefFa was not defined. Check the 4th Argument."
  exit
fi
# RefFa setting done


GATK=gatk4200

MemVol=24
JavaOptions='--java-options "-Xmx'"${MemVol}"'G"'


for i in ${CH_List}
do
  StoreDir=${WorkDir}/${BamTag}/Chr${i}/
  TargetVCF_Interval=${TargetVCFDir}/Chr${i}_merged.vcf
  mkdir -p ${StoreDir}
  CodeFile=${StoreDir}/../Mu2Call_${BamTag}_Chr${i}
  echo "#!/bin/bash" > ${CodeFile}
  echo ""            >> ${CodeFile}
  echo "module load openjdk/11.0.5-10" >> ${CodeFile}
  echo ""            >> ${CodeFile}
  echo "${GATK} ${JavaOptions} \\" >> ${CodeFile}
  echo "Mutect2 -R ${RefFa} \\" >> ${CodeFile}
  echo "-L ${i} \\" >> ${CodeFile}
  echo "--alleles ${TargetVCF_Interval} \\" >> ${CodeFile}
  echo "-I ${InputBam} \\" >> ${CodeFile}
  echo "-O ${StoreDir}/${BamTag}.Chr${i}.unfiltered.vcf" >> ${CodeFile}
  BsubM -n 2 -M ${MemVol} "bash ${CodeFile}" ${StoreDir}/../Mu2Call_Chr${i} >> ${StoreDir}/../JobID_Mu2Call_Chr${i}
  cat ${StoreDir}/../JobID_Mu2Call_Chr${i} >> JobID_Mu2ForceCall_onChrs_${timeID}
done

Input_UnFilVcfs=$(printf "${CH_List}" | xargs -d " " -i printf "-I ${WorkDir}/${BamTag}/Chr"{}"/${BamTag}.Chr"{}".unfiltered.vcf ")
MergeVCF_CMD="${GATK} MergeVcfs ${Input_UnFilVcfs} -O ${WorkDir}/${BamTag}/${BamTag}.ChrMerged.unfiltered.vcf"
echo "${MergeVCF_CMD}" > ${WorkDir}/${BamTag}/MergeVCF_CMD_${BamTag}

if ! [[ -e CodesList_toRunAfterMu2Call_Complete ]]
then
  echo "#!/bin/bash" > CodesList_toRunAfterMu2Call_Complete
  echo "" >> CodesList_toRunAfterMu2Call_Complete
  echo "module load openjdk/11.0.5-10" >> CodesList_toRunAfterMu2Call_Complete 
  echo "" >> CodesList_toRunAfterMu2Call_Complete
fi

echo "${WorkDir}/${BamTag}/MergeVCF_CMD_${BamTag}" >> CodesList_toRunAfterMu2Call_Complete
chmod +x ${WorkDir}/${BamTag}/MergeVCF_CMD_${BamTag}

exit
