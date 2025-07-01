#!/bin/bash

# $1: TagList, e.g.TagList_8269678_2010
# $2: Gender, F or M


FLAG_mm10="FALSE"
while getopts m OPT
do
  case $OPT in
    m) FLAG_mm10="TRUE" ;;
  esac
done

shift $((OPTIND - 1))

TagList=$1
Gender=$2

mkdir -p TagVCF_withFilMTX && cp ${TagList} TagVCF_withFilMTX/TagList && cd TagVCF_withFilMTX/
#ls ../Mu2Call_ToFMC_NoGNOMAD_gatk4200/*/*.merged.filtered.vcf | xargs -i realpath {} > Mu2Fil_VCFList
ls ../Mu2Call_ToFMC_NoGNOMAD_gatk4200/*/*.merged.filtered.vcf | xargs -i realpath {} > Mu2Fil_VCFList
ls ../FilMtx_From2ndMu2Call/FilMtx_Create/chr*/chr*_MTX_*.tsv > FilMTX_onChr_list.txt



if [[ ${FLAG_mm10} = "FALSE" ]]; then
  cat TagList | xargs -i BsubM -n 8 -M 40 "FilTaggingTo_Mu2filteredVCF.sh Mu2Fil_VCFList {} FilMTX_onChr_list.txt ${Gender}" FilTagging >> JobID_FilTagToVCF
else
  cat TagList | xargs -i BsubM -n 8 -M 40 "FilTaggingTo_Mu2filteredVCF.sh -m Mu2Fil_VCFList {} FilMTX_onChr_list.txt ${Gender}" FilTagging >> JobID_FilTagToVCF
fi
