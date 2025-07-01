#!/bin/bash

usage() {
  echo "Usage: $0 [-m] <VCFList> <TagList>"
  exit 1
}

reffa="hg19"

while getopts mh OPT
do
  case $OPT in
    m) reffa="mm10" ;;
    h) usage ;;
    *) usage ;;
  esac
done

shift $((OPTIND - 1))

# $1: BamList
# $2: TagList

BamList=$1
TagList=$2

mkdir -p FilMtx_From2ndMu2Call && cd FilMtx_From2ndMu2Call
mkdir -p UnfilteredVCF_on_1stMu2Call Mu2Call_ForceCallmVCF FilMtx_Create

cp ../${BamList} BamList
cp ../${TagList} TagList
ls ../Mu2Call_ToFMC_NoGNOMAD_gatk4200/*/*.merged.unfiltered.vcf | xargs realpath > UnfilteredVCF_on_1stMu2Call/UnfilteredVCFList_1stMu2Call
#ls ../Mu2Call_ToFMC_NoGNOMAD_gatk4300/*/*.merged.unfiltered.vcf | xargs realpath > UnfilteredVCF_on_1stMu2Call/UnfilteredVCFList_1stMu2Call

cp TagList UnfilteredVCF_on_1stMu2Call/.
cp BamList UnfilteredVCF_on_1stMu2Call/.
cp TagList FilMtx_Create/.

cd UnfilteredVCF_on_1stMu2Call/

if [ ${reffa} = "hg19" ] ; then
  echo "reffa = ${reffa}"
  BsubL -n 12 -M 60 "VCF_SampleMerger_toChrSep.sh -H UnfilteredVCFList_1stMu2Call TagList" VCFSampleMerge
elif [ ${reffa} = "mm10" ] ; then
  echo "reffa = ${reffa}"
  BsubL -n 12 -M 60 "VCF_SampleMerger_toChrSep.sh -m UnfilteredVCFList_1stMu2Call TagList" VCFSampleMerge
else
  echo "Species for reference is not defined."
fi
