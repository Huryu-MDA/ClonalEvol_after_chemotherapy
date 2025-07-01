#!/bin/bash

# $1: TagList
# $2: Spe ~ hg19, mm10, etc
# $3: Gender

if ! [[ $# -eq 3 ]]
then
  echo -e "Example of usage: \nVCF_concat_merge_to_FilMTX.sh TagList_816968_2010 hg19 F"
  exit 1
fi

TagList=$1
Spe=$2
Gender=$3

mkdir -p FilMtx_From2ndMu2Call && cd FilMtx_From2ndMu2Call
echo "Now in FilMtx_From2ndMu2Call ..."

if [[ -e JobID_MergeVCFonEachSamples ]]
then
  rm JobID_MergeVCFonEachSamples
fi

cat ../${TagList} | xargs -i rm -f Mu2Call_ForceCallmVCF/{}/{}.ChrMerged.unfiltered.vcf
cat ../${TagList} | xargs -i rm -f Mu2Call_ForceCallmVCF/{}/{}.ChrMerged.unfiltered.vcf.idx
rm -f MergeVCF_of_EachSamples_???.log

#cat CodesList_toRunAfterMu2Call_Complete | grep "MergeVCF" | xargs -i BsubM -n 3 -M 8 "module load openjdk/11.0.5-10 && {}" MergeVCF_of_EachSamples >> JobID_MergeVCFonEachSamples
cat CodesList_toRunAfterMu2Call_Complete | grep "MergeVCF" | xargs -i BsubM -n 3 -M 8 "module load openjdk/11.0.15+10 && {}" MergeVCF_of_EachSamples >> JobID_MergeVCFonEachSamples

WaitSignal=$(WaitSignalMaker.sh JobID_MergeVCFonEachSamples)

mkdir -p FilMtx_Create && cd FilMtx_Create/
echo "Now in FilMtx_From2ndMu2Call/FilMtx_Create ..."
#echo ${TagList}
cp ../../${TagList} ./TagList


echo "Now waiting Mu2_2ndcallVCF_merging ..."
Tag1=$(head -n 1 TagList)
while ! [[ -e ../Mu2Call_ForceCallmVCF/${Tag1}/${Tag1}.ChrMerged.unfiltered.vcf ]]
do
  sleep 60
done


#while ! [[ $(ls ../Mu2Call_ForceCallmVCF/*/*.ChrMerged.unfiltered.vcf | wc -l) -eq $(cat TagList | wc -l) ]]
#while ! [[ $(cat TagList | xargs -i ls ../Mu2Call_ForceCallmVCF/{}/{}.ChrMerged.unfiltered.vcf | wc -l) -eq $(cat TagList | wc -l) ]]
while ! [[ $(ls ../Mu2Call_ForceCallmVCF/*/*.ChrMerged.unfiltered.vcf | wc -l) -eq $(cat TagList | wc -l) ]]
do
  echo "Merged_VCF_from_Forced_Mu2_Call: "$(ls ../Mu2Call_ForceCallmVCF/*/*.ChrMerged.unfiltered.vcf | wc -l)
  echo "Expected Merged_VCF_from_Forced_Mu2_Call: "$(cat TagList | wc -l)
  echo "Waiting ..."
  sleep 60
done

ls ../Mu2Call_ForceCallmVCF/*/*.ChrMerged.unfiltered.vcf | xargs realpath > VCFList

CMD_VCFmerge="VCF_SampleMerger.sh VCFList TagList"
CMD_FilMTX="FilTagMatrixOnMergeVCFs.sh SampleMerged_VCFs.vcf ${Spe} ${Gender}"
BsubM -n 10 -M 20 "${CMD_VCFmerge} && ${CMD_FilMTX}" VCFmerge_FilMTX ${WaitSignal} > JobID_VCFmerge_FilMTX
