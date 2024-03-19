#!/bin/bash

module -q load python

# $1: VCFList
# $2: Tag
# $3: FilMTX_List (FilterMatrix_onChr_List)
# $4: Gender "M" or "F"

vcf_list=$1
Tag=$2
FilMTX_List=$3
Gender=$4

File_pickByTag_onlyExitStatus.sh ${vcf_list} ${Tag} || exit 1
vcf=$(File_pickByTag.sh ${vcf_list} ${Tag})

CH_List="1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\nX\nY"
echo -e $CH_List | xargs -i mkdir -p ${Tag}/chr{}
echo -e $CH_List | xargs -i -P 24 bash -c "grep -E '^#|^{}\s' ${vcf} > ${Tag}/chr{}/${Tag}_chr{}.vcf"

rm -f ${Tag}/chr*/OutLog  ${Tag}/chr*/ErrLog

for i in $(echo -e ${CH_List})
do
  FilMTX_CH=$(grep chr${i}_ ${FilMTX_List})
  FilMTX_CH=$(realpath ${FilMTX_CH})
  cd ${Tag}/chr${i} && \
  BsubS -n 1 -M 8 "cat ${Tag}_chr${i}.vcf | VarType_Tag.py | dist10BP_Tag.py | FilTag_ConcatToVCFBody.py ${FilMTX_CH} | LowCovTagOnSingleVCF.py -g ${Gender} -m 6 > 7FilTagged_${Tag}_chr${i}.merged.filtered.vcf" >> ../JobID_FilTagVCF && \
  cd ../..
done

cd ${Tag}
WaitSignal=$(WaitSignalMaker.sh  JobID_FilTagVCF)

#<concatenate separate_VCF>
BsubM -n 6 -M 40 "__postfilTaggedVCFonChrs_merge__.sh  ${Tag}"  postfilTaggedVCFonChrs_merging  ${WaitSignal} > JobID_postfilTaggedVCFonChrs_merging

waitsignal_for_file_compressing=$(WaitSignalMaker.sh  JobID_postfilTaggedVCFonChrs_merging)
sleep 60

# <Folder compressed and stored>
cd ..
BsubS -n 2 -M 30 "__taggedVCF_splitByChr_compress__.sh ${Tag}" split_vcf_dir_compression ${waitsignal_for_file_compressing} > JobID_split_VCF_on_chrs_compression

exit 0
