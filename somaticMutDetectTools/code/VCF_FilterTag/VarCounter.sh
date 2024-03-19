#!/bin/bash

module load python
#$1: VCFList
#$2: TagList

VCFList=$1
TagList=$2

# FilTag_Remed_VCF

cat ${TagList} | grep -v "^#" | xargs -i File_pickByTag_onlyExitStatus.sh ${VCFList} {} || exit 1

for i in $(cat ${TagList})
do
  TargetVCF=$(cat ${VCFList} | grep ${i})
  SaveVCF=7Fil_passed_${TargetVCF}
  cat ${TargetVCF} | grep "^#" > ${SaveVCF}
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" >> ${SaveVCF}
done

#exit

# SNV part

for i in $(cat ${TagList})
do
  TargetVCF=$(cat ${VCFList} | grep ${i})
  #echo ${TargetVCF}
  echo ${i} >> SNV_SampleNameCol.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> SNV_count_7Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> SNV_count_6Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> SNV_count_5Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> SNV_count_4Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> SNV_count_3Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | wc -l >> SNV_count_2Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | wc -l >> SNV_count_1Fil.tsv
  cat ${TargetVCF} | awk_Col.py 11 "SNV" | wc -l >> SNV_count_0Fil.tsv
done

paste SNV_SampleNameCol.tsv SNV_count_0Fil.tsv SNV_count_1Fil.tsv SNV_count_2Fil.tsv SNV_count_3Fil.tsv SNV_count_4Fil.tsv SNV_count_5Fil.tsv SNV_count_6Fil.tsv SNV_count_7Fil.tsv > VarCounts_SNV_onFil_1to7.tsv

# MVN part

for i in $(cat ${TagList})
do
  TargetVCF=$(cat ${VCFList} | grep ${i})
  echo ${i} >> MNV_SampleNameCol.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> MNV_count_7Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> MNV_count_6Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> MNV_count_5Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> MNV_count_4Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> MNV_count_3Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | wc -l >> MNV_count_2Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | wc -l >> MNV_count_1Fil.tsv
  cat ${TargetVCF} | awk_Col.py 11 "MNV" | wc -l >> MNV_count_0Fil.tsv
done

paste MNV_SampleNameCol.tsv MNV_count_0Fil.tsv MNV_count_1Fil.tsv MNV_count_2Fil.tsv MNV_count_3Fil.tsv MNV_count_4Fil.tsv MNV_count_5Fil.tsv MNV_count_6Fil.tsv MNV_count_7Fil.tsv > VarCounts_MNV_onFil_1to7.tsv

# INDEL part

for i in $(cat ${TagList})
do
  TargetVCF=$(cat ${VCFList} | grep ${i})
  echo ${i} >> INDEL_SampleNameCol.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> INDEL_count_7Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> INDEL_count_6Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> INDEL_count_5Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> INDEL_count_4Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> INDEL_count_3Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | wc -l >> INDEL_count_2Fil.tsv
  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | wc -l >> INDEL_count_1Fil.tsv
  cat ${TargetVCF} | awk_Col.py 11 "INDEL" | wc -l >> INDEL_count_0Fil.tsv
done

paste INDEL_SampleNameCol.tsv INDEL_count_0Fil.tsv INDEL_count_1Fil.tsv INDEL_count_2Fil.tsv INDEL_count_3Fil.tsv INDEL_count_4Fil.tsv INDEL_count_5Fil.tsv INDEL_count_6Fil.tsv INDEL_count_7Fil.tsv > VarCounts_INDEL_onFil_1to7.tsv

module load python/3.8.8-anaconda
#SNV_counts_To_LineBoxPlot.py VarCounts_SNV_onFil_1to7.tsv
SNV_counts_To_LineBoxBarPlot.py VarCounts_SNV_onFil_1to7.tsv

FileIO_postVarCount.sh ${TagList}

exit

# FilTag_Remed_VCF

#for i in $(cat ${TagList})
#do
#  TargetVCF=$(cat ${VCFList} | grep ${i})
#  SaveVCF=6Fil_passed_${TargetVCF}
#  cat ${TargetVCF} | grep "^#" > ${SaveVCF}
#  cat ${TargetVCF} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" >> ${SaveVCF}
#done
