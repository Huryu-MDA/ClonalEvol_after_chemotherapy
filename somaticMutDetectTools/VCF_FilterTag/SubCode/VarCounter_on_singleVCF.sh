#!/bin/bash

module load python
#$1: vcf_list
#$2: Tag

vcf_list=$1
Tag=$2

File_pickByTag_onlyExitStatus.sh ${vcf_list} ${Tag} || exit 1 
vcf=$(File_pickByTag.sh ${vcf_list} ${Tag})
vcf=$(realpath ${vcf})
vcf_output=7Fil_passed_$(basename ${vcf})

# Make the VCF that rows have passed all filters.
cat ${vcf} | grep "^#" > ${vcf_output}
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" >> ${vcf_output}

# SNV part

echo ${Tag} >> SNV_SampleNameCol.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> SNV_count_7Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> SNV_count_6Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> SNV_count_5Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> SNV_count_4Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> SNV_count_3Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | awk_Col.py 12 "None" | wc -l >> SNV_count_2Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "SNV" | wc -l >> SNV_count_1Fil.tsv
cat ${vcf} | awk_Col.py 11 "SNV" | wc -l >> SNV_count_0Fil.tsv

echo ${Tag} >> MNV_SampleNameCol.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> MNV_count_7Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> MNV_count_6Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> MNV_count_5Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> MNV_count_4Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> MNV_count_3Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | awk_Col.py 12 "None" | wc -l >> MNV_count_2Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "MNV" | wc -l >> MNV_count_1Fil.tsv
cat ${vcf} | awk_Col.py 11 "MNV" | wc -l >> MNV_count_0Fil.tsv

echo ${Tag} >> INDEL_SampleNameCol.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | awk_Col.py 17 "None" | wc -l >> INDEL_count_7Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | awk_Col.py 16 "None" | wc -l >> INDEL_count_6Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | awk_Col.py 15 "None" | wc -l >> INDEL_count_5Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | awk_Col.py 14 "None" | wc -l >> INDEL_count_4Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | awk_Col.py 13 "None" | wc -l >> INDEL_count_3Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | awk_Col.py 12 "None" | wc -l >> INDEL_count_2Fil.tsv
cat ${vcf} | awk_Col.py 7 "PASS" | awk_Col.py 11 "INDEL" | wc -l >> INDEL_count_1Fil.tsv
cat ${vcf} | awk_Col.py 11 "INDEL" | wc -l >> INDEL_count_0Fil.tsv

paste SNV_SampleNameCol.tsv SNV_count_0Fil.tsv SNV_count_1Fil.tsv SNV_count_2Fil.tsv SNV_count_3Fil.tsv SNV_count_4Fil.tsv SNV_count_5Fil.tsv SNV_count_6Fil.tsv SNV_count_7Fil.tsv > VarCounts_SNV_onFil_1to7.tsv
paste MNV_SampleNameCol.tsv MNV_count_0Fil.tsv MNV_count_1Fil.tsv MNV_count_2Fil.tsv MNV_count_3Fil.tsv MNV_count_4Fil.tsv MNV_count_5Fil.tsv MNV_count_6Fil.tsv MNV_count_7Fil.tsv > VarCounts_MNV_onFil_1to7.tsv
paste INDEL_SampleNameCol.tsv INDEL_count_0Fil.tsv INDEL_count_1Fil.tsv INDEL_count_2Fil.tsv INDEL_count_3Fil.tsv INDEL_count_4Fil.tsv INDEL_count_5Fil.tsv INDEL_count_6Fil.tsv INDEL_count_7Fil.tsv > VarCounts_INDEL_onFil_1to7.tsv

module load python/3.8.8-anaconda
#SNV_counts_To_LineBoxPlot.py VarCounts_SNV_onFil_1to7.tsv
SNV_counts_To_LineBoxBarPlot.py VarCounts_SNV_onFil_1to7.tsv

exit
