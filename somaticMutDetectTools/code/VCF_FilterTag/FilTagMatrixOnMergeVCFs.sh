#!/bin/bash

module load python
module load bcftools

TargetMergedVCFs=$1
Spe=$2
Gender=$3
timeID=$(timeID)

ref_fa_hg19="/rsrch3/scratch/reflib/iacs/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta"
ref_fa_mm10="/rsrch3/home/leuk-rsrch/huryu/Database/mm10/RefFa/mm10/fastas/mm10.fa"
CH_List="1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\nX\nY"

if [[ ${Spe} = "hg19" ]]
then
  ref_fa=${ref_fa_hg19}
elif [[ ${Spe} = "mm10" ]]
then
  ref_fa=${ref_fa_mm10}
else
  echo "Please check the arg1, arg2, arg3 for this command (FilTagMatrixOnMergeVCFs.sh)"
  exit
fi

# base_vcf_matrix_separation
cp ${TargetMergedVCFs} TargetMergedVCFs.vcf
echo -e $CH_List | xargs -i mkdir -p chr{}
echo -e $CH_List | xargs -i -P 24 bash -c 'grep -E "^#|^{}\s" TargetMergedVCFs.vcf > chr{}/chr{}_TargetMergedVCFs.vcf'

# remove the files derived from the previous execution. 
rm -f TargetMergedVCFs.vcf  chr*/OutLog  chr*/ErrLog  JobID_chrFilTagMatrix

# Filter matrix on each chromosome
for CH in $(echo -e ${CH_List})
do
  #cd chr${CH} && BsubS "FilTagMatrixOnMergeVCFs_onChr.sh TargetMergedVCFs.vcf ${Spe} ${Gender}" && cd ..
  cd chr${CH} && BsubS -n 6 -M 60 "FilTagMatrixOnMergeVCFs_onChr.sh -p chr${CH} chr${CH}_TargetMergedVCFs.vcf ${Spe} ${Gender}" >> ../JobID_chrFilTagMatrix && cd ..
done

WaitSignal=$(WaitSignalMaker.sh JobID_chrFilTagMatrix)

BsubS -n 6 -M 30 "__filMTXonChrs_merge__.sh" filMTXonChrs_merging ${WaitSignal} > JobID_filMTXonChrs_merging

exit 0
