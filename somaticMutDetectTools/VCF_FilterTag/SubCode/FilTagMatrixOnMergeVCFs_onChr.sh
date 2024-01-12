#!/bin/bash

module load python
module load bcftools

timeID=$(timeID)
output_filter_mtx="MTX_ChrToMu2F_INFO_VarType_GermLfil_VafFil_CovFil__${timeID}.tsv"

while getopts p: OPT
do
  case $OPT in
    p) output_filter_mtx="$OPTARG"_MTX_for_VCF_filter.tsv ;;
  esac
done

shift "$((OPTIND - 1))"

TargetMergedVCFs=$1
Spe=$2
Gender=$3

ref_fa_hg19="/rsrch3/scratch/reflib/iacs/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta"
ref_fa_mm10="/rsrch3/home/leuk-rsrch/huryu/Database/mm10/RefFa/mm10/fastas/mm10.fa"

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

bcftools norm -m- ${TargetMergedVCFs} | bcftools norm -f ${ref_fa} > Biallelic_${TargetMergedVCFs}
# bcftools query -f '%CHROM \t%POS \t%REF \t%ALT \t[ %AD\t]\n' Biallelic_${TargetMergedVCFs} | NullReadtoZeroReplace.py  > VarCovMtx_BiAlMergeVCF__${timeID}.tsv <- Older bcftools version??? Wasteful " " or "\t" seems to be the cause of error in the future codes.
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' Biallelic_${TargetMergedVCFs} | NullReadtoZeroReplace.py  > VarCovMtx_BiAlMergeVCF__${timeID}.tsv <- This will be used in the future. The below code, tee part can be removed.

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' Biallelic_${TargetMergedVCFs} | rev | cut -f 2- | rev | tee VarCovMtx_BiAlMergeVCF__${timeID}_NULL_UNREPLACED.tsv | NullReadtoZeroReplace.py  > VarCovMtx_BiAlMergeVCF__${timeID}.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' Biallelic_${TargetMergedVCFs} | rev | cut -f 2- | rev > VarCovMtx_BiAlMergeVCF__${timeID}_NULL_UNREPLACED.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' Biallelic_${TargetMergedVCFs} | NullReadtoZeroReplace.py | rev | cut -f 2- | rev > VarCovMtx_BiAlMergeVCF__${timeID}.tsv


GermlF_map.py VarCovMtx_BiAlMergeVCF__${timeID}.tsv > VarCovMtx_Germl_Somatic_FilTag__${timeID}.tsv
cat VarCovMtx_Germl_Somatic_FilTag__${timeID}.tsv | grep -v "Somatic mutation candidate" > VarCovMtx_GermlFil__${timeID}.tsv
cat Biallelic_${TargetMergedVCFs} | grep -v ^## | cut -f 1-8 > VCF_body_ChrToMu2F_INFO__${timeID}.tsv
cat VCF_body_ChrToMu2F_INFO__${timeID}.tsv | VarType_Tag.py | RefMtx_Tag.py VarCovMtx_GermlFil__${timeID}.tsv | VafF_CovF_Tag.py -lc 5 -mc 5 -g ${Gender} VarCovMtx_BiAlMergeVCF__${timeID}.tsv > ${output_filter_mtx}
